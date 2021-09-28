function Dinv = CalculateDinv(f, dy, k0, kx, tlineup_, tlinedown_, z0, wslot, Ny, walled)
    kx(kx == 0) = (1+1j)*eps;
    
    %% Determine integration path
    % Go to +-inf
    delta_ky = 0.01.*k0;
    lim1_ky = -5.*k0-1j.*delta_ky;
    lim2_ky = -lim1_ky;
    % Deform integration path around branch cuts.
    integrationpath_ky = [(-1-1j).*delta_ky, (1+1j).*delta_ky];
    
    
    Dinv = zeros(length(kx), Ny, Ny);
%     Ds = zeros(length(kx), Ny);
%     Ddowns = zeros(length(kx), 1);
%     Dfss = zeros(length(kx), Ny);
    parfor(kxii = 1:length(kx)) % parfor
        if(kxii == round(length(kx)/2))
        end
        
        D = zeros(1, Ny);
        for(dny = 0:Ny-1)
            D(dny+1) = 1 ./ (2*pi) .* fastintegral(...
                @(kyp) D_Integrand_kyp(f, dy, k0, kyp, kx(kxii), tlineup_, tlinedown_, z0, wslot, dny, walled), ...
                lim1_ky, lim2_ky, 'Waypoints', integrationpath_ky);
        end
        K = -1i*sqrt(-(k0.^2-kx(kxii).^2));
        Dfs = [-0.5/k0/z0*(k0^2-kx(kxii).^2).*besselj(0,wslot/4*K)*besselh(0,2,wslot/4*K), ... % self
               -0.5/k0/z0*(k0^2-kx(kxii).^2).*besselh(0,2,K*dy*(1:Ny-1))];                     % mutual
        if(walled)
            % Only add it to up, since the down contribution is calculated
            % in the sum below.
            D = D + Dfs;
        else
            % Add it for both up and down.
            D = D + 2.*Dfs;
        end
        
        % If there's walls, calculate the contribution from the
        % stratification below.
        if(walled)
            Ddown = 0;
            for(my = -20:20)
                Vtm = 1;
                Vte = 1;

                kym = -2*pi*my/dy;
                kr = sqrt(kx(kxii).^2 + kym.^2);
                isTE = 1;    ztedown = tlinedown_.GetInputImpedance(isTE, f, k0, kr);
                isTE = 0;    ztmdown = tlinedown_.GetInputImpedance(isTE, f, k0, kr);

                itedown = 1 ./ ztedown;
                itmdown = 1 ./ ztmdown;

                [Ghm] = SpectralGF.hm(z0, k0, kx(kxii), kym, Vtm, Vte, itmdown, itedown);
                Ghm_xx_down = Ghm.xx;

                Ddown = Ddown + 1/dy .* Ghm_xx_down .* besselj(0, -2*pi*my/dy*wslot/2);
            end
        else
            Ddown = 0;
        end
%         Ddowns(kxii) = Ddown;
%         Ds(kxii, :) = D;
        % Add Ddown to the self.
        D(1) = D(1) + Ddown;
        
        Dmat = complextoeplitz(D);
        Dinv(kxii, :, :) = pinv(Dmat);
    end
    %{
    for(ny = 1:Ny)
    [~, ~] = figureex(101);
        hAx = subplot(Ny, 1, ny);
        hold(hAx, 'on'); grid(hAx, 'on'); box(hAx, 'on');
%         repeatcolormap(hAx, 2);
        plot(hAx, real(kx./k0), real(Ds(:,ny)), 'x');
        plot(hAx, real(kx./k0), imag(Ds(:,ny)), 'o')
    end
    %}
    %{
    [~, hAx] = figureex(100);
        plot(hAx, real(kx./k0), real(Ddowns), 'x');
        plot(hAx, real(kx./k0), imag(Ddowns), 'o')
    %}
    %{
    [~, hAx] = figureex;
        hAx.ColorOrder = lines(Ny);
        plot(hAx, real(kx./k0), real(Ds));
        plot(hAx, real(kx./k0), imag(Ds), '--')
        plot(hAx, real(kx./k0), real(Dfss), 'x');
        plot(hAx, real(kx./k0), imag(Dfss), 'o')
    %}
end
%{
function v = D_Integrand_kyp_alternate(f, dy, k0, kyp, kx, tlineup, tlinedown, z0, wslot, nypp)
    Vtm = 1;
    Vte = 1;
    
    tlineupfs = FreeSpace();
    tlinedownfs = FreeSpace();
    kr = sqrt(kx.^2 + kyp.^2);
    
    % Use the real Green's function inside this interval, otherwise use FreeSpace
    % The ratio of abs(kyp)/k0 varies for different error levels
    %     Values determined for an er=3.0 ADL
    %     Error | 1e-1  | 1e-2  | 1e-3  | 1e-4  | 1e-5  | 1e-6  | 1e-7  | 1e-8  | 1e-9  | 1e-10
    %     ratio |  1.12 |  1.90 |  4.57 |  7.80 | 10.80 | 13.25 | 14.80 | 22.81 | 26.93 | 30.26 
    ratio = 15;
    iRealGF = (abs(kyp) < ratio*k0);
    zteup = zeros(size(kr));    ztmup = zteup;    ztedown = zteup;    ztmdown = zteup;
    
    isTE = 1;    zteup(iRealGF) = tlineup.GetInputImpedance(isTE, f, k0, kr(iRealGF));
    isTE = 0;    ztmup(iRealGF) = tlineup.GetInputImpedance(isTE, f, k0, kr(iRealGF));
    isTE = 1;    zteup(~iRealGF) = tlineupfs.GetInputImpedance(isTE, f, k0, kr(~iRealGF));
    isTE = 0;    ztmup(~iRealGF) = tlineupfs.GetInputImpedance(isTE, f, k0, kr(~iRealGF));

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kx, kyp, Vtm, Vte, itmup, iteup);
    Gxxup = Ghm.xx;

    if(~isempty(tlinedown))
        isTE = 1;    ztedown(iRealGF) = tlinedown.GetInputImpedance(isTE, f, k0, kr(iRealGF));
        isTE = 0;    ztmdown(iRealGF) = tlinedown.GetInputImpedance(isTE, f, k0, kr(iRealGF));
        isTE = 1;    ztedown(~iRealGF) = tlinedownfs.GetInputImpedance(isTE, f, k0, kr(~iRealGF));
        isTE = 0;    ztmdown(~iRealGF) = tlinedownfs.GetInputImpedance(isTE, f, k0, kr(~iRealGF));

        itedown = 1 ./ ztedown;
        itmdown = 1 ./ ztmdown;

        [Ghm] = SpectralGF.hm(z0, k0, kx, kyp, Vtm, Vte, itmdown, itedown);
        Gxxdown = Ghm.xx;
    else
        Gxxdown = 0;
    end

    Gxx = Gxxup + Gxxdown;

    v = Gxx.*besselj(0, kyp.*wslot./2).*exp(-1j.*kyp.*nypp.*dy);
end
%}
function integrand = D_Integrand_kyp(f, dy, k0, kyp, kx, tlineup, tlinedown, z0, wslot, dny, walled)
    [Gxxup, Gxxdown] = GreensFunction(f, k0, kx, kyp, tlineup, tlinedown, z0, walled);
    vup = Gxxup.*besselj(0, kyp.*wslot./2).*exp(-1j.*kyp.*dny.*dy);
    vdown = Gxxdown.*besselj(0, kyp.*wslot./2).*exp(-1j.*kyp.*dny.*dy);
    
    % Free-space
    fsup = FreeSpace();
    fsdown = [];
    if(~isempty(tlinedown))
        fsdown = fsup;
    end
    [Gxxfsup, Gxxfsdown] = GreensFunction(f, k0, kx, kyp, fsup, fsdown, z0, walled);
    
    vfsup = Gxxfsup.*besselj(0, kyp.*wslot./2).*exp(-1j.*kyp.*dny.*dy);
    vfsdown = Gxxfsdown.*besselj(0, kyp.*wslot./2).*exp(-1j.*kyp.*dny.*dy);
    
    vup(abs(kyp) > 5*k0) = vfsup(abs(kyp) > 5*k0);
    
    integrand = vup + vdown;
    
    integrand = integrand - vfsup - vfsdown;
end
function [Ghm_xx_up, Ghm_xx_down] = GreensFunction(f, k0, kx, ky, tlineup, tlinedown, z0, walled)
    Vtm = 1;
    Vte = 1;
    
    kr = sqrt(kx.^2 + ky.^2);
    isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);

    iteup = Vte ./ zteup;
    itmup = Vtm ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kx, ky, Vtm, Vte, itmup, iteup);
    Ghm_xx_up = Ghm.xx;

    if(~walled)
        isTE = 1;    ztedown = tlinedown.GetInputImpedance(isTE, f, k0, kr);
        isTE = 0;    ztmdown = tlinedown.GetInputImpedance(isTE, f, k0, kr);

        itedown = Vte ./ ztedown;
        itmdown = Vtm ./ ztmdown;

        [Ghm] = SpectralGF.hm(z0, k0, kx, ky, Vtm, Vte, itmdown, itedown);
        Ghm_xx_down = Ghm.xx;
    else
        Ghm_xx_down = 0;
    end
end





















