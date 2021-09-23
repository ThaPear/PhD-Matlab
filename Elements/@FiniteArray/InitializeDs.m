function InitializeDs(this, fs)
    newfs = setdiff(fs, this.D_fs);
    if(length(newfs) < 1)
        return;
    end

    dispex('Ds: Calculating for %i frequencies, %i slots, 2 integration paths.\n', length(newfs), this.Ny);
    tc = tic;
    
    Ny_ = this.Ny;

    dy = this.unitcell.dy;
    wslot = this.unitcell.wslot;
    walled = this.unitcell.walled;
    tlineup_ = this.tlineup;
    tlinedown_ = this.tlinedown;


    z0 = Constants.z0;
    c0 = Constants.c0;

    Nf = length(newfs);

    deformedpath = 0;

    %% Initialize progress bar
    hDataQueue = parallel.pool.DataQueue;
    hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
    afterEach(hDataQueue, @updateWaitbar);
    progress = -1; send(hDataQueue, nan);

    % deformedpath: 0 means straight integration path, 1 means deformed integration path
    for(deformedpath = 0)
        Dinvlist = cell(1,Nf);
        kxlist = cell(1,Nf);
        for(fi = 1:Nf)
%             dispex('Calculating %i.\n', ni);
            f = newfs(fi); %#ok<PFBNS> % Broadcast variable.

            %% Determine integration path
            k0 = 2*pi*f / c0;
            % Go to +-inf
            delta_ky = 0.01.*k0;
            lim1_ky = -5.*k0-1j.*delta_ky;
            lim2_ky = -lim1_ky;
            % Deform integration path around branch cuts.
            integrationpath_ky = [(-1-1j).*delta_ky, (1+1j).*delta_ky];

            %% Initial sample points
            if(~deformedpath)
                %% Straight path
                % Go to +-inf
                delta = 0.01.*k0;
                lim1 = -100.*k0-1j.*delta;
                lim2 = -lim1;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta];

                % Integration path (dotted)
                %               |
                %               |       N1
                %               |...............
                % --------------:---------------
                % .............:|\
                %      N1       | N2
                %               |

                N1 = 50; % Initial number of points on the horizontal sections
                N2 = 10;  % Initial number of points on the diagonal section through the origin
                newkx = [linspace(lim1, integrationpath(1), N1), ...
                         linspace(integrationpath(1), integrationpath(2), N2+2), ...
                         linspace(integrationpath(2), lim2, N1)];

                % Remove the duplicate elements
                newkx(N1) = [];
                newkx(end-N1) = [];
                realnewkx = real(newkx);
            else
                %% Deformed path
                % Go to -inf * 1j
                delta = 0.01*k0;
                lim1 = -5j.*k0-1.*delta;
                lim2 = -5j.*k0+1.*delta+1.5*k0;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];

                % Integration path (dotted)
                %                |      N3
                %            N2  |.............  N2
                %              \ :             :
                % --------------:+--------------:-
                %              : |              :
                %           N1 : |              : N1
                %              : |              :
                %              V |              V

                N1 = 100; % Initial number of points on the vertical sections
                N2 = 10;  % Initial number of points on the diagonal sections
                N3 = 30; % Initial number of points on the horizontal section
                newkx = [linspace(lim1,               integrationpath(1), N1),   ...
                         linspace(integrationpath(1), integrationpath(2), N2+2), ...
                         linspace(integrationpath(2), integrationpath(3), N3),   ...
                         linspace(integrationpath(3), integrationpath(4), N2+2), ...
                         linspace(integrationpath(4), lim2,               N1)];
                % Remove the duplicate elements
                newkx(N1+N2+2) = [];     % Same as newkx(N1+N2+3)
                newkx(N1) = [];          % Same as newkx(N1+1)
                newkx(end-N1-N2-2) = []; % Same as newkx(end-N1-N2-1)
                newkx(end-N1) = [];      % Same as newkx(end-N1+1)
                realnewkx = real(FiniteArray.UnfoldKVector(newkx, integrationpath));
            end
            
            [kxadapt, realkxadapt, Dinvadapt] = adaptivecalc(@(kx) adaptivecalcfun(f, dy, k0, kx, tlineup_, tlinedown_, z0, wslot, Ny_, walled, lim1_ky, lim2_ky, integrationpath_ky), newkx, realnewkx);
            
            %{
            % Plot all Dinv.
            hFig = figureex;
            for(ny = 1:Ny_)
                for(nyp = 1:Ny_)
                    hAx = subplot(Ny_, Ny_, ny+(nyp-1)*Ny_);
                    hold(hAx, 'on');
                    grid(hAx, 'on');
                    box(hAx, 'on');
                    plot(real(kxadapt)./k0, real(Dinvadapt(:, ny, nyp)));
                    plot(real(kxadapt)./k0, imag(Dinvadapt(:, ny, nyp)));
                end
            end
            %}
            
            Dinvlist{fi} = Dinvadapt;
            kxlist{fi} = kxadapt;
            send(hDataQueue, nan);
        end

        if(~deformedpath)
            % Store the frequency points in the this.D_fs vector.
            % Only once, so do it on the first iteration, where deformedpath is false
            this.D_fs = [this.D_fs, newfs];
            % Update progress bar.
            waitbar(1, hWaitbar, ...
                {sprintf('%.1f%% Precomputing D...', 50), ...
                 sprintf('Interpolating...'), ''});
        else
            % Update progress bar.
            waitbar(1, hWaitbar, ...
                {sprintf('%.1f%% Precomputing D...', 100), ...
                 sprintf('Finalizing...'), ''});
        end
        % Store the calculated Ds and their kx values in the this.Ds and this.D_kxs cells.
        for(fi = 1:Nf)
            f = newfs(fi);
            k0 = 2*pi*f/c0;
            fii = find(this.D_fs == newfs(fi), 1);

            % The D vector is pre-interpolated on a fixed set of points to speed up later
            % interpolation.
            tempkxs = kxlist{fi};
            tempDs = Dinvlist{fi};
            
            if(~deformedpath)
                realtempkxs = real(tempkxs);
            else
                delta = 0.01*k0;
                integrationpath = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];
                realtempkxs = real(FiniteArray.UnfoldKVector(tempkxs, integrationpath));
            end

            for(ny = 0:Ny_-1)
                for(nyp = 0:Ny_-1)
                    this.Dinv_interpolants{deformedpath+1, fii, ny+1, nyp+1} = griddedInterpolant(realtempkxs, tempDs(:, ny+1, nyp+1));
                end
            end
            this.Dinvs{deformedpath+1, fii} = tempDs;
            this.Dinv_kxs{deformedpath+1, fii} = realtempkxs;
        end
        if(~deformedpath)
            progress = 0;
            dispex('Ds: Completed straight integration path in %s.\n', fancyduration(toc(tc)));
            halfdt = toc(tc);
        else
            dispex('Ds: Completed deformed integration path in %s.\n', fancyduration(toc(tc) - halfdt));
        end
    end
    
    function updateWaitbar(~)
        progress = progress + 1;
        waitbar(progress/Nf, hWaitbar, ...
            {sprintf('%.1f%% Precomputing D...', progress/Nf/2*100+deformedpath*50), ...
             sprintf('%i/%i iterations done.', progress, Nf), ...
             sprintf('%i/%i integration paths done.', deformedpath, 2)});
    end
    delete(hWaitbar);
    
    dt = toc(tc);
    dispex('Ds: Completed in %s.\n', fancyduration(dt));
end
function Dinv = adaptivecalcfun(f, dy, k0, kx, tlineup_, tlinedown_, z0, wslot, Ny, walled, lim1_ky, lim2_ky, integrationpath_ky)
    kx(kx == 0) = (1+1j)*eps;
    Dinv = zeros(length(kx), Ny, Ny);
    Ds = zeros(length(kx), Ny);
    Ddowns = zeros(length(kx), 1);
%     Dfss = zeros(length(kx), Ny);
    parfor(kxii = 1:length(kx)) % parfor
        if(kxii == round(length(kx)/2))
        end
        
        D = zeros(1, Ny);
        for(dny = 0:Ny-1)
            D(dny+1) = 1 ./ (2*pi) .* fastintegral(...
                @(kyp) D_Integrand_kyp_old(f, dy, k0, kyp, kx(kxii), tlineup_, tlinedown_, z0, wslot, dny, walled), ...
                lim1_ky, lim2_ky, 'Waypoints', integrationpath_ky);
        end
        K = -1i*sqrt(-(k0.^2-kx(kxii).^2));
        Dfs = [-0.5/k0/z0*(k0^2-kx(kxii).^2).*besselj(0,wslot/4*K)*besselh(0,2,wslot/4*K), ... % self
               -0.5/k0/z0*(k0^2-kx(kxii).^2).*besselh(0,2,K*dy*(1:Ny-1))];                     % mutual
        if(walled)
            D = D + Dfs; % Only up
        else
            D = D + 2.*Dfs; % Times two for up & down
        end
        
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
function v = D_Integrand_kyp(f, dy, k0, kyp, kx, tlineup, tlinedown, z0, wslot, nypp)
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
function integrand = D_Integrand_kyp_old(f, dy, k0, kyp, kx, tlineup, tlinedown, z0, wslot, dny, walled)
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





















