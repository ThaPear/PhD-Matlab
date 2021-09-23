function Iny = CurrentSpectrum(this, fs, excitation, kx, ny)
    this.InitializeDs(fs);
    this.InitializeZMatrix(fs);

    dedge_ = this.dedge;
    Nx_ = this.Nx;
    Ny_ = this.Ny;
    dx = this.unitcell.dx;
%     dy = this.unitcell.dy;
    wslot = this.unitcell.wslot;
    dslot = this.unitcell.dslot;
    zfeed_ = this.zfeed;

    % Do not deform for now.
    deformedpath = 0;

    % Add zeroes to the excitation matrix for the terminations
    excitation_ = zeros(Nx_+2, Ny_);
    excitation_(2:end-1, :) = excitation;

    basispositions = [-dedge_, (0:Nx_-1)*dx, (Nx_-1)*dx+dedge_];
    x = linspace(min(basispositions), max(basispositions), 20);


    for(fi = 1:length(fs)) % PARFOR?
        f = fs(fi);
        lambda = Constants.c0/f;
        
        D_fi = find(this.D_fs == f, 1);

        g = 5/3 * sqrt(wslot * lambda);
        
        Ffeed = @(kx) sinc(kx .* dslot ./ (2*pi));
        Fright = @(kx) (besselj(0, kx .* g./2) - 1j .* struve(0, kx.*g./2) ...
                        - 2/pi .* sinc(kx.*g./(4*pi)).*exp(-1j.*kx.*g./4)) .* exp(1j.*kx.*g./2);
        % The basis function of the second termination is the reverse of that of the first
        % termination.
        Fleft = @(kx) Fright(-kx);
        
        %% Copied from GetInputImpedance
            excitation_(2:end-1, :) = excitation(:, :, fi);
            f = fs(fi);
            % Retrieve the Z matrix for this frequency
            fii = find(this.Z_fs == f, 1);
            Zmat_ = this.Zmat{fii};
            % Build a matrix with zfeed on the diagonal, but only for non-termination elements
            Zl = repmat([0, ones(1, Nx_)*zfeed_, 0], 1, Ny_);
            Zlmat = diag(Zl);
            % Add the loads to the mutual impedances
            Zp = Zmat_ + Zlmat;
            % Norton
            % Determine current through the elements
            i = Zp\(Zlmat * excitation_(:));
            % Determine voltage on the elements themselves.
            v = Zmat_ * i;
        % End of copy
        i = reshape(i, Nx_+2, Ny_);
        v = reshape(v, Nx_+2, Ny_);
        
        basispositions = [-dedge_, (0:Nx_-1)*dx, (Nx_-1)*dx+dedge_];
        basisfunctions = cell(1,Nx_+2);
        basisfunctions{1} = Fleft;
        basisfunctions(2:Nx_+1) = {Ffeed};
        basisfunctions{end} = Fright;
        
        Iny = zeros(size(kx));
        for(nyp = 0:Ny_-1)
            for(nxp = -1:Nx_)
                Fp = basisfunctions{nxp+2};
                xp = basispositions(nxp+2);

                Dinv_interpolant = this.Dinv_interpolants{deformedpath+1, D_fi, ny+1, nyp+1};
                % TODO: If path is deformed this won't work.
                interpolationkx = real(kx);
                Dinv = Dinv_interpolant(interpolationkx);
    
                % 1x1 slot in free-space closed-form Dinv
%                 k0 = 2*pi*fs/3e8;
%                 z0 = 120*pi;
%                 K = -1i*sqrt(-(k0.^2-kx.^2));
%                 Dinv_fs = 1./(2*-0.5/k0/z0*(k0^2-kx.^2).*besselj(0,wslot/4*K).*besselh(0,2,wslot/4*K));

                Fpkx = Fp(kx);
                for(kxi = 1:length(kx))
%                     Vny(kxi) = Vny(kxi) + -Dinv(kxi) .* (i(nxp+2,nyp+1) - 1/zfeed_ * v(nxp+2, nyp+1)) .* Fp(kx(kxi)) .* exp(1j .* kx(kxi) .* xp);
                    Iny(kxi) = Iny(kxi) + -(i(nxp+2,nyp+1)) .* Fpkx(kxi) .* exp(1j .* kx(kxi) .* xp);
                end
            end
        end
    end
end
        
%         % Determine the D values for the given kxs.
%         fii = find(this.D_fs == f, 1);
%         % If the path is straight, there's no need to unfold.
%         interpolationkx = real(kx);
%         % Interpolate the D vectors at the desired kx values.
%         D = zeros(length(kx), Ny_);
%         for(nypp = 0:Ny_-1)
%     %         D(:, nypp+1) = interp1(this.D_kxs{typei, fii, nypp+1}, this.Ds{typei, fii, nypp+1}, real(kx));
%         end
%         
% 
%         %% Basis functions
%         Ffeed = @(kx) sinc(kx .* dslot ./ (2*pi));
%         Fright = @(kx) ( besselj(0, kx .* g./2) - 1j .* struve(0, kx.*g./2) ...
%                         - 2/pi .* sinc(kx.*g./(4*pi)).*exp(-1j.*kx.*g./4)) .* exp(1j.*kx.*g./2);
%         % The basis function of the second termination is the reverse of that of the first
%         % termination.
%         Fleft = @(kx) Fright(-kx); 
% 
%         Fs = cell(1,Nx_+2);
%         Fs{1} = Fleft;
%         Fs(2:Nx_+1) = {Ffeed};
%         Fs{end} = Fright;
% 
%         sumVJ = 0;
%         for(ny = 0:Ny_-1)
%             V_kx_ny = zeros(size(kx));
%             for(kxi = 1:length(kx)) % PARFOR?
%                 if(Ny_ == 1) % Indexing D(kxi, :) fails when Ny_ is 1.
%                     Dparam = D(kxi);
%                 else
%                     Dparam = D(kxi, :);
%                 end
%                 V_kx_ny(kxi) = dy ./ (2*pi) .* fastintegral(@(ky) V_Integrand_ky(kx(kxi), ky, Nx_, Ny_, ny, dy, Dparam, Fs, basispositions, i), -pi/dy, pi./dy);
%             end
% 
%             sumVJ = sumVJ + V_kx_ny .* besselj(0, ky .* wslot ./ 2) .* exp(1j .* ky .* ny .* dy);
%         end
% 
%         M(fi, :) = sumVJ;
%     end
% end
