classdef FiniteArrayX < handle
    properties
        unitcell % Unit cell in the finite array
        
        tlineup   % Stratification above the array
        tlinedown % Stratification below the array
        
        Nx % Number of unit cells.
        ax % Excitation amplitude vector
        dedge % Distance to the slot termination.
        zfeed % Impedance of the feeding network.
        
        numM % Number of modes to sum (-numM:numM). Default -10:10.
        
        Zmat % Z-matrices
        Z_fs % Frequencies corresponding to the matrices in Zmat
    end
    methods
        function this = FiniteArrayX(unitcell, tlineup, tlinedown, Nx, ax, dedge, zfeed)
            this.unitcell = unitcell;
            this.tlineup = tlineup;
            this.tlinedown = tlinedown;
            this.Nx = Nx;
            this.dedge = dedge;
            this.zfeed = zfeed;
            
            if(nargin < 3 || length(ax) ~= Nx)
                error('Invalid number of weights specified');
            end
            this.ax = ax;
            
            this.numM = 10;
            
            this.Zmat = [];
            this.Z_fs = [];
        end
        function Zas = GetInputImpedance(this, fs, th, ph)
            this.numM = 10;
            
            %% Modes
            my = [-this.numM:this.numM];
            
            Nf = length(fs);
            Nx_ = this.Nx;
            ax_ = this.ax;
            dedge_ = this.dedge;
            Zfeed_ = this.zfeed;
            
            dx = this.unitcell.dx;
            dy = this.unitcell.dy;
            wslot = this.unitcell.wslot;
            dslot = this.unitcell.dslot;
            tlineup_ = this.tlineup;
            tlinedown_ = this.tlinedown;
            walled = this.unitcell.walled;
            
            z0 = Constants.z0;
            c0 = Constants.c0;
            
%             this = parallel.pool.Constant(this);
            %% Initialize progress bar
            hDataQueue = parallel.pool.DataQueue;
            hWaitbar = waitbar(0, {'0% Initializing...', ''});
            afterEach(hDataQueue, @updateWaitbar);
            progress = -1; send(hDataQueue, nan);
            
            Zas = zeros(Nx_, Nf);
            for(fi = 1:Nf)
                f = fs(fi);
                
                [k0, ~, ky0, ~] = k(f, 1, th, ph);
                
                lambda = c0/f;
                g = 5/3 * sqrt(wslot * lambda);
                
                Ffeed = @(kx) sinc(kx .* dslot ./ (2*pi));
                Fleft = @(kx)  ( besselj(0,  kx .* g./2) - 1j .* struve(0,  kx.*g./2) ...
                                - 2/pi .* sinc( kx.*g./(4*pi)).*exp(-1j.* kx.*g./4)) .* exp(1j.* kx.*g./2);
                Fright = @(kx) ( besselj(0, -kx .* g./2) - 1j .* struve(0, -kx.*g./2) ...
                                - 2/pi .* sinc(-kx.*g./(4*pi)).*exp(-1j.*-kx.*g./4)) .* exp(1j.*-kx.*g./2);
                
                %% Self impedances
                % Go to +-inf
                delta = 0.01.*k0;
                lim1 = -100.*k0-1j.*delta;
                lim2 = -lim1;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta];
                
                % Feed to itself
                xi = 1 .* dx;
                xj = 1 .* dx; % 1 result
                Z1n = zeros(1,Nx_);
                Z1n(1) = -1./(2*pi) .* integral(...
                    @(kx) Z_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, ...
                                                 Ffeed, Ffeed, tlineup_, tlinedown_, xi, xj), ...
                    lim1, lim2, 'Waypoints', integrationpath);
                
                % Edge to itself
                xi = 1 .* dx - dedge_;
                xj = 1 .* dx - dedge_; % 1 result
                Zrr = -1./(2*pi) .* integral(...
                    @(kx) Z_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, ...
                                                 Fleft, Fleft, tlineup_, tlinedown_, xi, xj), ...
                    lim1, lim2, 'Waypoints', integrationpath);
                
                %% Mutual impedances
                % Go to -inf * 1j
                delta = 0.01*k0;
                lim1_2 = -20j.*k0-1.*delta;
                lim2_2 = -20j.*k0+1.*delta+1.5*k0;
                % Deform integration path around branch cuts.
                integrationpath_2 = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];
                
                % Feed to feed
                xi = 1 .* dx;
                for(j = 2:Nx_) % N-1 results
                    xj = j .* dx;
                    if(abs(xi - xj) < 2*g)
                        Z1n(j) = -1./(2*pi) .* integral(...
                            @(kx) Z_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, ...
                                                         Ffeed, Ffeed, tlineup_, tlinedown_, xi, xj), ...
                            lim1, lim2, 'Waypoints', integrationpath);
                    else
                        Z1n(j) = -1./(2*pi) .* integral(...
                            @(kx) Z_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, ...
                                                         Ffeed, Ffeed, tlineup_, tlinedown_, xi, xj), ...
                            lim1_2, lim2_2, 'Waypoints', integrationpath_2);
                    end
                end
                % Edge to feed
                xi = 1 .* dx - dedge_;
                Zrn = zeros(1,Nx_);
                for(j = 1:Nx_) % N results
                    xj = j .* dx;
                    if(abs(xi - xj) < 2*g)
                        Zrn(j) = -1./(2*pi) .* integral(...
                            @(kx) Z_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, ...
                                                         Fleft, Ffeed, tlineup_, tlinedown_, xi, xj), ...
                            lim1, lim2, 'Waypoints', integrationpath);
                    else
                        Zrn(j) = -1./(2*pi) .* integral(...
                            @(kx) Z_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, ...
                                                         Fleft, Ffeed, tlineup_, tlinedown_, xi, xj), ...
                            lim1_2, lim2_2, 'Waypoints', integrationpath_2);
                    end
                end
                % Edge to edge
                xi = 1 .* dx - dedge_;
                xj = Nx_ .* dx + dedge_; % 1 result
                if(abs(xi - xj) < 2*g)
                    Zrl = -1./(2*pi) .* integral(...
                        @(kx) Z_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, ...
                                                     Fleft, Fright, tlineup_, tlinedown_, xi, xj), ...
                        lim1, lim2, 'Waypoints', integrationpath);
                else
%                     Zrl_oldpath = -1./(2*pi) .* integral(...
%                         @(kx) Z_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, ...
%                                                      Fedge, Fedge2, tlineup_, tlinedown_, xi, xj), ...
%                         lim1, lim2, 'Waypoints', integrationpath);
                    Zrl = -1./(2*pi) .* integral(...
                        @(kx) Z_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, ...
                                                     Fleft, Fright, tlineup_, tlinedown_, xi, xj), ...
                        lim1_2, lim2_2, 'Waypoints', integrationpath_2);
                end
                
                %% Fill toeplitz Z matrix
                Zv = [Zrr, fliplr(Zrn), Zrl];
                Z = toeplitz(real(Zv)) + 1j .* toeplitz(imag(Zv));
                Z(2:end-1, 2:end-1) = toeplitz(real(Z1n)) + 1j .* toeplitz(imag(Z1n));
                
                %% Construct Zprime matrix
                % Diagonal matrix with ZL.
                ZL_mat = Zfeed_ .* eye(Nx_+2);
                % Edge elements have no load impedance.
                ZL_mat([1 Nx_+2], [1 Nx_+2]) = 0;
                Zp = Z + ZL_mat;
                
                %% Define voltages
                v = [0 ax_ 0].'; % Column vector.
                % Multiply matrices.
                vant = pinv(Zp) * Z * v;
                iant = Zp\v; % identical to iant = pinv(Zp) * v;
                
                % Store z-matrix.
                this.Zmat(fi, :, :) = Z;
                this.Z_fs(fi) = f;
                
                %%
                Zn = vant./iant;
                
                Zas(:, fi) = Zn(2:end-1);
                send(hDataQueue, nan);
            end
            function updateWaitbar(~)
                progress = progress + 1;
                waitbar(progress/Nf, hWaitbar, ...
                    {sprintf('%.1f%% Calculating Z Matrix...', progress/Nf*100), ...
                     sprintf('%i/%i frequencies done.', progress, Nf)});
            end
            delete(hWaitbar);
        end
        function [x, v] = Voltage(this, f, th, ph)
            fi = find(this.Z_fs == f, 1);
            if(isempty(fi))
                this.GetInputImpedance(f);
            end
            
            Nx_ = this.Nx;
            ax_ = this.ax;
            dedge_ = this.dedge;
            Zfeed_ = this.zfeed;
%             
            dx = this.unitcell.dx;
            dy = this.unitcell.dy;
            wslot = this.unitcell.wslot;
            dslot = this.unitcell.dslot;
            tlineup_ = this.tlineup;
            tlinedown_ = this.tlinedown;
            walled = this.unitcell.walled;
            
            my = [-this.numM:this.numM];
            
            z0 = Constants.z0;
            c0 = Constants.c0;
            [k0, ~, ky0, ~] = k(f, 1, th, ph);
            lambda = c0/f;
            g = 5/3 * sqrt(wslot * lambda);
            
            % Retrieve Z-matrix for this frequency.
            Z = squeeze(this.Zmat(fi, :, :));
            ZL_mat = Zfeed_ .* eye(Nx_+2);
            
            ZL_mat([1 Nx_+2], [1 Nx_+2]) = 0;
            Zp = Z + ZL_mat;
            
            % Solve equivalent circuit
            v = [0 ax_ 0].'; % Column vector.
            % Multiply matrices.
            vant = pinv(Zp) * Z * v;
            iant = Zp\v; % identical to iant = pinv(Zp) * v;
            
            % Go to +-inf
            delta = 0.01.*k0;
            lim1 = -100.*k0-1j.*delta;
            lim2 = -lim1;
            % Deform integration path around branch cuts.
            integrationpath = [(-1-1j).*delta, (1+1j).*delta];
            
            Ffeed = @(kx) sinc(kx .* dslot ./ (2*pi));
            Fleft = @(kx)  ( besselj(0,  kx .* g./2) - 1j .* struve(0,  kx.*g./2) ...
                            - 2/pi .* sinc( kx.*g./(4*pi)).*exp(-1j.* kx.*g./4)) .* exp(1j.* kx.*g./2);
            Fright = @(kx) ( besselj(0, -kx .* g./2) - 1j .* struve(0, -kx.*g./2) ...
                            - 2/pi .* sinc(-kx.*g./(4*pi)).*exp(-1j.*-kx.*g./4)) .* exp(1j.*-kx.*g./2);
                        
            xs = [1*dx-dedge_, (1:Nx_)*dx, Nx_*dx+dedge_];
            Fs = cell(1,Nx_+2);
            Fs{1} = Fleft;
            Fs(2:end-1) = {Ffeed};
            Fs{end} = Fright;
            
            x = linspace(min(xs)-g, max(xs)+g, 100);
            
            for(xi = 1:length(x))
                v(xi) = 1/(2*pi) .* fastintegral(@(kx) V_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, ...
                                                                   walled, Fs, tlineup_, tlinedown_, ...
                                                                   xs, Nx_, iant, x(xi), dslot), ...
                        lim1, lim2, 'Waypoints', integrationpath);
            end
        end
    end
end

function v = Z_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, Fi, Fj, tlineup, tlinedown, xi, xj)
    global kxs;
    kxs = [kxs, kx];

    %% Solve Slot equivalent circuit.
    Vte = 1;
    Vtm = 1;
    
    kym = ky0 - 2*pi*my/dy;
    [kxmat, kymat] = meshgrid(kx, kym);

    % Up stratification
    kr = sqrt((kxmat).^2 + (kymat).^2);
    isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kxmat, kymat, Vtm, Vte, itmup, iteup);
    GFxxup = Ghm.xx;
    Dinfup = 1./dy .* sum(GFxxup .* besselj(0, kymat .* wslot ./ 2), 1);
    

    % Down stratification
    if(walled)
        % Walls suppress ky0.
        kym =     - 2*pi*my/dy;
        [kxmat, kymat] = meshgrid(kx, kym);
        kr = sqrt((kxmat).^2 + (kymat).^2);
    end
    isTE = 1;    ztedown = tlinedown.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmdown = tlinedown.GetInputImpedance(isTE, f, k0, kr);

    itedown = 1 ./ ztedown;
    itmdown = 1 ./ ztmdown;

    [Ghm] = SpectralGF.hm(z0, k0, kxmat, kymat, Vtm, Vte, itmdown, itedown);
    GFxxdown = Ghm.xx;
    Dinfdown = 1./dy .* sum(GFxxdown .* besselj(0, kymat .* wslot ./ 2), 1);
    
    D = Dinfup + Dinfdown;
    
    global Ds;
    Ds = [Ds, D];
    
    v = (Fi(kx) .* Fj(-kx) ./ D) .* exp(-1j .* kx .* abs(xi-xj));
end

function v = V_Integrand(f, k0, z0, kx, ky0, my, dy, wslot, walled, Fs, tlineup, tlinedown, xs, Nx_, i, x, dslot)
    %% Solve Slot equivalent circuit.
    Vte = 1;
    Vtm = 1;
    
    kym = ky0 - 2*pi*my/dy;
    [kxmat, kymat] = meshgrid(kx, kym);

    % Up stratification
    kr = sqrt((kxmat).^2 + (kymat).^2);
    isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kxmat, kymat, Vtm, Vte, itmup, iteup);
    GFxxup = Ghm.xx;
    Dinfup = 1./dy .* sum(GFxxup .* besselj(0, kymat .* wslot ./ 2), 1);
    

    % Down stratification
    if(walled)
        % Walls suppress ky0.
        kym =     - 2*pi*my/dy;
        [kxmat, kymat] = meshgrid(kx, kym);
        kr = sqrt((kxmat).^2 + (kymat).^2);
    end
    isTE = 1;    ztedown = tlinedown.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmdown = tlinedown.GetInputImpedance(isTE, f, k0, kr);

    itedown = 1 ./ ztedown;
    itmdown = 1 ./ ztmdown;

    [Ghm] = SpectralGF.hm(z0, k0, kxmat, kymat, Vtm, Vte, itmdown, itedown);
    GFxxdown = Ghm.xx;
    Dinfdown = 1./dy .* sum(GFxxdown .* besselj(0, kymat .* wslot ./ 2), 1);
    
    D = Dinfup + Dinfdown;
    
    sumiF = 0;
    for(nxp = 0:Nx_+1)
        sumiF = sumiF + i(nxp+1) .* Fs{nxp+1}(kx) .* exp(1j .* kx .* xs(nxp+1));
    end
%     nxp = 2;
%     sumiF = sumiF + i(nxp+1) .* sinc(kx.*dslot./(2*pi)) .* exp(1j .* kx .* xs(nxp+1));
    
    v = -1 ./ D .* sumiF .* exp(-1j .* kx .* x);
%     v = sumiF .* exp(-1j .* kx .* x);
    
%     figureex;
%     plot(real(kx), abs(sumiF));
end
















