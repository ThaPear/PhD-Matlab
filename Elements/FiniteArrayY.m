classdef FiniteArrayY
    properties
        unitcell % Unit cell in the finite array
        
        tlineup   % Stratification above the array
        tlinedown % Stratification below the array
        
        Ny % Number of unit cells.
        ay % Excitation amplitude vector
        zfeed % Impedance of the feeding network.
        
        numM % Number of modes to sum (-numM:numM). Default -10:10.
    end
    methods
        function this = FiniteArrayY(unitcell, tlineup, tlinedown, Ny, ay, zfeed)
            this.unitcell = unitcell;
            this.tlineup = tlineup;
            this.tlinedown = tlinedown;
            this.Ny = Ny;
            this.zfeed = zfeed;
            
            if(nargin < 3 || length(ay) ~= Ny)
                error('Invalid number of weights specified');
            end
            this.ay = ay;
            
            this.numM = 20;
        end
        function Zas = GetInputImpedance(this, fs, th, ph)
            mx = [-this.numM:this.numM].';
            
            dx = this.unitcell.dx;
            dy = this.unitcell.dy;
            wslot = this.unitcell.wslot;
            dslot = this.unitcell.dslot;
            
            z0 = Constants.z0;
            ay_ = this.ay;
            nyp = 1:this.Ny;
            
            Nf = length(fs);
            Ny_ = this.Ny;
            zfeed_ = this.zfeed;

            tlineup_ = this.tlineup;
            tlinedown_ = this.tlinedown;
            %% Initialize progress bar
            hDataQueue = parallel.pool.DataQueue;
            hWaitbar = waitbar(0, {'0% Initializing...', ''});
            afterEach(hDataQueue, @updateWaitbar);
            progress = -1; send(hDataQueue, nan);
            
            Zas = zeros(this.Ny, Nf);
            parfor(fi = 1:Nf)
                f = fs(fi);
                Z = zeros(1, Ny_);

                [k0, kx0, ~, ~] = k(f, 1, th, ph);
                kxm = kx0 - 2*pi*mx/dx;

                %% Calculate the integral for D
                delta = 0.01.*k0;
                lim1 = -200.*k0-1j.*delta;
                lim2 = -lim1;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta];
                % Loop version - faster
                Dmy = zeros(length(kxm), length(nyp));
                for(kxmi = 1:length(kxm))
                    for(nypi = 1:length(nyp))
                        Dmy(kxmi, nypi) = 1/(2*pi) .* integral(...
                            @(ky) D_Integrand(f, dy, k0, ky, kxm(kxmi), tlineup_, tlinedown_, z0, wslot, nyp(nypi)), ...
                            lim1, lim2, 'Waypoints', integrationpath);
                    end
                end
                % ArrayValued version
%                 Dmy = 1/(2*pi) .* integral(...
%                     @(ky) D_Integrand(f, dy, k0, ky, kxm, tlineup, tlinedown, z0, wslot, nyp), ...
%                     lim1, lim2, 'Waypoints', integrationpath, ...
%                     'ArrayValued', 1);

                %% Calculate the integral for Z
                lim1 = -pi/dy;
                lim2 = -lim1;
                % Loop version
%                 int = zeros(length(kxm), length(nyp));
%                 for(kxmi = 1:length(kxm))
%                     for(nypi = 1:length(nyp))
%                         int(kxmi, nypi) = integral(...
%                             @(ky) Z_Integrand2(ky, dy, 1, nyp(nypi), nyp.', Dmy(kxmi, :).'), ...
%                             lim1, lim2);
%                     end
%                 end
                % ArrayValued version - equally fast but simpler.
                int = integral(...
                    @(ky) Z_Integrand(ky, dy, 1, nyp, Dmy), ...
                    lim1, lim2, ...
                    'ArrayValued', 1);
                
                %% Perform the sum for Z & build Z matrix
                Zv = 1./dx .* sum(-sinc(kxm.*dslot./(2*pi)).^2 .* dy./(2*pi) .* int, 1);
                
                % Build the toeplitz matrix
                Z = toeplitz(real(Zv)) + 1j .* toeplitz(imag(Zv));
                % Add ZL
                Zp = Z + eye(Ny_) .* zfeed_;
                
                %% Solve the equivalent circuit and find Zactive
                i = ay_;
                % Solve for v
                v = Zp\(zfeed_*Z*i.'); % Identical to v = pinv(Zp) * (zfeed_*Z*i.');
                % Calculate input impedance
                Yas = i.' ./ v - 1./zfeed_;
                Zas(:, fi) = 1./Yas;
                
                send(hDataQueue, nan);
            end
            function updateWaitbar(~)
                progress = progress + 1;
                waitbar(progress/Nf, hWaitbar, ...
                    {sprintf('%.1f%% Calculating Z Matrix...', progress/Nf*100), ...
                     sprintf('%i/%i frequencies done.', progress, Nf)});
            end
            delete(hWaitbar);
            
            %% Old method (Thesis Riccardo)
%             Zas = zeros(this.Ny, length(fs));
%             for(fi = 1:length(fs))
%                 Zs = zeros(1, Ny_);
%                 for(ny = 1:Ny_)
% %                 dispex('Calculating unit cell %i.\n', ny);
%                     f = fs(fi);
%                     
%                     [k0, kx0, ~, ~] = k(f, 1, th, ph);
%                     
%                     delta = 0.01.*k0;
%                     lim1 = -50.*k0-1j.*delta;
%                     lim2 = -lim1;
%                     
%                     [gridnyp, gridmx] = meshgrid(nyp, mx_lin);
% 
%                     % Up stratification
%                     int = 1 ./ (2*pi) .* integral(...
%                         @(ky) Z_Integrand(f, k0, z0, kx0, gridmx, ky, dx, dy, ...
%                                           wslot, tlineup, tlinedown, ny, gridnyp), ...
%                         lim1, lim2, 'Waypoints', [(-1-1j).*delta, (1+1j).*delta], ...
%                         'ArrayValued', 1);
%                     
%                     D = sum(ay_./ay_(ny).*int, 2);
%                     
%                     kxm = kx0 - 2.*pi.*mx_lin./dx;
%                     Z = -1./dx .* sum(sinc(kxm .* dslot ./ (2*pi)).^2 ./ D.');
%                     Zs(ny) = Z;
%                 end
%                 Zas(:, fi) = Zs;
%             end
        end
    end
end

function v = D_Integrand(f, dy, k0, ky, kxm, tlineup, tlinedown, z0, wslot, nyp)
    Vtm = 1;
    Vte = 1;

    kr = sqrt((kxm).^2 + (ky ).^2);
    isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kxm, ky, Vtm, Vte, itmup, iteup);
    Gxxup = Ghm.xx;

    isTE = 1;    zdownte = tlinedown.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    zdowntm = tlinedown.GetInputImpedance(isTE, f, k0, kr);

    itedown = 1 ./ zdownte;
    itmdown = 1 ./ zdowntm;

    [Ghm] = SpectralGF.hm(z0, k0, kxm, ky, Vtm, Vte, itmdown, itedown);
    Gxxdown = Ghm.xx;

    Gxx = Gxxup + Gxxdown;

    v = Gxx.*besselj(0, ky.*wslot./2).*exp(-1j.*ky.*dy.*(nyp-1));
end

function v = Z_Integrand(ky, dy, ny, nyp, Dmy)
    s = sum(Dmy .* exp(1j .* (nyp-1) .* ky .* dy), 2);
    v = exp(-1j .* ky .* dy .* abs(ny-nyp)) ./ s;
end
%% Non-arrayvalued version
function v = Z_Integrand2(ky, dy, ny, nyp, nypp, Dmy)
    s = sum(Dmy .* exp(1j .* (nypp-1) .* ky .* dy), 1);
    v = exp(-1j .* ky .* dy .* abs(ny-nyp)) ./ s;
end
%% Old method (Thesis Riccardo)
% function v = Z_Integrand(f, k0, z0, kx0, gridmx, ky, dx, dy, wslot, tlineup, tlinedown, ny, gridnyp)
%     Vtm = 1;
%     Vte = 1;
% 
%     kxm = kx0 - 2*pi*mx/dx;
%     kr = sqrt((kxm).^2 ...
%             + (ky ).^2);
%     isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
%     isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);
% 
%     iteup = 1 ./ zteup;
%     itmup = 1 ./ ztmup;
% 
%     [Ghm] = SpectralGF.hm(z0, k0, kxm, ky, Vtm, Vte, itmup, iteup);
%     GFxxup = Ghm.xx;
%     
%     isTE = 1;    zdownte = tlinedown.GetInputImpedance(isTE, f, k0, kr);
%     isTE = 0;    zdowntm = tlinedown.GetInputImpedance(isTE, f, k0, kr);
% 
%     itedown = 1 ./ zdownte;
%     itmdown = 1 ./ zdowntm;
% 
%     [Ghm] = SpectralGF.hm(z0, k0, kxm, ky, Vtm, Vte, itmdown, itedown);
%     GFxxdown = Ghm.xx;
% 
%     v = besselj(0, ky.*wslot./2) .* (GFxxup + GFxxdown) .* exp(1j .* ky .* (nyp - ny) .* dy);
% end

