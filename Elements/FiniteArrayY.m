classdef FiniteArrayY < handle
    properties
        unitcell % Unit cell in the finite array
        
        tlineup   % Stratification above the array
        tlinedown % Stratification below the array
        
        Ny % Number of unit cells.
        ay % Excitation amplitude vector
        zfeed % Impedance of the feeding network.
        
        numM % Number of modes to sum (-numM:numM). Default -10:10.
        
        Dmys
        Dmy_fs
        Dmy_th
        Dmy_ph
        
        Zmat
        Zmat_fs
        Zmat_th
        Zmat_ph
        
        ints
        ints_fs
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
        function InitializeDs(this, fs, th, ph)
            if(~isempty(this.Dmy_th) && ~isempty(this.Dmy_ph) && (this.Dmy_th ~= th || this.Dmy_ph ~= ph))
                % If the angle is different, simply wipe it.
                this.Dmys = [];
                this.Dmy_fs = [];
                this.Dmy_th = [];
                this.Dmy_ph = [];
            end
            newfs = setdiff(fs, this.Dmy_fs);
            if(length(newfs) < 1)
                return;
            end
            
            dispex('Ds: Calculating for %i frequencies, %i slots, %i modes.\n', length(newfs), this.Ny, this.numM*2+1);
            tc = tic;
            
            dx = this.unitcell.dx;
            dy = this.unitcell.dy;
            wslot = this.unitcell.wslot;
            % dslot = this.unitcell.dslot;
            walled = this.unitcell.walled;
            
            z0 = Constants.z0;
            
            mxs = [-this.numM:this.numM].';
            nyps = 1:this.Ny;
            Nf = length(newfs);

            tlineup_ = this.tlineup;
            tlinedown_ = this.tlinedown;
            
            [mximat, nypimat, fimat] = ndgrid(1:length(mxs), 1:length(nyps), 1:Nf);
            mximat = mximat(:);
            nypimat = nypimat(:);
            fimat = fimat(:);
            
            Ni = length(fimat);
            
            %% Initialize progress bar
            hDataQueue = parallel.pool.DataQueue;
            hWaitbar = waitbar(0, {'0% Initializing...', ''});
            afterEach(hDataQueue, @updateWaitbar);
            progressStep = 10;
            progress = -progressStep; send(hDataQueue, nan);
            
            Dmys_ = zeros(1, Ni);
            parfor(ni = 1:Ni) % parfor
                f = newfs(fimat(ni)); %#ok<PFBNS> Broadcast variable.
                mx = mxs(mximat(ni)); %#ok<PFBNS> Broadcast variable.
                nyp = nyps(nypimat(ni)); %#ok<PFBNS> Broadcast variable.

                [k0, kx0, ~, ~] = k(f, 1, th, ph);
                kxm = kx0 - 2*pi*mx/dx;
                
                delta = 0.01*k0;
                kxm = kxm + 1j .* delta .* sign(kxm);

                %% Calculate the integral for D
                delta = 0.01.*k0;
                lim1 = -100.*k0-1j.*delta;
                lim2 = -lim1;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta];
                % Loop version - faster
                Dmys_(ni) = -1/(2*pi) .* integral(...
                    @(ky) D_Integrand(f, dy, k0, ky, kxm, tlineup_, tlinedown_, z0, wslot, nyp, walled), ...
                    lim1, lim2, 'Waypoints', integrationpath);
                if(abs(Dmys_(ni)) < eps)
                    dispex('Nope\n');
                end
                % ArrayValued version
%                 Dmy = 1/(2*pi) .* integral(...
%                     @(ky) D_Integrand(f, dy, k0, ky, kxm, tlineup, tlinedown, z0, wslot, nyp), ...
%                     lim1, lim2, 'Waypoints', integrationpath, ...
%                     'ArrayValued', 1);
                if(mod(ni, progressStep) == 0)
                    send(hDataQueue, nan);
                end
            end
            this.Dmy_fs = [this.Dmy_fs, newfs];
            Dmat = reshape(Dmys_, length(mxs), length(nyps), length(newfs));
            for(fi = 1:length(newfs))
                this.Dmys(:, :, this.Dmy_fs == newfs(fi)) = Dmat(:, :, fi);
            end
            this.Dmy_th = th;
            this.Dmy_ph = ph;
            
            function updateWaitbar(~)
                progress = progress + progressStep;
                waitbar(progress/Ni, hWaitbar, ...
                    {sprintf('%.1f%% Calculating D...', progress/Ni*100), ...
                     sprintf('%i/%i iterations done.', progress, Ni)});
            end
            delete(hWaitbar);
            
            dt = toc(tc);
            dispex('Ds: Completed in %s.\n', fancyduration(dt));
        end
        function InitializeZMatrix(this, fs, th, ph)
            if(~isempty(this.Zmat_th) && ~isempty(this.Zmat_ph) && (this.Zmat_th ~= th || this.Zmat_ph ~= ph))
                % If the angle is different, simply wipe it.
                this.Zmat = [];
                this.Zmat_fs = [];
                this.Zmat_th = [];
                this.Zmat_ph = [];
            end
            
            newfs = setdiff(fs, this.Zmat_fs);
            if(length(newfs) < 1)
                return;
            end
            this.InitializeDs(newfs, th, ph);
            
            dispex('Z matrix: Calculating for %i frequencies, %i slots.\n', length(newfs), this.Ny);
            tc = tic;
            
            dx = this.unitcell.dx;
            dy = this.unitcell.dy;
            wslot = this.unitcell.wslot;
            dslot = this.unitcell.dslot;
            
            walled = this.unitcell.walled;
            
            z0 = Constants.z0;
            
            tlinedown_ = this.tlinedown;
            
            mxs = [-this.numM:this.numM].';
            mys = mxs;
            nyps = 1:this.Ny;
            Nf = length(newfs);
            
            %% Initialize progress bar
            hDataQueue = parallel.pool.DataQueue;
            hWaitbar = waitbar(0, {'0% Initializing...', ''});
            afterEach(hDataQueue, @updateWaitbar);
            progressStep = 1;
            progress = -progressStep; send(hDataQueue, nan);
            
            % TODO: Can unpack this and use fi to index, which allows slicing of the variable.
            Dmys_ = this.Dmys;
            Dmy_fs_ = this.Dmy_fs;
            
            ints = zeros(length(mxs), this.Ny, Nf);
            Zmat_ = zeros(this.Ny, this.Ny, Nf);
            % TODO: Can be sliced into fi, kxm, and nyp for more parallelization.
            for(fi = 1:Nf) % parfor
                f = newfs(fi);
                Dmy = Dmys_(:, :, Dmy_fs_ == f); %#ok<PFBNS> Broadcast variable

                [k0, kx0, ~, ~] = k(f, 1, th, ph);
                kxm = kx0 - 2*pi*mxs/dx;
                
                %% Calculate Ddown
                Ddown = 0;
                if(walled)
                    for(myi = 1:length(mys))
                        my = mys(myi);
                        kym =  - 2*pi*my/dy;

                        Vte = 1;
                        Vtm = 1;

                        kr = sqrt(kxm.^2 + kym.^2);
                        isTE = 1;    ztedown = tlinedown_.GetInputImpedance(isTE, f, k0, kr);
                        isTE = 0;    ztmdown = tlinedown_.GetInputImpedance(isTE, f, k0, kr);

                        itedown = 1 ./ ztedown;
                        itmdown = 1 ./ ztmdown;

                        [Ghm] = SpectralGF.hm(z0, k0, kxm, kym, Vtm, Vte, itmdown, itedown);
                        Ghm_xx_down = Ghm.xx;

                        Ddown = Ddown + -1/dy .* Ghm_xx_down .* besselj(0, -2*pi*my/dy*wslot/2);
                    end
                end

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
                    @(ky) Z_Integrand(ky, dy, 1, nyps, Dmy, Ddown), ...
                    lim1, lim2, ...
                    'ArrayValued', 1);
                
                ints(:, :, fi) = int;
                
                %% Perform the sum for Z & build Z matrix
                Zv = 1./dx .* sum(sinc(kxm.*dslot./(2*pi)).^2 .* dy./(2*pi) .* int, 1);
                
                % Build the toeplitz matrix
                Zmat_(:, :, fi) = toeplitz(real(Zv)) + 1j .* toeplitz(imag(Zv));
                
                send(hDataQueue, nan);
            end
            this.ints_fs = [this.ints_fs, newfs];
            for(fi = 1:length(newfs))
                this.ints(:, :, this.ints_fs == newfs(fi)) = ints(:, :, fi);
            end
            this.Zmat_fs = [this.Zmat_fs, newfs];
            for(fi = 1:length(newfs))
                this.Zmat(:, :, this.Zmat_fs == newfs(fi)) = Zmat_(:, :, fi);
            end
            this.Zmat_th = th;
            this.Zmat_ph = ph;
            
            function updateWaitbar(~)
                progress = progress + progressStep;
                waitbar(progress/Nf, hWaitbar, ...
                    {sprintf('%.1f%% Calculating Z Matrix...', progress/Nf*100), ...
                     sprintf('%i/%i frequencies done.', progress, Nf)});
            end
            delete(hWaitbar);
            
            dt = toc(tc);
            dispex('Z matrix: Completed in %s.\n', fancyduration(dt));
            
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
        
        function Zas = GetInputImpedance(this, fs, th, ph)
            this.InitializeDs(fs, th, ph);
            this.InitializeZMatrix(fs, th, ph);
            
            dispex('Active Z: Calculating for %i frequencies, %i slots.\n', length(fs), this.Ny);
            tc = tic;
            
            Nf = length(fs);
            
            Zmat_ = this.Zmat;
            Zmat_fs_ = this.Zmat_fs;
            
            Ny_ = this.Ny;
            ay0_ = this.ay;
            
            dy = this.unitcell.dy;
            
            
            zfeed_ = this.zfeed;
            
            Zas = zeros(this.Ny, Nf);
            parfor(fi = 1:Nf) % parfor
                f = fs(fi);
                
                [~, ~, ky0, ~] = k(f, 1, th, ph);
                ay_ = ay0_ .* exp(-1j .* ky0 .* (1:Ny_) .* dy);
                
                Z = Zmat_(:, :, Zmat_fs_ == f); %#ok<PFBNS> Broadcast variable
                
                % Add ZL
                Zp = Z + eye(Ny_) .* zfeed_;
                
                %% Solve the equivalent circuit and find Zactive
                i = ay_;
                % Solve for v
                v = Zp\(zfeed_*Z*i.'); % Identical to v = pinv(Zp) * (zfeed_*Z*i.');
                % Calculate input impedance
                Yas = i.' ./ v - 1./zfeed_;
                Zas(:, fi) = 1./Yas;
            end
            
            dt = toc(tc);
            dispex('Active Z: Completed in %s.\n', fancyduration(dt));
        end
        function BuildCST(this, project, parentcomponent)
            if(nargin < 3 || isempty(parentcomponent))
                parentcomponent = '';
            else
                if(~strcmp(parentcomponent(end), '/'))
                    parentcomponent = [parentcomponent, '/'];
                end
            end
            componentname = [parentcomponent, 'Array'];
            

            boundary = project.Boundary();
            boundary.Ymin('open');
            boundary.Ymax('open');

            project.StoreParameter('Ny', this.Ny);

            this.unitcell.BuildCST(project, [componentname, '/UnitCell']);

            wcs = project.WCS();
            wcs.Enable(); 
            wcs.Store('Pre-Slot'); 
            wcs.RotateWCS('u', 180); 

            % Build down-stratification. 
            this.tlinedown.BuildCST(project, [componentname, '/Stratification']); 

            wcs.Restore('Pre-Slot'); 

            % Build up-stratification. 
            this.tlineup.BuildCST(project, [componentname, '/Stratification']); 

            wcs.Restore('Pre-Slot'); 
            wcs.Delete('Pre-Slot'); 
            wcs.Disable(); 

            transform = project.Transform();

            transform.Reset();
                transform.Name([componentname, '/UnitCell']);
                transform.Vector(0, 'dy', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Ny-1');
            transform.Transform('Shape', 'Translate');

            transform.Reset();
                transform.Name([componentname, '/Stratification']);
                transform.Vector(0, '-dy', 0);
                transform.MultipleObjects(0);
                transform.GroupObjects(0);
                transform.Repetitions('1');
            transform.Transform('Shape', 'Translate');
            
            transform.Reset();
                transform.Name([componentname, '/Stratification']);
                transform.Vector(0, 'dy', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Ny+1');
            transform.Transform('Shape', 'Translate');
            
            brick = project.Brick();
            brick.Reset();
                brick.Name('SlotPlane');
                brick.Component(componentname);
                brick.Xrange('-dx/2', 'dx/2');
                brick.Yrange('-3/2*dy', '-dy/2');
                brick.Zrange(0, 0);
                brick.Material('PEC');
            brick.Create();
            
            transform.Reset();
                transform.Name([componentname, ':SlotPlane']);
                transform.Vector(0, '(Ny+1)*dy', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('1');
            transform.Transform('Shape', 'Translate');
%             brick.Reset();
%                 brick.Name('SlotPlane2');
%                 brick.Component(componentname);
%                 brick.Xrange('-dx/2', 'dx/2');
%                 brick.Yrange('Ny*dy-dy/2', 'Ny*dy+dy/2');
%                 brick.Zrange(0, 1);
%                 brick.Material('PEC');
%             brick.Create();

            transform.Reset();
                transform.Name('port1 (SlotFeed)');
                transform.Vector(0, 'dy', 0);
                transform.MultipleObjects(1);
                transform.Repetitions('Ny-1');
            transform.Transform('Port', 'Translate');
        end
    end
end

function v = D_Integrand(f, dy, k0, ky, kxm, tlineup, tlinedown, z0, wslot, nyp, walled)
    Vtm = 1;
    Vte = 1;

    kr = sqrt((kxm).^2 + (ky ).^2);
    isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kxm, ky, Vtm, Vte, itmup, iteup);
    Gxxup = Ghm.xx;

    if(~walled)
        isTE = 1;    ztedown = tlinedown.GetInputImpedance(isTE, f, k0, kr);
        isTE = 0;    ztmdown = tlinedown.GetInputImpedance(isTE, f, k0, kr);

        itedown = 1 ./ ztedown;
        itmdown = 1 ./ ztmdown;

        [Ghm] = SpectralGF.hm(z0, k0, kxm, ky, Vtm, Vte, itmdown, itedown);
        Gxxdown = Ghm.xx;
    else
        Gxxdown = 0;
    end

    Gxx = Gxxup + Gxxdown;

    v = Gxx.*besselj(0, ky.*wslot./2).*exp(-1j.*ky.*dy.*(nyp-1));
    
%     % Assume anything that is NaN can be approximated by 0.
%     v(isnan(v)) = 0;
end

function v = Z_Integrand(ky, dy, ny, nyp, Dmy, Ddown)
    s = sum(Dmy .* exp(1j .* (nyp-1) .* ky .* dy), 2) + Ddown;
    v = exp(-1j .* ky .* dy .* abs(ny-nyp)) ./ s;
end
%% Non-arrayvalued version
function v = Z_Integrand2(ky, dy, ny, nyp, nypp, Dmy, Ddown)
    s = sum(Dmy .* exp(1j .* (nypp-1) .* ky .* dy), 1) + Ddown;
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
%     isTE = 1;    ztedown = tlinedown.GetInputImpedance(isTE, f, k0, kr);
%     isTE = 0;    ztmdown = tlinedown.GetInputImpedance(isTE, f, k0, kr);
% 
%     itedown = 1 ./ ztedown;
%     itmdown = 1 ./ ztmdown;
% 
%     [Ghm] = SpectralGF.hm(z0, k0, kxm, ky, Vtm, Vte, itmdown, itedown);
%     GFxxdown = Ghm.xx;
% 
%     v = besselj(0, ky.*wslot./2) .* (GFxxup + GFxxdown) .* exp(1j .* ky .* (nyp - ny) .* dy);
% end

