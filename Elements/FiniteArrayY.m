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
        
        Ymat
        Ymat_fs
        Ymat_th
        Ymat_ph
        
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
            if(th == 0)
                error('Cannot have theta 0.');
            end
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
%             D = zeros(1, Ni);
            parfor(ni = 1:Ni) % parfor
%                 dispex('%i', ni);
                f = newfs(fimat(ni)); %#ok<PFBNS> Broadcast variable.
                mx = mxs(mximat(ni)); %#ok<PFBNS> Broadcast variable.
                nyp = nyps(nypimat(ni)); %#ok<PFBNS> Broadcast variable.

                [k0, kx0, ~, ~] = k(f, 1, th, ph);
                kxm = kx0 - 2*pi*mx/dx;

                %% Calculate the integral for D
                delta = 0.01.*k0;
                lim1 = -500.*k0-1j.*delta;
                lim2 = -lim1;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta];
                
                Dmy = 1/(2*pi) .* integral(...
                    @(ky) D_Integrand(f, dy, k0, ky, kxm, tlineup_, tlinedown_, z0, wslot, nyp, walled), ...
                    lim1, lim2, 'Waypoints', integrationpath);
                
                K = -1i*sqrt(-(k0.^2-kxm.^2));
                if(nyp == 1)
                    Dfs = -0.5/k0/z0*(k0^2-kxm.^2).*besselj(0,wslot/4*K)*besselh(0,2,wslot/4*K);
                else
                    Dfs = -0.5/k0/z0*(k0^2-kxm.^2).*besselh(0,2,K*dy*(nyp-1));
                end
                if(~walled)
                    Dfs = Dfs .* 2;
                end
                Dmys_(ni) = Dmy + Dfs;
                
                % Free-space analytical D
%                 K = -1i.*sqrt(-(k0.^2-kxm.^2));
%                 if nyp==1
%                     Dmys_(ni) = -0.5./k0./z0.*(k0.^2-kxm.^2).*besselj(0,wslot./4.*K)*besselh(0,2,wslot./4.*K);
%                 else
%                     Dmys_(ni) = -1./k0./z0.*(k0.^2-kxm.^2).*besselh(0,2,K.*dy.*(nyp-1));
%                 end

                if(mod(ni, progressStep) == 0)
                    send(hDataQueue, nan);
                end
            end
            %% Store D matrix
            this.Dmy_fs = [this.Dmy_fs, newfs];
            Dmat = reshape(Dmys_, length(mxs), length(nyps), length(newfs));
            % Multiply all non-self Ds by two.
            Dmat(:, 2:end, :) = Dmat(:, 2:end, :) .* 2;
            for(fi = 1:length(newfs))
                this.Dmys(:, :, this.Dmy_fs == newfs(fi)) = Dmat(:, :, fi);
            end
            this.Dmy_th = th;
            this.Dmy_ph = ph;
            

            %{
            % Compare to free-space analytical D
                Dmat2 = reshape(D, length(mxs), length(nyps), length(newfs));
                figureex;
                    repeatcolormap(2);
                    plot(real(Dmat(:,1,1)));
                    plot(imag(Dmat(:,1,1)), '--');
                    plot(real(Dmat2(:,1,1)));
                    plot(imag(Dmat2(:,1,1)), '--');
                figureex;
                    repeatcolormap(2);
                    plot(real(Dmat(:,2,1)));
                    plot(imag(Dmat(:,2,1)), '--');
                    plot(real(Dmat2(:,2,1)));
                    plot(imag(Dmat2(:,2,1)), '--');
            %}
            
            %% Waitbar update
            function updateWaitbar(~)
                progress = progress + progressStep;
                waitbar(progress/Ni, hWaitbar, ...
                    {sprintf('%.1f%% Calculating D...', progress/Ni*100), ...
                     sprintf('%i/%i iterations done.', progress, Ni)});
            end
            delete(hWaitbar);
            %%
            dt = toc(tc);
            dispex('Ds: Completed in %s.\n', fancyduration(dt));
        end
        function InitializeYMatrix(this, fs, th, ph)
            if(th == 0)
                error('Cannot have theta 0.');
            end
            if(~isempty(this.Ymat_th) && ~isempty(this.Ymat_ph) && (this.Ymat_th ~= th || this.Ymat_ph ~= ph))
                % If the angle is different, simply wipe it.
                this.Ymat = [];
                this.Ymat_fs = [];
                this.Ymat_th = [];
                this.Ymat_ph = [];
            end
            
            newfs = setdiff(fs, this.Ymat_fs);
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
            zfeed_ = this.zfeed;
            
            tlineup_ = this.tlineup;
            tlinedown_ = this.tlinedown;
            
            Ny_ = this.Ny;
            
            mxs = [-this.numM:this.numM];
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
            
            Ymat_ = zeros(this.Ny, this.Ny, Nf);
            % TODO: Can be sliced into fi, kxm, and nyp for more parallelization.
            parfor(fi = 1:Nf) % parfor
                f = newfs(fi);
                Dmy = Dmys_(:, :, Dmy_fs_ == f).'; %#ok<PFBNS> Broadcast variable

                [k0, kx0, ky0, ~] = k(f, 1, th, ph);
                kxm = kx0 - 2*pi*mxs/dx;
                
                %% Calculate infinite case
                Dinfdown = 0;
                Dinfup = 0;
                for(myi = 1:length(mys))
                    my = mys(myi);
                    kym = ky0 - 2*pi*my/dy;
                    
                    if(myi == 21)
                        aaaa = 1;
                    end

                    Vte = 1;
                    Vtm = 1;

                    % Up
                    kr = sqrt(kxm.^2 + kym.^2);
                    isTE = 1;    zteup = tlineup_.GetInputImpedance(isTE, f, k0, kr);
                    isTE = 0;    ztmup = tlineup_.GetInputImpedance(isTE, f, k0, kr);

                    iteup = Vte ./ zteup;
                    itmup = Vtm ./ ztmup;

                    [Ghm] = SpectralGF.hm(z0, k0, kxm, kym, Vtm, Vte, itmup, iteup);
                    Ghm_xx_up = Ghm.xx;
                    
                    Dinfup = Dinfup + 1/dy .* Ghm_xx_up .* besselj(0, -2*pi*my/dy*wslot/2);

                    % Down
                    if(walled)
                        kym = - 2*pi*my/dy;
                    end
                    kr = sqrt(kxm.^2 + kym.^2);
                    isTE = 1;    ztedown = tlinedown_.GetInputImpedance(isTE, f, k0, kr);
                    isTE = 0;    ztmdown = tlinedown_.GetInputImpedance(isTE, f, k0, kr);

                    itedown = Vte ./ ztedown;
                    itmdown = Vtm ./ ztmdown;

                    [Ghm] = SpectralGF.hm(z0, k0, kxm, kym, Vtm, Vte, itmdown, itedown);
                    Ghm_xx_down = Ghm.xx;

                    Dinfdown = Dinfdown + 1/dy .* Ghm_xx_down .* besselj(0, -2*pi*my/dy*wslot/2);
                end
%                 Dinf = Dinfdown + Dinfup;
%                 Zinf   =  -1./dx .* sum(sinc(kxm .* dslot ./ (2.*pi)).^2 ./ (Dinf));
%                 Vinf = -sinc(kxm*dslot/2/pi) ./ Dinf ./ Zinf;

                %% Determine Ddown
                if(walled)
                    Ddown = Dinfdown;
                else
                    Ddown = zeros(size(kxm));
                end

                % Only add Ddown to the self.
                mulDdown = zeros(Ny_, 1);
                mulDdown(1) = 1;

                %% Calculate single-slot case
                Dsingle = Dmy(1,:) + Ddown;
                Zsingle = 1 ./ dx .* sum(-sinc(kxm .* dslot/(2*pi)).^2 ./ Dsingle);
                Vsingle = 1 ./ Zsingle .* -sinc(kxm .* dslot/(2*pi)) ./ Dsingle;
                % Check normalization, aa should be 1+0j
%                 aa = 1 ./ dx .* sum(Vsingle .* sinc(kxm .* dslot / (2*pi)));
                % Subtract the voltage drop due to the resistance zfeed.
%                 Yl = 1 ./ zfeed_;
%                 Ysingle = 1 ./ Zsingle;
%                 a = 1 ./(-1./dx .* sum(sinc(kxm .* dslot./(2*pi)) ./ Dsingle) .* (Ysingle./(Ysingle+Yl)));
%                 Vsingle = a .* -sinc(kxm .* dslot./(2.*pi)) ./ Dsingle;
    
                Yterms = zeros(Ny_, length(kxm));
                for(ik = 1:length(kxm))
                    Yterms(:, ik) = (Dmy(:, ik) + Ddown(ik) .* mulDdown) .* Vsingle(ik) .* sinc(kxm(ik) .* dslot/(2*pi));
%                     Yterms(:, ik) = (Dmy(:, ik) + Ddown(ik) .* mulDdown) .* Vinf(ik) .* sinc(kxm(ik) .* dslot/(2*pi));
                end
                Yv = zeros(1,Ny_);
                for(ny = 1:Ny_)
                    Yv(ny) = 1 ./ (1./dx .* sum(sinc(kxm .* dslot/(2*pi)).^2)) .* ...
                        (-1 ./ dx .* sum(Yterms(ny,:)));
                end
                
                % Build the toeplitz matrix
                Ymat_(:, :, fi) = toeplitz(real(Yv)) + 1j .* toeplitz(imag(Yv));
                
                send(hDataQueue, nan);
            end
            this.Ymat_fs = [this.Ymat_fs, newfs];
            for(fi = 1:length(newfs))
                this.Ymat(:, :, this.Ymat_fs == newfs(fi)) = Ymat_(:, :, fi);
            end
            this.Ymat_th = th;
            this.Ymat_ph = ph;
            
            function updateWaitbar(~)
                progress = progress + progressStep;
                waitbar(progress/Nf, hWaitbar, ...
                    {sprintf('%.1f%% Calculating Z Matrix...', progress/Nf*100), ...
                     sprintf('%i/%i frequencies done.', progress, Nf)});
            end
            delete(hWaitbar);
            
            dt = toc(tc);
            dispex('Z matrix: Completed in %s.\n', fancyduration(dt));
        end
        
        function Zas = GetInputImpedance(this, fs, th, ph)
            this.InitializeDs(fs, th, ph);
            this.InitializeYMatrix(fs, th, ph);
            
            dispex('Active Z: Calculating for %i frequencies, %i slots.\n', length(fs), this.Ny);
            tc = tic;
            
            Nf = length(fs);
            
            Ymat_ = this.Ymat;
            Ymat_fs_ = this.Ymat_fs;
            
            Ny_ = this.Ny;
            ay0_ = this.ay;
            
            dy = this.unitcell.dy;
            
            
            zfeed_ = this.zfeed;
            
            Zas = zeros(this.Ny, Nf);
            parfor(fi = 1:Nf) % parfor
                f = fs(fi);
                
                [~, ~, ky0, ~] = k(f, 1, th, ph);
                ay_ = ay0_ .* exp(-1j .* ky0 .* (1:Ny_) .* dy);
                
                Y = Ymat_(:, :, Ymat_fs_ == f); %#ok<PFBNS> Broadcast variable
                
                % Add Yl
                Yl = eye(Ny_) ./ zfeed_;
                Yp = Y + Yl;
                
                %% Solve the equivalent circuit and find Zactive
                % Solve for v
                v = Yp \ (ay_.'); % Identical to v = pinv(Yp) * (i.');
                i = Yp \ Y * ay_.';
                % Calculate input impedance
                Zas(:, fi) = v ./ i;
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

            project.NextCommandConditional('Ny>1');
            transform.Reset();
                transform.Name([componentname, '/UnitCell']);
                transform.Vector(0, 'dy', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Ny-1');
            transform.Transform('Shape', 'Translate');

            transform.Reset();
                transform.Name([componentname, '/Stratification']);
                transform.Vector(0, '-dy*2', 0);
                transform.MultipleObjects(0);
                transform.GroupObjects(0);
                transform.Repetitions('1');
            transform.Transform('Shape', 'Translate');
            
            transform.Reset();
                transform.Name([componentname, '/Stratification']);
                transform.Vector(0, 'dy', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Ny+3');
            transform.Transform('Shape', 'Translate');
            
            brick = project.Brick();
            brick.Reset();
                brick.Name('SlotPlane');
                brick.Component(componentname);
                brick.Xrange('-dx/2', 'dx/2');
                brick.Yrange('-5/2*dy', '-dy/2');
                brick.Zrange(0, 0);
                brick.Material('PEC');
            brick.Create();
            
            transform.Reset();
                transform.Name([componentname, ':SlotPlane']);
                transform.Vector(0, '(Ny+2)*dy', 0);
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

            project.NextCommandConditional('Ny>1');
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

    vup = Gxxup.*besselj(0, ky.*wslot./2).*exp(-1j.*ky.*dy.*(nyp-1));
    vdown = Gxxdown.*besselj(0, ky.*wslot./2).*exp(-1j.*ky.*dy.*(nyp-1));
    
    %% Free-space
    fs = FreeSpace();
    Vtm = 1;
    Vte = 1;

    kr = sqrt((kxm).^2 + (ky ).^2);
    isTE = 1;    zteup = fs.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = fs.GetInputImpedance(isTE, f, k0, kr);

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kxm, ky, Vtm, Vte, itmup, iteup);
    Gxxup = Ghm.xx;

    if(~walled)
        isTE = 1;    ztedown = fs.GetInputImpedance(isTE, f, k0, kr);
        isTE = 0;    ztmdown = fs.GetInputImpedance(isTE, f, k0, kr);

        itedown = 1 ./ ztedown;
        itmdown = 1 ./ ztmdown;

        [Ghm] = SpectralGF.hm(z0, k0, kxm, ky, Vtm, Vte, itmdown, itedown);
        Gxxdown = Ghm.xx;
    else
        Gxxdown = 0;
    end
    
    vfsup = Gxxup.*besselj(0, ky.*wslot./2).*exp(-1j.*ky.*dy.*(nyp-1));
    vfsdown = Gxxdown.*besselj(0, ky.*wslot./2).*exp(-1j.*ky.*dy.*(nyp-1));
    %%
    
%     v(isnan(v)) = vfs(isnan(v));
%     v(abs(ky) > 100*k0) = vfs(abs(ky) > 100*k0);
    
    vup(abs(ky) > 1*k0) = vfsup(abs(ky) > 1*k0);
    
    v = vup + vdown - vfsup - vfsdown;
    
%     % Assume anything that is NaN can be approximated by 0.
%     v(isnan(v)) = 0;
end

