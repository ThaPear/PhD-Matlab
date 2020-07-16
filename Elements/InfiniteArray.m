classdef InfiniteArray < handle
    properties
        unitcell  % Unit cell in the infinite array
        
        tlineup   % Stratification above the array
        tlinedown % Stratification below the array
        
        numM     % Number of floquet modes to sum (-numM:numM), default 20 (41 modes)
    end
    methods
        function this = InfiniteArray(unitcell, tlineup, tlinedown)
            this.unitcell = unitcell;
            this.tlineup = tlineup;
            this.tlinedown = tlinedown;
            
            this.numM = 20;
        end
        function Zas = GetInputImpedance(this, fs, th, ph)
            %% Modes
            mx_lin = [-this.numM:this.numM];
            my_lin = mx_lin;
            [mx, my] = meshgrid(mx_lin, my_lin);
            
            dx_ = this.unitcell.dx;
            dy_ = this.unitcell.dy;
            wslot_ = this.unitcell.wslot;
            dslot_ = this.unitcell.dslot;
            walled_ = this.unitcell.walled;
            
            z0 = Constants.z0;
            
%             this = parallel.pool.Constant(this);
            
            Zas = zeros(size(fs));
            for(fi = 1:length(fs))
                f = fs(fi);

                %% Propagation constants.
                [k0, kx0, ky0, ~] = k(f, 1, th, ph);

                kxm = kx0 - 2.*pi.*mx./dx_;
                kxm_lin = kxm(1,:);
                
                %% Solve Slot equivalent circuit.
                Vte = 1;
                Vtm = 1;
                
                % Up stratification
                kym = ky0 - 2.*pi.*my./dy_;
                kr = sqrt((kx0 - 2*pi*mx/dx_).^2 ...
                        + (ky0 - 2*pi*my/dy_).^2);
                isTE = 1;    zupte = this.tlineup.GetInputImpedance(isTE, f, k0, kr);
                isTE = 0;    zuptm = this.tlineup.GetInputImpedance(isTE, f, k0, kr);
                
                iteup = 1 ./ zupte;
                itmup = 1 ./ zuptm;
                
                [Ghm] = SpectralGF.hm(z0, k0, kxm, kym, Vtm, Vte, itmup, iteup);
                GFxxup = Ghm.xx;
                Dinfup = 1./dy_ .* sum(GFxxup .* besselj(0, kym .* wslot_ ./ 2), 1);

                % Down stratification
                if(walled_)
                    kym =     - 2.*pi.*my./dy_;
                    kr = sqrt((kx0 - 2*pi*mx/dx_).^2 ...
                            + (    - 2*pi*my/dy_).^2);
                end
                isTE = 1;    zdownte = this.tlinedown.GetInputImpedance(isTE, f, k0, kr);
                isTE = 0;    zdowntm = this.tlinedown.GetInputImpedance(isTE, f, k0, kr);
                
                itedown = 1 ./ zdownte;
                itmdown = 1 ./ zdowntm;
                
                [Ghm] = SpectralGF.hm(z0, k0, kxm, kym, Vtm, Vte, itmdown, itedown);
                GFxxdown = Ghm.xx;
                Dinfdown = 1./dy_ .* sum(GFxxdown .* besselj(0, kym .* wslot_ ./ 2), 1);
                
                % Impedance
                Za   =  -1./dx_ .* sum(sinc(kxm_lin .* dslot_ ./ (2.*pi)).^2 ./ (Dinfup + Dinfdown));
                Zas(fi) = Za;
            end
        end
        function [krhos, k0, Dinfups, Dinfdowns, Iteups, Itmups, Itedowns, Itmdowns] = GetDinf(this, f)
            %             f = 10.6e9;
            th = eps*pi/180;
            ph = eps*pi/180;
            [k0, kx0, ky0, ~] = k(f, 1, th, ph);
%             ky0 = 0;
            
            kxs = (0:0.001:2)*k0;
%             kxs = (1:0.000001:1.3)*k0;
            
            dx_ = this.unitcell.dx;
            dy_ = this.unitcell.dy;
            wslot_ = this.unitcell.wslot;
            dslot_ = this.unitcell.dslot;
            walled_ = this.unitcell.walled;
            
            Dinfups = zeros(size(kxs));
            Dinfdowns = zeros(size(kxs));
            Iteups = zeros(size(kxs));
            Itmups = zeros(size(kxs));
            Itedowns = zeros(size(kxs));
            Itmdowns = zeros(size(kxs));
            
            
            mx_lin = [-this.numM:this.numM];
            my_lin = mx_lin;
            [mx, my] = meshgrid(mx_lin, my_lin);
            my = my(:,1);
            
            z0 = Constants.z0;
            tic;
            parfor(kxi = 1:length(kxs))
                kx = kxs(kxi);
                % Up stratification

                kxm = kx;
                kxm_lin = kxm(1,:);
                
                %% Solve Slot equivalent circuit.
                Vte = 1;
                Vtm = 1;
                
                % Up stratification
                kym = ky0 - 2.*pi.*my./dy_;
                kr = sqrt((kx).^2 ...
                        + (ky0 - 2*pi*my/dy_).^2);
                isTE = 1;    zupte = this.tlineup.GetInputImpedance(isTE, f, k0, kr);
                isTE = 0;    zuptm = this.tlineup.GetInputImpedance(isTE, f, k0, kr);
                
                iteup = 1 ./ zupte;
                itmup = 1 ./ zuptm;
                
                [Ghm] = SpectralGF.hm(z0, k0, kxm, kym, Vtm, Vte, itmup, iteup);
                GFxxup = Ghm.xx;
                Dinfup = -1./dy_ .* sum(GFxxup .* besselj(0, kym .* wslot_ ./ 2), 1);
                
                % Down
                if(walled_)
                    kym =     - 2.*pi.*my./dy_;
                    kr = sqrt((kx).^2 ...
                            + (    - 2*pi*my/dy_).^2);
                end
                
                isTE = 1;    zdownte = this.tlinedown.GetInputImpedance(isTE, f, k0, kr);
                isTE = 0;    zdowntm = this.tlinedown.GetInputImpedance(isTE, f, k0, kr);
                
                itedown = 1 ./ zdownte;
                itmdown = 1 ./ zdowntm;
                
                [Ghm] = SpectralGF.hm(z0, k0, kxm, kym, Vtm, Vte, itmdown, itedown);
                GFxxdown = Ghm.xx;
                Dinfdown = -1./dy_ .* sum(GFxxdown .* besselj(0, kym .* wslot_ ./ 2), 1);
                
                Dinfups(kxi) = Dinfup;
                Dinfdowns(kxi) = Dinfdown;
                Iteups(kxi) = iteup(my == 0);
                Itedowns(kxi) = itedown(my == 0);
                Itmups(kxi) = itmup(my == 0);
                Itmdowns(kxi) = itmdown(my == 0);
                Itmups(kxi) = GFxxup(my == 0);
                Itmdowns(kxi) = GFxxdown(my == 0);
                
                
            end
            toc;
            krhos = sqrt((kxs).^2 ...
                       + (ky0).^2);
%             figureex(100); plot(kxs/k0, abs(iteups));
%             figureex(101); plot(kxs/k0, abs(itmups));
%             figureex(102); plot(kxs/k0, abs(itedowns));
%             figureex(103); plot(kxs/k0, abs(itmdowns));
%             figureex(104); plot(kxs/k0, abs(Dinfups)); title('dinfup');
%             figureex(105); plot(kxs/k0, abs(Dinfdowns)); title('dinfdown');
%             figureex(106); plot(kxs/k0, abs(Dinfdowns+Dinfups)); title('sum');
        end
        function BuildCST(this, project, parentcomponent)
            this.unitcell.BuildCST(parentcomponent);
            
            % Set boundary conditions to periodic.
            boundary = project.Boundary();
            boundary.Xmin('periodic');  boundary.Xmax('periodic');
            boundary.Ymin('periodic');  boundary.Ymax('periodic');
            
            wcs.Enable();
            wcs.Store('Pre-Slot');
            wcs.RotateWCS('u', 180);
            % Build down-stratification.
            this.tlinedown.BuildCST(project, parentcomponent);
            
            wcs.Restore('Pre-Slot');
            wcs.Delete('Pre-Slot');
            
            % Build up-stratification.
            this.tlineup.BuildCST(project, parentcomponent);
            
            wcs.Disable();
        end
    end
end