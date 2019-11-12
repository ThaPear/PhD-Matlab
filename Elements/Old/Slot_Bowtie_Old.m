classdef Slot
    properties
        dx              % Unit cell dimensions
        dy              % Unit cell dimensions
        wslot           % Width of slot
        dslot           % Length of feed
        
        numM            % Number of floquet modes to sum (-numM:numM), default 10
        
        tlineup         % Stratification above the slot.
        tlinedown       % Stratification below the slot.
        
        walled          % Vertical walls in y-direction?
        
        bowtieratio     % Ratio of wslot and the width of the bowtie
    end
    methods
        function this = Slot(dx, dy, wslot, dslot, tlineup, tlinedown, walled, bowtieratio)
            %            <------------ p ------------>
            %           /-------------------------------/
            %          /         /           /         /|
            %         /         /           /         /||
            %        /         /-----------/  -      /|||
            %       /         /           /  /      /||||
            %      /         /           /  /dslot /||||/
            %     /         /           /  /      /||||/
            %    /         /-----------/  -      /||||/
            %   /         /           /         /|||<--- if walled
            %  /         /<- wslot ->/         /||||/
            % /-------------------------------/||||/
            % ||||/                           ||||/
            % |||/                            |||/
            % ||/                             ||/
            % |/                              |/
            % /                               /
            
            % |-------------------------------|
            % |         |           |         |
            % |         |           |         |
            % |         |           |         |
            % |          \         /          |
            % |           \       /           |
            % |            |     | -          |
            % |            |     | |          |
            % |            |<-*->| dbowtie    | *wbowtie
            % |            |     | |          |
            % |            |     | -          |
            % |           /       \           |
            % |          /         \          |
            % |         |           |         |
            % |         |           |         |
            % |         |<- wslot ->|         |
            % |-------------------------------|
            
            this.dx = dx;
            this.dy = dy;
            this.wslot = wslot;
            this.dslot = dslot;
            this.tlineup = tlineup;
            this.tlinedown = tlinedown;
            this.walled = walled;
            if(nargin < 8)
                bowtieratio = 1;
            end
            this.bowtieratio = bowtieratio;
            
            this.numM = 10;
        end
        
        function Zas = GetInputImpedance(this, fs, th, ph)
            %% Modes
            mx_lin = [-this.numM:this.numM];
            my_lin = mx_lin;
            [mx, my] = meshgrid(mx_lin, my_lin);
            
            z0 = Constants.z0;
            
            Zas = zeros(size(fs));
            parfor(fi = 1:length(fs))
                f = fs(fi);

                %% Propagation constants.
                [k0, kx0, ky0, ~] = k(f, 1, th, ph);

                kxm = kx0 - 2.*pi.*mx./this.dx;
                kxm_lin = kxm(1,:);
                
                %% Solve Slot equivalent circuit.
                Vte = 1;
                Vtm = 1;
                
                % Up stratification
                kym = ky0 - 2.*pi.*my./this.dy;
                kr = sqrt((kx0 - 2*pi*mx/this.dx).^2 ...
                        + (ky0 - 2*pi*my/this.dy).^2);
                isTE = 1;    zupte = this.tlineup.GetInputImpedance(isTE, f, k0, kr);
                isTE = 0;    zuptm = this.tlineup.GetInputImpedance(isTE, f, k0, kr);
                
                iteup = 1 ./ zupte;
                itmup = 1 ./ zuptm;
                
                [Ghm] = SpectralGF.hm(z0, k0, kxm, kym, Vtm, Vte, itmup, iteup);
                GFxxup = Ghm.xx;
                Dinfup = -1./this.dy .* sum(GFxxup .* besselj(0, kym .* this.wslot ./ 2), 1);

                % Down stratification
                if(this.walled)
                    kym =     - 2.*pi.*my./this.dy;
                    kr = sqrt((kx0 - 2*pi*mx/this.dx).^2 ...
                            + (    - 2*pi*my/this.dy).^2);
                end
                isTE = 1;    zdownte = this.tlinedown.GetInputImpedance(isTE, f, k0, kr);
                isTE = 0;    zdowntm = this.tlinedown.GetInputImpedance(isTE, f, k0, kr);
                
                itedown = 1 ./ zdownte;
                itmdown = 1 ./ zdowntm;
                
                [Ghm] = SpectralGF.hm(z0, k0, kxm, kym, Vtm, Vte, itmdown, itedown);
                GFxxdown = Ghm.xx;
                Dinfdown = -1./this.dy .* sum(GFxxdown .* besselj(0, kym .* this.wslot ./ 2), 1);
                
                % Impedance
                Za   =  1./this.dx .* sum(sinc(kxm_lin .* this.dslot ./ (2.*pi)).^2 ./ (Dinfup + Dinfdown));
                Zas(fi) = Za;
            end
        end
        
        function h = GetHeight(~)
            h = 0;
        end
        
        function BuildCST(this, project)
            lambda0 = str2double(project.RestoreParameter('lambda0'))/1e3;
            
            project.StoreParameter('wslot', [num2str(this.wslot/lambda0, '%.10g'), '*lambda0']);
            project.StoreParameter('dslot', [num2str(this.dslot/lambda0, '%.10g'), '*lambda0']);
            project.StoreParameter('dx', [num2str(this.dx/lambda0, '%.10g'), '*lambda0']);
            project.StoreParameter('dy', [num2str(this.dy/lambda0, '%.10g'), '*lambda0']);
            project.MakeSureParameterExists('zslot', 50);
            project.MakeSureParameterExists('bowtieratio', this.bowtieratio);
            
            component = project.Component();
            polygon3d = project.Polygon3D();
            transform = project.Transform();
            brick = project.Brick();
            port = project.Port();
            discretefaceport = project.DiscreteFacePort();
            wcs = project.WCS();
            extrudecurve = project.ExtrudeCurve();
            
            component.New('Slot');
            
            % Draw the slot metal.
            polygon3d.Reset();
            polygon3d.Curve('Slot');
            polygon3d.Name('Metal');
            polygon3d.Point('-dx/2'   , '-dy/2'   , '0');
            polygon3d.Point( 'dx/2'   , '-dy/2'   , '0');
            polygon3d.Point( 'dx/2'   , '-wslot/2', '0');
            polygon3d.Point( 'dslot/2', '-wslot/2', '0');
            polygon3d.Point( 'dslot/2/bowtieratio', '-wslot/2/bowtieratio', '0');
            polygon3d.Point('-dslot/2/bowtieratio', '-wslot/2/bowtieratio', '0');
            polygon3d.Point('-dslot/2', '-wslot/2', '0');
            polygon3d.Point('-dx/2'   , '-wslot/2', '0');
            polygon3d.Point('-dx/2'   , '-dy/2'   , '0');
            polygon3d.Create();
            
            extrudecurve.Reset();
            extrudecurve.Curve('Slot');
            extrudecurve.Component('Slot');
            extrudecurve.Name('Metal');
            extrudecurve.Material('PEC');
            extrudecurve.Thickness(0);
            extrudecurve.DeleteProfile(1);
            extrudecurve.Create();
            
            % Mirror the metal in x = 0.
            transform.Reset();
            transform.Name('Slot:Metal');
            transform.PlaneNormal(0, 1, 0);
            transform.Center(0, 0, 0);
            transform.Origin('Free');
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Transform('Shape', 'Mirror');
            
%             % Create one side of the slot's metal.
%             brick.Reset();
%             brick.Component('Slot');
%             brick.Name('Metal');
%             brick.Xrange('-dx/2', 'dx/2');
%             brick.Yrange('-dy/2', '-wslot/2');
%             brick.Zrange(0, 0);
%             brick.Material('PEC');
%             brick.Create();
%             
%             % Create the other side of the slot's metal.
%             brick.Name('Metal2');
%             brick.Yrange('wslot/2', 'dy/2');
%             brick.Create();
%             
%             % Merge the two metal plates.
%             solid.Add('Slot:Metal', 'Slot:Metal2');
            
            % Create the vacuum plate for the feed.
            
            brick.Reset();
            brick.Component('Slot');
            brick.Name('FeedPort');
            brick.Xrange('-dslot/2/bowtieratio', 'dslot/2/bowtieratio');
            brick.Yrange('-wslot/2/bowtieratio', 'wslot/2/bowtieratio');
            brick.Zrange(0, 0);
            brick.Material('Vacuum');
            brick.Create();
            
            pick = project.Pick();
            pick.PickEdgeFromId('Slot:FeedPort', 2, 2);
            pick.PickEdgeFromId('Slot:FeedPort', 4, 4);
            
            % Define the port.
            Nports = port.StartPortNumberIteration();
            portnumber = Nports + 1;
            discretefaceport.Reset();
            discretefaceport.PortNumber(portnumber);
            discretefaceport.Label('SlotFeed');
            discretefaceport.Impedance('zslot');
            discretefaceport.SetP1(1);
            discretefaceport.SetP2(1);
            discretefaceport.CenterEdge(1);
            discretefaceport.Monitor(1);
            discretefaceport.Create();
            
%             obj.tlineup.BuildCST(project);
            wcs.Enable();
            wcs.Store('Pre-Slot');
            wcs.RotateWCS('u', 180);
            
            % Build walls if needed.
            if(this.walled)
                project.MakeSureParameterExists('dwall', 0);
                brick.Reset();
                brick.Name('Walls');
                brick.Component('Slot');
                brick.Xrange('-dx/2',  'dx/2');
                brick.Yrange('-dy/2', '-dy/2+dwall');
                brick.Zrange('0', this.tlinedown.GetHeight()*1e3);
                brick.Material('PEC');
                brick.Create();

                transform.Reset();
                transform.Name(['Slot:Walls']);
                transform.Material('PEC');
                transform.Vector(0, 'dy-dwall', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Transform('Shape', 'Translate');
            end
            % Build down-stratification.
            this.tlinedown.BuildCST(project);
            
            wcs.Restore('Pre-Slot');
            wcs.Delete('Pre-Slot');
            
            % Build up-stratification.
            this.tlineup.BuildCST(project);
            
            wcs.Disable();
        end
        
    end
end