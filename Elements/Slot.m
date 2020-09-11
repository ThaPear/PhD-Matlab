classdef Slot
    properties
        dx              % Unit cell dimensions
        dy              % Unit cell dimensions
        wslot           % Width of slot
        dslot           % Length of feed
        
        walled          % Vertical walls in y-direction?
    end
    methods
        function this = Slot(dx, dy, wslot, dslot, walled)
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
            
            if(nargin < 5)
                walled = 1;
            end
            
            this.dx = dx;
            this.dy = dy;
            this.wslot = wslot;
            this.dslot = dslot;
            this.walled = walled;
        end
        function h = GetHeight(~)
            h = 0;
        end
        function BuildCST(this, project, parentcomponent)
            if(nargin < 3 || isempty(parentcomponent))
                parentcomponent = '';
            else
                if(~strcmp(parentcomponent(end), '/'))
                    parentcomponent = [parentcomponent, '/'];
                end
            end
            componentname = [parentcomponent, 'Slot'];
            
            project.StoreParameter('slot_width', this.wslot*1e3);
            project.StoreParameter('slot_feedlength', this.dslot*1e3);
            project.StoreParameter('slot_bowtie', ['slot_feedlength']);
            project.StoreParameter('slot_bowtie_outer', ['slot_bowtie']);
            project.StoreParameter('slot_feedwidth', ['slot_width']);
            project.StoreParameter('dx', this.dx*1e3);
            project.StoreParameter('dy', this.dy*1e3);
            project.MakeSureParameterExists('slot_impedance', 80);
            
            component = project.Component();
            polygon3d = project.Polygon3D();
            transform = project.Transform();
            brick = project.Brick();
            port = project.Port();
            discretefaceport = project.DiscreteFacePort();
            wcs = project.WCS();
            extrudecurve = project.ExtrudeCurve();
            
            wcs.Enable();
            
            component.New(componentname);
            
            % Draw the slot metal.
            polygon3d.Reset();
            polygon3d.Curve('Slot');
            polygon3d.Name('Metal');
            polygon3d.Point('-dx/2'   , '-dy/2'   , '0');
            polygon3d.Point( 'dx/2'   , '-dy/2'   , '0');
            polygon3d.Point( 'dx/2'   , '-slot_width/2', '0');
            polygon3d.Point( 'slot_bowtie_outer/2', '-slot_width/2', '0');
            polygon3d.Point( 'slot_bowtie/2', '-slot_feedwidth/2', '0');
            polygon3d.Point('-slot_bowtie/2', '-slot_feedwidth/2', '0');
            polygon3d.Point('-slot_bowtie_outer/2', '-slot_width/2', '0');
            polygon3d.Point('-dx/2'   , '-slot_width/2', '0');
            polygon3d.Point('-dx/2'   , '-dy/2'   , '0');
            polygon3d.Create();
            
            extrudecurve.Reset();
            extrudecurve.Curve('Slot');
            extrudecurve.Component(componentname);
            extrudecurve.Name('Metal');
            extrudecurve.Material('PEC');
            extrudecurve.Thickness(0);
            extrudecurve.DeleteProfile(1);
            extrudecurve.Create();
            
            % Mirror the metal in x = 0.
            transform.Reset();
            transform.Name([componentname, ':Metal']);
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
%             brick.Yrange('-dy/2', '-slot_width/2');
%             brick.Zrange(0, 0);
%             brick.Material('PEC');
%             brick.Create();
%             
%             % Create the other side of the slot's metal.
%             brick.Name('Metal2');
%             brick.Yrange('slot_width/2', 'dy/2');
%             brick.Create();
%             
%             % Merge the two metal plates.
%             solid.Add('Slot:Metal', 'Slot:Metal2');
            
            % Create the vacuum plate for the feed.
            
            brick.Reset();
            brick.Component(componentname);
            brick.Name('FeedPort');
            brick.Xrange('-slot_feedlength/2', 'slot_feedlength/2');
            brick.Yrange('-slot_feedwidth/2', 'slot_feedwidth/2');
            brick.Zrange(0, 0);
            brick.Material('Vacuum');
            brick.Create();
            
            pick = project.Pick();
            pick.PickEdgeFromId([componentname, ':FeedPort'], 2, 2);
            pick.PickEdgeFromId([componentname, ':FeedPort'], 4, 4);
            
            % Define the port.
            Nports = port.StartPortNumberIteration();
            portnumber = Nports + 1;
            discretefaceport.Reset();
            discretefaceport.PortNumber(portnumber);
            discretefaceport.Label('SlotFeed');
            discretefaceport.Impedance('slot_impedance');
            discretefaceport.SetP1(1);
            discretefaceport.SetP2(1);
            discretefaceport.CenterEdge(1);
            discretefaceport.Monitor(1);
            discretefaceport.Create();
            
            % Build walls if needed.
            if(this.walled)
                project.MakeSureParameterExists('dwall', 0);
                project.MakeSureParameterExists('hback', 1);
                brick.Reset();
                brick.Name('Walls');
                brick.Component(componentname);
                brick.Xrange('-dx/2',  'dx/2');
                brick.Yrange('-dy/2', '-dy/2+dwall');
                brick.Zrange('0', '-hback');
                brick.Material('PEC');
                brick.Create();

                transform.Reset();
                transform.Name([componentname, ':Walls']);
                transform.Material('PEC');
                transform.Vector(0, 'dy-dwall', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Transform('Shape', 'Translate');
            end
        end
        
    end
end