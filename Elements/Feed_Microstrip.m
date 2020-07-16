classdef Feed_Microstrip
    properties
        depth
        width
        length
        
        
    end
    properties(SetAccess = protected)
        slot
        dualpol
    end
    methods
        %     \                  /
        %      \                /
        %       \              /
        %        \            /           feed port
        %         \          /                |
        %   end   |          |                V
        %   of  ->|----------|-----------------    -
        %   feed  |          |                |    |
        %         |          |< -- length -- >|  width
        %         |          |                |    |
        %         |----------|-----------------    -
        %         |          |
        %         /          \
        %        /            \
        %       /              \
        %      /                \
        %     /                  \
        
        function this = Feed_Microstrip(slot, depth, width, length)
            this.depth = depth;
            this.width = width;
            this.length = length;
            
            this.slot = slot;
            this.dualpol = contains(lower(class(slot)), 'dualpol');
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
            
            brick = project.Brick();
            solid = project.Solid();
            wcs = project.WCS();
            port = project.Port();
            pick = project.Pick();
            discretefaceport = project.DiscreteFacePort();
            
            s0 = Globals.slot_s0;
            
            project.StoreParameter('feed_depth', this.depth*1e3);
            project.StoreParameter('feed_width', this.width*1e3);
            project.StoreParameter('feed_length', this.length*1e3);
            
            
            if(this.dualpol)
                wcs.Enable();
                wcs.Store('Pre-Feed');
                wcs.MoveWCS('local', [num2str(s0), ' * dx'], [num2str(s0), ' * dy'], 0);
                
                Nports = port.StartPortNumberIteration();
                
                % Delete existing ports.
                solid.Delete([componentname, ':FeedPort2']);
                port.Delete(Nports);
                solid.Delete([componentname, ':FeedPort1']);
                port.Delete(Nports-1);
                
                Nports = Nports - 2;
                
                % X-feed microstrip
                brick.Reset();
                brick.Name('FeedX');
                brick.Component(componentname);
                brick.Material('PEC');
                brick.Xrange('-dx/2-feed_width/2', '-dx/2+feed_width/2');
                brick.Yrange('wfeed/2', '-wfeed/2-feed_length');
                brick.Zrange('-feed_depth', '-feed_depth');
                brick.Create();
                
                % X-feed short
                brick.Reset();
                brick.Name('FeedShortX');
                brick.Component(componentname);
                brick.Material('PEC');
                brick.Xrange('-dx/2-feed_width/2', '-dx/2+feed_width/2');
                brick.Yrange('wfeed/2', 'wfeed/2');
                brick.Zrange('-feed_depth', '0');
                brick.Create();
                
                % X-feed port
                brick.Reset();
                brick.Name('FeedPortX');
                brick.Component(componentname);
                brick.Material('Vacuum');
                brick.Xrange('-dx/2-feed_width/2', '-dx/2+feed_width/2');
                brick.Yrange('-wfeed/2-feed_length', '-wfeed/2-feed_length');
                brick.Zrange('-feed_depth', '0');
                brick.Create();
                
                % Define the port.
                pick.PickEdgeFromId([componentname, ':FeedPortX'], 2, 2);
                pick.PickEdgeFromId([componentname, ':FeedPortX'], 4, 4);
                discretefaceport.Reset();
                discretefaceport.PortNumber(Nports+1);
                discretefaceport.Label('SlotFeedPortX');
                discretefaceport.Impedance('zslot');
                discretefaceport.SetP1(1);
                discretefaceport.SetP2(1);
                discretefaceport.CenterEdge(1);
                discretefaceport.Monitor(1);
                discretefaceport.Create();
                
                % Y-feed microstrip
                brick.Reset();
                brick.Name('FeedY');
                brick.Component(componentname);
                brick.Material('PEC');
                brick.Xrange('wfeed/2', '-wfeed/2-feed_length');
                brick.Yrange('-dy/2-feed_width/2', '-dy/2+feed_width/2');
                brick.Zrange('-feed_depth', '-feed_depth');
                brick.Create();
                
                % Y-feed short
                brick.Reset();
                brick.Name('FeedShortY');
                brick.Component(componentname);
                brick.Material('PEC');
                brick.Xrange('wfeed/2', 'wfeed/2');
                brick.Yrange('-dy/2-feed_width/2', '-dy/2+feed_width/2');
                brick.Zrange('-feed_depth', '0');
                brick.Create();
                
                % Y-feed port
                brick.Reset();
                brick.Name('FeedPortY');
                brick.Component(componentname);
                brick.Material('Vacuum');
                brick.Xrange('-wfeed/2-feed_length', '-wfeed/2-feed_length');
                brick.Yrange('-dy/2-feed_width/2', '-dy/2+feed_width/2');
                brick.Zrange('-feed_depth', '0');
                brick.Create();
                
                % Define the port.
                pick.PickEdgeFromId([componentname, ':FeedPortY'], 3, 3);
                pick.PickEdgeFromId([componentname, ':FeedPortY'], 1, 1);
                discretefaceport.Reset();
                discretefaceport.PortNumber(Nports+2);
                discretefaceport.Label('SlotFeedPortY');
                discretefaceport.Impedance('zslot');
                discretefaceport.SetP1(1);
                discretefaceport.SetP2(1);
                discretefaceport.CenterEdge(1);
                discretefaceport.Monitor(1);
                discretefaceport.Create();
            else
                breakpoint;
            end
            
        end
    end
end