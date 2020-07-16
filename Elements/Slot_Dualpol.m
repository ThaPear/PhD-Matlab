classdef Slot_Dualpol < Slot
    properties
        wfeed
    end
    methods
        function this = Slot_Dualpol(dx, dy, wslot, dslot, tlineup, tlinedown, wfeed)
            this@Slot(dx, dy, wslot, dslot, tlineup, tlinedown, 0, 1);
            this.wfeed = wfeed;
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
            
            project.StoreParameter('wslot', this.wslot*1e3);
            project.StoreParameter('lfeed', this.dslot*1e3/2);
            project.StoreParameter('wfeed', ['wslot / 2']);
            project.StoreParameter('dx', this.dx*1e3);
            project.StoreParameter('dy', this.dy*1e3);
            project.MakeSureParameterExists('zslot', 80);
            
            component = project.Component();
            polygon3d = project.Polygon3D();
            transform = project.Transform();
            solid = project.Solid();
            brick = project.Brick();
            port = project.Port();
            discretefaceport = project.DiscreteFacePort();
            wcs = project.WCS();
            extrudecurve = project.ExtrudeCurve();
            
            if(Globals.exists('slot_s0'))
                s0 = Globals.slot_s0;
            else
                s0 = 0;
            end
            
            project.MakeSureParameterExists('slot_s0', s0);
            
            wcs.Enable();
            wcs.Store('Pre-Slot');
            wcs.MoveWCS('local', 'slot_s0 * dx', 'slot_s0 * dx', 0);
            
            component.New(componentname);
            
            % Draw the slot metal.
            polygon3d.Reset();
            polygon3d.Curve('Slot');
            polygon3d.Name('Metal');
            polygon3d.Point('-wslot/2', '-wslot/2', '0');
            polygon3d.Point('-wfeed/2', '-dy/2+lfeed/2', '0');
            polygon3d.Point('-wfeed/2', '-dy/2-lfeed/2', '0');
            polygon3d.Point('-wslot/2', '-dy+wslot/2', '0');
            polygon3d.Point('-dx/2+lfeed/2', '-dy+wfeed/2', '0');
            polygon3d.Point('-dx/2-lfeed/2', '-dy+wfeed/2', '0');
            polygon3d.Point('-dx+wslot/2', '-dy+wslot/2', '0');
            polygon3d.Point('-dx+wfeed/2', '-dy/2-lfeed/2', '0');
            polygon3d.Point('-dx+wfeed/2', '-dy/2+lfeed/2', '0');
            polygon3d.Point('-dx+wslot/2', '-wslot/2', '0');
            polygon3d.Point('-dx/2-lfeed/2', '-wfeed/2', '0');
            polygon3d.Point('-dx/2+lfeed/2', '-wfeed/2', '0');
            polygon3d.Point('-wslot/2', '-wslot/2', '0');
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
            transform.Vector('0', 'dy', '0');
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Transform('Shape', 'Translate');
            
            transform.Vector('dx', '0', '0');
            transform.Transform('Shape', 'Translate');
            
            % Cut the metal to unit cell size.
            brick.Reset();
            brick.Component(componentname);
            brick.Name('UnitCell');
            brick.Xrange(['dx * (-1/2 - slot_s0)'], ['dx * (1/2 - slot_s0)']);
            brick.Yrange(['dy * (-1/2 - slot_s0)'], ['dy * (1/2 - slot_s0)']);
            brick.Zrange(-0.1, 0.1);
            brick.Material('Transparent');
            brick.Create();
            
            solid.Intersect([componentname, ':Metal'], [componentname, ':UnitCell']);
            
            
            % Make feed vacuum.
            brick.Reset();
            brick.Component(componentname);
            brick.Name('FeedPort1');
            brick.Xrange('-dx/2-lfeed/2', '-dx/2+lfeed/2');
            brick.Yrange('-wfeed/2', 'wfeed/2');
            brick.Zrange(0, 0);
            brick.Material('Vacuum');
            brick.Create();
            
            % Other feed vacuum
            brick.Reset();
            brick.Component(componentname);
            brick.Name('FeedPort2');
            brick.Xrange('-wfeed/2', 'wfeed/2');
            brick.Yrange('-dy/2-lfeed/2', '-dy/2+lfeed/2');
            brick.Zrange(0, 0);
            brick.Material('Vacuum');
            brick.Create();
            
            Nports = port.StartPortNumberIteration();
            
            % Define the port.
            pick = project.Pick();
            pick.PickEdgeFromId([componentname, ':FeedPort1'], 2, 2);
            pick.PickEdgeFromId([componentname, ':FeedPort1'], 4, 4);
            discretefaceport.Reset();
            discretefaceport.PortNumber(Nports+1);
            discretefaceport.Label('SlotFeedPort1');
            discretefaceport.Impedance('zslot');
            discretefaceport.SetP1(1);
            discretefaceport.SetP2(1);
            discretefaceport.CenterEdge(1);
            discretefaceport.Monitor(1);
            discretefaceport.Create();
            
            % Define the port.
            pick = project.Pick();
            pick.PickEdgeFromId([componentname, ':FeedPort2'], 1, 1);
            pick.PickEdgeFromId([componentname, ':FeedPort2'], 3, 3);
            discretefaceport.Reset();
            discretefaceport.PortNumber(Nports+2);
            discretefaceport.Label('SlotFeedPort2');
            discretefaceport.Impedance('zslot');
            discretefaceport.SetP1(1);
            discretefaceport.SetP2(1);
            discretefaceport.CenterEdge(1);
            discretefaceport.Monitor(1);
            discretefaceport.Create();
            
            
            wcs.Restore('Pre-Slot');
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