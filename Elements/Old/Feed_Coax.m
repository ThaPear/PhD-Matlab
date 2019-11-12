classdef Feed_Coax
    properties
        
    end
    methods
        function this = Feed_Coax()
            
        end
        function BuildCST(this, project)
            cylinder = project.Cylinder();
            component = project.Component();
            brick = project.Brick();
            solid = project.Solid();
            port = project.Port();
            pick = project.Pick();
            discretefaceport = project.DiscreteFacePort();
            polygon3d = project.Polygon3D();
            blendcurve = project.BlendCurve();
            tracefromcurve = project.TraceFromCurve();
            transform = project.Transform();
            
            lambda0 = Constants.c0/11e9;
            
            shield_ncyl = 5;
            shield_angles = linspace(0, 180, shield_ncyl);
            ratio = Coaxial.GetRadiusRatio(50, 3.66);
            
%             lms = 3.58e-3 - core_radius - core_ms_radius;

            project.MakeSureParameterExists('dms', '0.25');
%             project.MakeSureParameterExists('wms', 'dms * 28/25'); % 80 ohm
%             project.MakeSureParameterExists('zms', 50);
%             project.MakeSureParameterExists('res_w', 'dms * 39/25'); % 63 Ohm
%             
            project.StoreParameter('feed_core_height', 'hback-dms');
            project.StoreParameter('feed_core_radius', ['max(hback/12, 0.25/2)']);
%             project.StoreParameter('feed_core_ms_radius', 0.1);%'res_w/2 / sin(atn(res_w/2 / core_radius)) - core_radius + 0.001');
% %             project.StoreParameter('core_ms_radius', 'max(0.001, wms/2 - core_radius)');
%             project.StoreParameter('feed_core_tophole_inner', 'core_radius+0.1');
%             project.StoreParameter('feed_core_tophole_outer', 'core_radius+0.3');
% 
%             
            project.StoreParameter('feed_shield_radius', ['max(hback/12, 0.25/2)']);
            project.StoreParameter('feed_shield_height', 'hback');
            project.StoreParameter('feed_shield_distance', [num2str(ratio), ' * feed_core_radius + feed_shield_radius + 0.1']);
%             
            project.StoreParameter('feed_shield_Nvias', 5);
            project.StoreParameter('feed_shield_startangle', 90);
            project.StoreParameter('feed_shield_arc', '180');
%             
%             project.StoreParameter('feed_lstraight', 1);
%             project.StoreParameter('feed_l2', 'feed_ltotal - feed_lstraight - feed_core_height');
%             project.StoreParameter('feed_ltotal', '0.25 * lambda0 / sqr(3.66) * 9.5/5.5');
%             project.StoreParameter('feed_angle', 45);
            
%             polygon3d.Reset();
%             polygon3d.Name('FeedMS');
%             polygon3d.Curve('FeedMS');
%             polygon3d.Point(
%             polygon3d.Point(
%             polygon3d.Point(
%             polygon3d.Point(
%             polygon3d.Point(
%             polygon3d.Point(
            
            project.StoreParameter('via_dx', -2);
            project.StoreParameter('via_dy', -2);
            project.StoreParameter('shieldangle', 45);
            
            s0 = Globals.slot_s0;
            
            cylinder.Reset();
            cylinder.Component('Vias');
            cylinder.Name('Core');
            cylinder.Axis('z');
            cylinder.Material('PEC');
            cylinder.Xcenter([num2str(s0) '*dx+via_dx']);
            cylinder.Ycenter([num2str(s0) '*dy-dy/2+via_dy']);
            cylinder.Zrange('-hback', '-hback+feed_core_height');
            cylinder.Innerradius(0);
            cylinder.Outerradius('feed_core_radius');
            cylinder.Create();
            
            cylinder.Reset();
            cylinder.Component('Vias');
            cylinder.Name('Shield');
            cylinder.Axis('z');
            cylinder.Material('PEC');
            cylinder.Xcenter([num2str(s0) '*dx+via_dx+feed_shield_distance']);
            cylinder.Ycenter([num2str(s0) '*dy-dy/2+via_dy']);
            cylinder.Zrange('-hback', '-hback+feed_shield_height');
            cylinder.Innerradius(0);
            cylinder.Outerradius('feed_shield_radius');
            cylinder.Create();
            
            transform.Reset();
            transform.Name('Vias:Shield');
            transform.Origin('Free');
            transform.Center([num2str(Globals.slot_s0) '*dx+via_dx'], [num2str(Globals.slot_s0) '*dy-dy/2+via_dy'], 0);
            transform.Angle(0, 0, 'feed_shield_startangle');
            transform.Transform('Shape', 'Rotate');
            
            transform.Reset();
            transform.Name('Vias:Shield');
            transform.Origin('Free');
            transform.Center([num2str(Globals.slot_s0) '*dx+via_dx'], [num2str(Globals.slot_s0) '*dy-dy/2+via_dy'], 0);
            transform.Angle(0, 0, 'feed_shield_arc / (feed_shield_Nvias-1)');
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Repetitions('feed_shield_Nvias-1')
            transform.Transform('Shape', 'Rotate');
            
            transform.Reset();
            transform.Name('Vias:Core');
            transform.AddName('Vias:Shield');
            transform.Origin('Free');
            transform.Center(0,0,0);
            transform.PlaneNormal(1, -1, 0);
            transform.MultipleObjects(1);
            transform.Transform('Shape', 'Mirror');
            
            
            
            
            
            
        end
    end
end
    