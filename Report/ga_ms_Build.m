ea_CST_Build;
project.StartBulkMode();
wcs = project.WCS();
wcs.Reset();
wcs.Enable();
wcs.MoveWCS('global', '(slot_s0-0.5)*dx', '(slot_s0)*dy', 0);

solid = project.Solid();

solid.Delete('Slot:FeedPort1');
solid.Delete('Slot:FeedPort2');

if(~coax)
%     dms = 'dms';
%     wms = 0.13; % (0.32@0.254=80ohm@20GHz)
%     lms = 'cavity_width/2-slot_feedwidth/2';%+feed_shield_distance-feed_shield_radius-0.1';
%     wms2 = 'feed_ms_width';
%     lms2 = 0.4;%'feed_shield_distance-feed_shield_radius-0.1';
%     
%     patch_angle = 55;
%     patch_l1 = 0.6;
%     patch_capacitance = 0.2e-12;
%     patch_area = 0.7;
    % patch_area = [num2str(patch_capacitance, '%.15g'), '*dms/(2.2*', num2str(Constants.ep0, '%.15g'), ')*1e3'];
    CST.BuildCoaxFeed(project, dms, wms, lms, wms2, lms2, ...
                    [], [], [], ...
                    [], [], [], [], [], ...
                    patch_area, patch_angle, patch_l1);
else
%     dms = 'dms';
%     wms = 0.13; % (0.32@0.254=80ohm@20GHz)
%     lms = 0.5;
%     wms2 = 'feed_ms_width';
%     lms2 = 0.4;%'feed_shield_distance-feed_shield_radius-0.1';
%     core_radius = 0.1125;
%     core_tophole_radius = 'feed_core_radius';
%     core_transition_radius = 0.05;
%     shield_radius = 'feed_core_radius';
%     shield_distance = 0.53;
%     shield_startangle = ['90 - feed_shield_totalangle/2 - feed_ms_buried_angle'];
%     shield_totalangle = 90;
%     shield_Nvias = 3;
%     % cylinder_height = 0.916;
% %     cylinder_Nvias = [];
% %     cylinder_angleoffset = 30;
% %     cylinder_connector_radius = shield_radius + 0.05;
%     patch_angle = 70;
%     patch_l1 = 1;
% %     patch_capacitance = 0.2e-12;
% %     patch_area = [num2str(patch_capacitance, '%.15g'), '*dms/(2.2*', num2str(Constants.ep0, '%.15g'), ')*1e3'];
    CST.BuildCoaxFeed(project, dms, wms, lms, wms2, lms2, ...
                    core_radius, core_tophole_radius, core_transition_radius, ...
                    shield_radius, shield_distance, shield_startangle, shield_totalangle, shield_Nvias, ...
                    patch_area, patch_angle, patch_l1);
end

if(exist('shortcore', 'var') && shortcore)
    project.StoreParameter('feed_core_height', 'hback-feed_ms_depth');
end


port = project.Port();
port.Delete(1);
port.Delete(2);

pick = project.Pick();
if(~coax)
    %% X-feed
    pick.PickEdgeFromId('Feed/Xfeed:Microstrip_Buried', '3', '4');
    pick.PickFaceFromId('Slot:Metal', '10');
    discretefaceport = project.DiscreteFacePort();
    discretefaceport.Reset();
    discretefaceport.PortNumber(1);
    discretefaceport.Label('XFeed');
    discretefaceport.Impedance('slot_impedance');
    discretefaceport.SetP1(1);
    discretefaceport.SetP2(1);
    discretefaceport.CenterEdge(1);
    discretefaceport.Monitor(1);
    discretefaceport.Create();

    %% Y-feed
    pick.PickEdgeFromId('Feed/Yfeed:Microstrip_Buried', '3', '4');
    pick.PickFaceFromId('Slot:Metal', '10');
    discretefaceport = project.DiscreteFacePort();
    discretefaceport.Reset();
    discretefaceport.PortNumber(2);
    discretefaceport.Label('YFeed');
    discretefaceport.Impedance('slot_impedance');
    discretefaceport.SetP1(1);
    discretefaceport.SetP2(1);
    discretefaceport.CenterEdge(1);
    discretefaceport.Monitor(1);
    discretefaceport.Create();
else
    solid.Delete('Vias:Via');
    solid.Delete('Vias:Via_1');
    
    %% Final operations on X-feed.
    if(shortcore)
        solid.Delete('Feed/Xfeed:TopHole');
    else
        solid.Subtract('Slot:Metal', 'Feed/Xfeed:TopHole');
    end
    solid.Subtract('BackingReflector:Metal', 'Feed/Xfeed:GroundHole');

    pick.PickEdgeFromId('BackingReflector:Metal', 6, 6);

    % Define X-waveguide port.
    port.Reset();
    port.PortNumber(2);
    port.Label('Xwaveguide');
    port.NumberOfModes(1);
    port.Coordinates('Picks');
    port.Orientation('negative');
    port.PortOnBound(1);
    port.ClipPickedPortToBound(0)
    port.Create();

    % Rotate the port by 45 degrees to avoid overlap.
    transform.Reset();
    transform.Name('port2 (Xwaveguide)');
    transform.Origin('ShapeCenter');
    transform.Center(0, 0, 0);
    transform.Angle(0, 0, 45);
    transform.MultipleObjects(0);
    transform.GroupObjects(0);
    transform.Repetitions(1);
    transform.MultipleSelection(0);
    transform.Transform('Port', 'Rotate');

    %% Final operations on Y-feed.
    if(shortcore)
        solid.Delete('Feed/Yfeed:TopHole');
    else
        solid.Subtract('Slot:Metal', 'Feed/Yfeed:TopHole');
    end
    solid.Subtract('BackingReflector:Metal', 'Feed/Yfeed:GroundHole');

    pick.PickEdgeFromId('BackingReflector:Metal', 8, 8);

    % Define Y-waveguide port.
    port.Reset();
    port.PortNumber(4);
    port.Label('Ywaveguide');
    port.NumberOfModes(1);
    port.Coordinates('Picks');
    port.Orientation('negative');
    port.PortOnBound(1);
    port.ClipPickedPortToBound(0)
    port.Create();

    % Rotate the port by 45 degrees to avoid overlap.
    transform.Reset();
    transform.Name('port4 (Ywaveguide)');
    transform.Origin('ShapeCenter');
    transform.Center(0, 0, 0);
    transform.Angle(0, 0, 45);
    transform.MultipleObjects(0);
    transform.GroupObjects(0);
    transform.Repetitions(1);
    transform.MultipleSelection(0);
    transform.Transform('Port', 'Rotate');
end

project.EndBulkMode();