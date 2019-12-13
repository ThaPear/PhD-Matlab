% close all;
SetupPath;
clear;
clear global;

newproject = 0;
ok = questdlg('Build in new project?','','Active project','New project','New project');
if(strcmp(ok, 'New project'))
    newproject = 1;
end

if(newproject)
    [project, dsproject] = CST.InitializePeriodicProject();
    project.StoreParameter('dx', 4.35);
    project.StoreParameter('dy', 'dx');
    project.StoreParameter('dms', 0.127);
    project.StoreParameter('slot_feedwidth', 0.1);
    project.StoreParameter('slot_feedlength', 0.1);
    project.StoreParameter('hback', 0.9);
    project.StoreParameter('slot_s0', 0.3);
    project.StoreParameter('zslot', 80);
    project.StoreParameter('openboundary_distance', 1);
else
    project = CST.Application.Active3D();
    dsproject = CST.Application.ActiveDS();
end

brick = project.Brick();
material = project.Material();
solid = project.Solid();
wcs = project.WCS();
pick = project.Pick();
discretefaceport = project.DiscreteFacePort();
port = project.Port();
transform = project.Transform();

%% Construct the unit cell.
if(newproject)
    brick.Reset();
    brick.Component('UnitCell');
    brick.Name('SlotPlane');
    brick.Material('PEC');
    brick.Xrange('-dx/2', 'dx/2');
    brick.Yrange('-dy/2', 'dy/2');
    brick.Zrange(0,0);
    brick.Create();

    brick.Reset();
    brick.Component('UnitCell');
    brick.Name('GroundPlane');
    brick.Material('PEC');
    brick.Xrange('-dx/2', 'dx/2');
    brick.Yrange('-dy/2', 'dy/2');
    brick.Zrange('-hback', '-hback');
    brick.Create();
    % Dielectric material.
    material.Reset();
    materialname = '2.2';
    material.Name(materialname);
    material.Folder('Generated');
    material.Colour(0.5, 0.5, 0.5);
    material.Epsilon(2.2);
    material.Transparency(0.5);
    material.Create();

    brick.Reset();
    brick.Component('UnitCell');
    brick.Name('Dielectric');
    brick.Material('Generated/2.2');
    brick.Xrange('-dx/2', 'dx/2');
    brick.Yrange('-dy/2', 'dy/2');
    brick.Zrange('-hback', 0);
    brick.Create();

    cavity = Cavity_Dualpol(1.65e-3, 3.65e-3);
    cavity.BuildCST(project);
    vias = Vias_Dualpol(1.3e-3, 0.35e-3);
    vias.BuildCST(project);

    solid.Subtract('UnitCell:Dielectric', 'Cavity:Cavity');
    solid.Subtract('UnitCell:Dielectric', 'Cavity:CavityDiag');
end

%% Construct the feed.
wcs.Reset();
wcs.Enable();
wcs.MoveWCS('global', '(slot_s0-0.5)*dx', '(slot_s0)*dy', 0);

dms = 'dms';
wms = 0.13; % (0.32@0.254=80ohm@20GHz)
lms = 'cavity_width/2-slot_feedwidth/2';%+feed_shield_distance-feed_shield_radius-0.1';
wms2 = 'feed_ms_width';
lms2 = 0.4;%'feed_shield_distance-feed_shield_radius-0.1';
core_radius = 0.1125;
core_tophole_radius = 'feed_core_radius';
core_transition_radius = 0.05;
shield_radius = 'feed_core_radius';
shield_distance = 0.6;
shield_startangle = ['90 - feed_shield_totalangle/2 - feed_ms_buried_angle'];
shield_totalangle = 70;
shield_Nvias = 3;
% cylinder_height = 0.916;
cylinder_Nvias = [];
cylinder_angleoffset = 30;
cylinder_connector_radius = shield_radius + 0.05;
patch_angle = 55;
patch_l1 = 0.6;
patch_capacitance = 0.2e-12;

patch_area = [num2str(patch_capacitance, '%.15g'), '*dms/(2.2*', num2str(Constants.ep0, '%.15g'), ')*1e3'];

CST.BuildCoaxFeed(project, dms, wms, lms, wms2, lms2, ...
                core_radius, core_tophole_radius, core_transition_radius, ...
                shield_radius, shield_distance, shield_startangle, shield_totalangle, shield_Nvias, ...
                patch_area, patch_angle, patch_l1);

%% Final operations on X-feed.
if(newproject)
    solid.Subtract('UnitCell:SlotPlane', 'Feed/Xfeed:TopHole');
    solid.Subtract('UnitCell:GroundPlane', 'Feed/Xfeed:GroundHole');
    if(isnumeric(patch_area) && patch_area == 0)
        solid.Delete('Feed/Xfeed:Short');
    else
        solid.Delete('Feed/Xfeed:Patch');
    end
    solid.Delete('Feed/Xfeed:SlotFeed');

    % Define X-microstrip port.
    pick.PickEdgeFromId('Feed/Xfeed:Microstrip', 1, 2);
    pick.PickFaceFromId('UnitCell:SlotPlane', 6);
    discretefaceport.Reset();
    discretefaceport.PortNumber(1);
    discretefaceport.Label('Xfeed');
    discretefaceport.Impedance('zslot');
    discretefaceport.SetP1(1);
    discretefaceport.SetP2(1);
    discretefaceport.CenterEdge(1);
    discretefaceport.Monitor(1);
    discretefaceport.Create();
    
    pick.PickEdgeFromId('UnitCell:GroundPlane', 6, 6);
else
    solid.Subtract('Slot:Metal', 'Feed/Xfeed:TopHole');
    solid.Subtract('BackingReflector:Metal', 'Feed/Xfeed:GroundHole');
    
    pick.PickEdgeFromId('BackingReflector:Metal', 6, 6);
end

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
if(newproject)
    solid.Subtract('UnitCell:SlotPlane', 'Feed/Yfeed:TopHole');
    solid.Subtract('UnitCell:GroundPlane', 'Feed/Yfeed:GroundHole');
    if(isnumeric(patch_area) && patch_area == 0)
        solid.Delete('Feed/Yfeed:Short');
    else
        solid.Delete('Feed/Yfeed:Patch');
    end
    solid.Delete('Feed/Yfeed:SlotFeed');

    % Define Y-microstrip port.
    pick.PickEdgeFromId('Feed/Yfeed:Microstrip', 1, 2);
    pick.PickFaceFromId('UnitCell:SlotPlane', 10);
    discretefaceport.Reset();
    discretefaceport.PortNumber(3);
    discretefaceport.Label('Yfeed');
    discretefaceport.Impedance('zslot');
    discretefaceport.SetP1(1);
    discretefaceport.SetP2(1);
    discretefaceport.CenterEdge(1);
    discretefaceport.Monitor(1);
    discretefaceport.Create();

    pick.PickEdgeFromId('UnitCell:GroundPlane', 8, 8);
else
    solid.Subtract('Slot:Metal', 'Feed/Yfeed:TopHole');
    solid.Subtract('BackingReflector:Metal', 'Feed/Yfeed:GroundHole');
    
    pick.PickEdgeFromId('BackingReflector:Metal', 8, 8);
end

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


if(newproject)
    CST.BuildSchematic(dsproject, 80, inf);
    CST.AddSchematicTask(dsproject);
end

project.Rebuild();