% dms = 0.254e-3;
% wms = 0.23e-3;
% lms = 1.5e-3;
% core_radius = 0.1e-3;
% core_tophole_radius = 0.1e-3;
% core_transition_radius = 0.05e-3;
% shield_radius = 0.1e-3;
% shield_distance = 0.5e-3;
% shield_startangle = ['90 - shield_totalangle/2'];
% shield_totalangle = 100;
% shield_Nvias = 3;
% cylinder_height = 0.916e-3;
% cylinder_Nvias = [];
% cylinder_angleoffset = 30;
% cylinder_connector_radius = shield_radius + 0.05e-3;
% project = CST.Application.Active3D();
% CST.BuildCoaxFeed_ViaCylinder(project, dms, wms, lms, ...
%                 core_radius, core_tophole_radius, core_transition_radius, ...
%                 shield_radius, shield_distance, shield_startangle, shield_totalangle, shield_Nvias, ...
%                 cylinder_height, cylinder_Nvias, cylinder_angleoffset, cylinder_connector_radius)


function BuildCoaxFeed(project, dms, wms, lms, wms_buried, lms_buried, ...
                core_radius, core_tophole_radius, core_transition_radius, ...
                shield_radius, shield_distance, shield_startangle, shield_totalangle, shield_Nvias, ...
                patch_area, patch_angle, patch_l1)
    % It is assumed that the WCS is correctly positioned in the center of the feeding point on the
    % slot plane.
    % TODO: Implement cylinder_Nvias
    
    if(nargin == 0)
        warning('No arguments provided.');
        return;
    end
    
    %% Create relevant parameters in the project.
    % Parameters that should already exist.
    wfeedname = 'slot_feedwidth';
    lfeedname = 'slot_feedlength';
    hbackname = 'hback';
    % slot_s0
    % dx, dy
    
    % New parameters
    dmsname = 'feed_ms_depth';                      project.StoreParameter(dmsname, dms);
    wmsname = 'feed_ms_width';                      project.StoreParameter(wmsname, wms);
    lmsname = 'feed_ms_length';                     project.StoreParameter(lmsname, lms);
    thms2name = 'feed_ms_buried_angle';             project.StoreParameter(thms2name, 0);
    wms2name = 'feed_ms_buried_width';              project.StoreParameter(wms2name, wms_buried);
    lms2name = 'feed_ms_buried_length';             project.StoreParameter(lms2name, lms_buried);
    rcorename = 'feed_core_radius';                 project.StoreParameter(rcorename, core_radius);
    rtopholename = 'feed_core_tophole_radius';      project.StoreParameter(rtopholename, core_tophole_radius);
    rcore2msname = 'feed_core_transition_radius';   project.StoreParameter(rcore2msname, core_transition_radius);
    hcorename = 'feed_core_height';                 project.StoreParameter(hcorename, 'hback');
    rshieldname = 'feed_shield_radius';             project.StoreParameter(rshieldname, shield_radius);
    dshieldname = 'feed_shield_distance';           project.StoreParameter(dshieldname, shield_distance);
    thstartname = 'feed_shield_startangle';         project.StoreParameter(thstartname, shield_startangle);
    thtotalname = 'feed_shield_totalangle';         project.StoreParameter(thtotalname, shield_totalangle);
    nvianame = 'feed_shield_Nvias';                 project.StoreParameter(nvianame, shield_Nvias);
    groundholeradiusname = 'feed_groundhole_radius';project.StoreParameter(groundholeradiusname, [dshieldname, '-', rshieldname, '-0.1']);
    
    apatchname = 'feed_patch_area';                 project.StoreParameter(apatchname, patch_area);
    wmsfeedname = 'feed_ms_feedgap_width';          project.StoreParameter(wmsfeedname, ['0.8*', lfeedname]);
    thpatchname = 'feed_patch_angle';               project.StoreParameter(thpatchname, patch_angle);
    l1patchname = 'feed_patch_l1';                  project.StoreParameter(l1patchname, patch_l1);
    l2patchname = 'feed_patch_l2';                  project.StoreParameter(l2patchname, ...
        ['(', apatchname, '-', wmsfeedname, '*', l1patchname, '*cosd(', thpatchname, ')-(', l1patchname, ')^2*cosd(', thpatchname, ')*sind(', thpatchname, '))/(2*', l1patchname, '*sind(', thpatchname, ')+', wmsfeedname, ')']);
%         (patch_area-feed_ms_w*patch_l*cosd(patch_angle)-(patch_l)^2*cosd(patch_angle)*sind(patch_angle))/(2*patch_l*sind(patch_angle)+feed_ms_w)

    %% CST interface objects
%     arc = project.Arc();
    brick = project.Brick();
    component = project.Component();
    cylinder = project.Cylinder();
    tracefromcurve = project.TraceFromCurve();
    transform = project.Transform();
    wcs = project.WCS();
    polygon3d = project.Polygon3D();
    extrudecurve = project.ExtrudeCurve();
    
    %% Create microstrip - old
    % Microstrip
%     brick.Reset();
%     brick.Component('Feed/Xfeed');
%     brick.Name('Microstrip');
%     brick.Material('PEC');
%     brick.Xrange(['-', wmsname, '/2'], [wmsname, '/2']);
%     brick.Yrange(['-', wfeedname, '/2 - ', lmsname], ['-', wfeedname, '/2']);
%     brick.Zrange(['-', dmsname], ['-', dmsname]);
%     brick.Create();
    % Microstrip buried
%     brick.Reset();
%     brick.Component('Feed/Xfeed');
%     brick.Name('Microstrip_Buried');
%     brick.Material('PEC');
%     brick.Xrange(['-', wms2name, '/2'], [wms2name, '/2']);
%     brick.Yrange(['-', wfeedname, '/2 - (', lmsname, ' + ', lms2name, ')'], ['-', wfeedname, '/2 - ', lmsname]);
%     brick.Zrange(['-', dmsname], ['-', dmsname]);
%     brick.Create();
    
    %% Create microstrip
    % Microstrip
    m1_points = {{  '0', ...
                    ['-', wfeedname, '/2'], ...
                    ['-', dmsname]}, ...
                 {  '0', ...
                    ['-', wfeedname, '/2 - ', lmsname], ...
                    ['-', dmsname]}, ...
                };
    polygon3d.Reset();
    polygon3d.Curve('Feed');
    polygon3d.Name('Microstrip');
    for(i = 1:length(m1_points))
        polygon3d.Point(m1_points{i}{1}, m1_points{i}{2}, m1_points{i}{3});
    end
    polygon3d.Create();
    
    tracefromcurve.Reset();
    tracefromcurve.Component('Feed/Xfeed');
    tracefromcurve.Name('Microstrip');
    tracefromcurve.Material('PEC');
    tracefromcurve.Curve('Feed:Microstrip');
    tracefromcurve.Thickness(0)
    tracefromcurve.Width(wmsname);
    tracefromcurve.DeleteCurve(1);
    tracefromcurve.RoundEnd(['" & (', thms2name, ' <> 0) & "']);
    tracefromcurve.Create();
    
    % Microstrip buried
    m2_points = {{  '0', ...
                    ['-', wfeedname, '/2 - ', lmsname], ...
                    ['-', dmsname]}, ...
                 {  ['sin(', thms2name, '*pi/180)*', lms2name], ...
                    ['-', wfeedname, '/2 - ', lmsname, ' - cos(', thms2name, '*pi/180)*', lms2name], ...
                    ['-', dmsname]}, ...
                };
    polygon3d.Reset();
    polygon3d.Curve('Feed');
    polygon3d.Name('Microstrip');
    for(i = 1:length(m1_points))
        polygon3d.Point(m2_points{i}{1}, m2_points{i}{2}, m2_points{i}{3});
    end
    polygon3d.Create();
    
    tracefromcurve.Reset();
    tracefromcurve.Component('Feed/Xfeed');
    tracefromcurve.Name('Microstrip_Buried');
    tracefromcurve.Material('PEC');
    tracefromcurve.Curve('Feed:Microstrip');
    tracefromcurve.Thickness(0)
    tracefromcurve.Width(wms2name);
    tracefromcurve.DeleteCurve(1);
    tracefromcurve.RoundStart(['" & (', thms2name, ' <> 0) & "']);
    tracefromcurve.Create();
    
    % Create feed part
    brick.Reset();
    brick.Name('SlotFeed');
    brick.Component('Feed/Xfeed');
    brick.Material('PEC');
    brick.Xrange(['-', wmsfeedname, '/2'], [wmsfeedname, '/2']);
    brick.Yrange(['-', wfeedname, '/2'], [wfeedname, '/2']);
    brick.Zrange(['-', dmsname], ['-', dmsname]);
    brick.Create();
    
    if(patch_area == 0)
        %% Create short for the feed
        brick.Reset();
        brick.Name('Short');
        brick.Component('Feed/Xfeed');
        brick.Material('PEC');
        brick.Xrange(['-', wmsfeedname, '/2'], [wmsfeedname, '/2']);
        brick.Yrange([wfeedname, '/2'], [wfeedname, '/2']);
        brick.Zrange(['-', dmsname], '0');
        brick.Create();
    else
        %% Create patch for the feed
        x1 = [wmsfeedname, '/2'];
        x2 = [x1, ' + ', l1patchname, '*sind(', thpatchname, ')'];
        y1 = [wfeedname, '/2'];
        y2 = [y1, ' + ', l1patchname, '*cosd(', thpatchname, ')'];
        y3 = [y2, ' + ', l2patchname];
        z1 = ['-', dmsname];
        % Draw outline
        polygon3d.Reset();
%         polygon3d.Version(10);
        polygon3d.Name('Patch');
        polygon3d.Curve('Feed');
        polygon3d.Point(x1, y1, z1);
        polygon3d.Point(x2, y2, z1);
        polygon3d.Point(x2, y3, z1);
        polygon3d.Point(['-(', x2, ')'], y3, z1);
        polygon3d.Point(['-(', x2, ')'], y2, z1);
        polygon3d.Point(['-(', x1, ')'], y1, z1);
        polygon3d.Point(x1, y1, z1);
        polygon3d.Create();
        
        % Extrude to fill
        extrudecurve.Reset();
        extrudecurve.Name('Patch');
        extrudecurve.Component('Feed/Xfeed');
        extrudecurve.Material('PEC');
        extrudecurve.Thickness('0.0');
%         extrudecurve.Twistangle('0.0');
%         extrudecurve.Taperangle('0.0');
        extrudecurve.DeleteProfile('True');
        extrudecurve.Curve('Feed:Patch');
        extrudecurve.Create();
    end
    
    wcs.MoveWCS('local', m2_points{end}{1}, m2_points{end}{2}, 0);
    
    %% Create coax core
    cylinder.Reset();
    cylinder.Component('Feed/Xfeed');
    cylinder.Name('Core');
    cylinder.Material('PEC');
    cylinder.Axis('z');
    cylinder.Xcenter(0);
    cylinder.Ycenter(0);%['-', wfeedname, '/2 - (', lmsname, ' + ', lms2name, ')']);
    cylinder.Zrange(['-', hbackname], ['-', hbackname, ' + ', hcorename]);
    cylinder.Innerradius(0);
    cylinder.Outerradius(rcorename);
    cylinder.Create();
    
    %% Create microstrip-to-core transition
    cylinder.Reset();
    cylinder.Component('Feed/Xfeed');
    cylinder.Name('Core2ms');
    cylinder.Material('PEC');
    cylinder.Axis('z');
    cylinder.Xcenter(0);
    cylinder.Ycenter(0);%['-', wfeedname, '/2 - (', lmsname, ' + ', lms2name, ')']);
    cylinder.Zrange(['-', dmsname], ['-', dmsname]);
    cylinder.Innerradius(0);
    cylinder.Outerradius([rcorename, ' + ', rcore2msname]);
    cylinder.Create();
    
    %% Create coax shield
    % Create first shield via
    cylinder.Reset();
    cylinder.Component('Feed/Xfeed');
    cylinder.Name('Shield');
    cylinder.Material('PEC');
    cylinder.Axis('z');
    cylinder.Xcenter(dshieldname);
    cylinder.Ycenter(0);%['-', wfeedname, '/2 - (', lmsname, ' + ', lms2name, ')']);
    cylinder.Zrange(['-', hbackname], 0);
    cylinder.Innerradius(0);
    cylinder.Outerradius(rshieldname);
    cylinder.Create();
    
    % Rotate shield via to shield start position
    transform.Reset();
    transform.Name('Feed/Xfeed:Shield');
    transform.Angle(0, 0, ['-', thstartname]);
    transform.Origin('Free');
    transform.Center(0, 0, 0);%['-', wfeedname, '/2 - (', lmsname, ' + ', lms2name, ')'], 0);
    transform.Repetitions(1);
    transform.MultipleObjects(0);
    transform.GroupObjects(0);
    transform.Transform('Shape', 'Rotate');
    
    % Copy shield via to form shield
    transform.Reset();
    transform.Name('Feed/Xfeed:Shield');
    transform.Angle(0, 0, ['-', thtotalname, ' / (', nvianame, ' - 1)']);
    transform.Origin('Free');
    transform.Center(0, 0, 0);%['-', wfeedname, '/2 - (', lmsname, ' + ', lms2name, ')'], 0);
    transform.Repetitions([nvianame, ' - 1']);
    transform.MultipleObjects(1);
    transform.GroupObjects(0);
    transform.Transform('Shape', 'Rotate');

    %% Create ground hole
    cylinder.Reset();
    cylinder.Name('GroundHole');
    cylinder.Component('Feed/Xfeed');
    cylinder.Material('PEC');
    cylinder.Outerradius([groundholeradiusname]);
    cylinder.Innerradius('0.0');
    cylinder.Axis('z');
    cylinder.Zrange(['-', hbackname], ['-', hbackname]);
    cylinder.Xcenter('0');
    cylinder.Ycenter(0);%['-', wfeedname, '/2-(', lmsname, ' + ', lms2name, ')']);
    cylinder.Segments('0');
    cylinder.Create();

    %% Create top hole
    cylinder.Reset();
    cylinder.Name('TopHole');
    cylinder.Component('Feed/Xfeed');
    cylinder.Material('PEC');
    cylinder.Outerradius([rcorename, ' + ', rtopholename]);
    cylinder.Innerradius(rcorename);
    cylinder.Axis('z');
    cylinder.Zrange(0, 0);
    cylinder.Xcenter('0');
    cylinder.Ycenter(0);%['-', wfeedname, '/2-(', lmsname, ' + ', lms2name, ')']);
    cylinder.Segments('0');
    cylinder.Create();


    %% Rotate X-feed by 90 to form the Y-feed
%     wcs.Disable();
%     component.New('Feed/Yfeed');
%     transform.Reset();
%     transform.Name('Feed/Xfeed');
%     transform.Angle(0, 0, -90);
%     transform.Origin('Free');
%     transform.Center('(slot_s0-0.5)*dx', '(slot_s0-0.5)*dy', 0);
%     transform.Repetitions(1);
%     transform.MultipleObjects(1);
%     transform.GroupObjects(0);
%     transform.Destination('Feed/Yfeed');
%     transform.Material('');
%     transform.Transform('Shape', 'Rotate');
    
    %% Mirror X-feed in 0,0 to form the Y-feed
    wcs.Disable();
    component.New('Feed/Yfeed');
    transform.Reset();
    transform.Name('Feed/Xfeed');
    transform.Origin('Free');
    transform.Center(0, 0, 0);
    transform.PlaneNormal(1, -1, 0);
    transform.MultipleObjects(1);
    transform.GroupObjects(0);
    transform.Repetitions(1);
    transform.MultipleSelection(0);
    transform.Destination('Feed/Yfeed');
    transform.Material('');
    transform.Transform('Shape', 'Mirror');
end