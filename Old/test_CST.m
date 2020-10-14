close all;
clearvars -except project dsproject;
SetupPath;
clc;

if(~exist('project', 'var'))
    project = CST.InitializeUnitCellProject();
    dsproject = CST.Application.ActiveDS();
    CST.BuildSchematic(dsproject, 80, 1e-12);
    CST.AddSchematicTask(dsproject);
end

% project.StoreGlobalDataValue('matlabfcn', 'ReportInformation "TEST"');
% project.RunScript(GetFullPath('CST Interface\Bas\test.bas'))

% block = dsproject.Block();
% block.Reset();
% block.Name('CAP1');
% block.StartPropertyIteration();
% propertyname = ' ';
% while(~isempty(propertyname))
%     [propertyname, type, value] = block.GetNextProperty();
%     [propertyname, ' - ', type, ' - ', value]
% end

% simulationtask = dsproject.SimulationTask();
% simulationtask.Name('SPara1');
% [amplitude, fmin, fmax] = simulationtask.GetGaussProperties('1')

% boundary = project.Boundary();
% [xmin, xmax, ymin, ymax, zmin, zmax] = boundary.GetCalculationBox();
% [valid, theta, phi, direction] = boundary.GetUnitCellScanAngle()

% fdsolver = project.FDSolver();
% [phistart, phiend, nphisteps, thetastart, thetaend, nthetasteps, einctheta, eincphi, activation] = fdsolver.GetRCSSweepProperties();

% pick = project.Pick();
% [bool, x, y, z] = pick.GetPickpointCoordinatesByIndex(0)
% [shapename, edgeid, vertexid] = pick.GetPickedEdgeByIndex(0)
% [shapename, faceid] = pick.GetPickedFaceByIndex(0)

plot1d = project.Plot1D();
project.SelectTreeItem('1D Results\S-Parameters\SZmax(1),Zmax(1)');
% project.SelectTreeItem('1D Results\S-Parameters\S1,1');
plot1d.SetCurveLimit(1, 25)
plot1d.Plot();
% [enabled, curvelimit] = plot1d.GetCurveLimit()

% plot1d.DeleteAllBackGroundShapes();
% plot1d.AddThinBackGroundLine(28, -200, 28, 0);
% plot1d.Plot();


% [bool, x, y, z, vxre, vyre, vzre, vxim, vyim, vzim] = project.GetFieldVector()

%%
% project = CST_Util.InitializeProject();
% project.StoreParameter('dx', 13.5);
% project.StoreParameter('dy', 13.5);
% project.StoreParameter('hback', 6.5);
% project.StoreParameter('slot_s0', 0.25);
% 
% % cav = CST.Cavity(5e-3, 11e-3);
% % cav.BuildCST(project);
% 
% f0 = 10e9;
% lambda0 = Constants.c0 / f0;
% 
% %% Figure 6
% %     Requires theta = 60, f0 = 5e9.
% p  = 13.5e-3/4;
% %               0-1    1-2    2-3    3-4    4-5    5-N
% ds = lambda0 * [0.006, 0.012, 0.012, 0.012, 0.012, 0.006];
% ss = p       * [       0.000, 0.000, 0.000, 0.000  0.000];
% %               1      2      3      4      5
% ws = lambda0 * [0.010, 0.015, 0.020, 0.025, 0.030];
% 
% 
% erhosts = ones(size(ds)) * 1;
% slab = ADS(p, ds, ss, ws, erhosts);
% slab.BuildCST(project);
% slab.BuildCST(project);


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
% CST.BuildCoaxFeed(project, dms, wms, lms, ...
%                 core_radius, core_tophole_radius, core_transition_radius, ...
%                 shield_radius, shield_distance, shield_startangle, shield_totalangle, shield_Nvias, ...
%                 cylinder_height, cylinder_Nvias, cylinder_angleoffset, cylinder_connector_radius)

% [f, parameters, data] = CST.LoadData('Temp/mat.mat');

%%
% project = CST.Application.Active3D();
% 
% resultname = '1D Results\Port Information\Line Impedance\1(1)';
% 
% basefilename = 'Temp/coax_impedancelookup';
% exportfileextension = '.txt';
% 
% Zs = zeros(1, 3300);
% success = project.SelectTreeItem(resultname);
% for(i = 1:25:3300)
%     ret = project.ResultNavigatorRequest('set selection', num2str(i:i+24));
%     asciiexport = project.ASCIIExport();
%     exportfilename = [basefilename, '_', num2str(i, '%04i'), '-', num2str(i+24, '%04i')];
%     exportfilename = GetFullPath(exportfilename);
%     asciiexport.Reset();
%     asciiexport.FileName([exportfilename, exportfileextension]);
%     asciiexport.Execute();
%     [f, parameters, values] = CST.LoadData([exportfilename, exportfileextension]);
%     Zs(i:i+24) = values.';
% end
% 
% 
% 
% %%
% parameters = readtable('Temp/result_navigator.csv', 'Delimiter', ';');
% feed_shield_totalangle = [cellfun(@str2double, parameters.feed_shield_totalangle)];
% feed_shield_Nvias = [cellfun(@str2double, parameters.feed_shield_Nvias)];
% feed_shield_distance = [cellfun(@str2double, parameters.feed_shield_distance)];
% 
% viadistance = 2 .* feed_shield_distance .* sind(feed_shield_totalangle ./ feed_shield_Nvias ./ 2);
% % Zs(viadistance < 0.3) = nan;
% 
% for(n = 2:6)
%     runIDs = [cellfun(@str2double, parameters.x3DRunID)];
%     runIDs = runIDs(feed_shield_Nvias == n);
% 
%     Zsmat = [];
%     xs = [];
%     ys = [];
%     for(i = runIDs.')
%         xs(feed_shield_totalangle(i)/5-3) = feed_shield_totalangle(i);
%         ys(floor(feed_shield_distance(i)*20)-4) = feed_shield_distance(i);
%         Zsmat(feed_shield_totalangle(i)/5-3, floor(feed_shield_distance(i)*20)-4) = Zs(i);
%     end
%     figureex; hold off; imagesc(ys, xs, Zsmat);
%     title(['Nvias = ', num2str(n)]);
%     colorbar;
%     ax = gca;    map = ax.Colormap;    map(1,:) = 1;    ax.Colormap = map; % Make 0 (or nan) white.
%     alignplot(gcf, 4, 2, n-1, [], 2);
% end