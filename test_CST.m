close all;
clear;
SetupPath;
clc;


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

[f, parameters, data] = CST.LoadData('Temp/mat.mat');
