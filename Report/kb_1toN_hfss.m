SetupPath;
clear;
close all;

folder = 'h:\Server\Exports\';

% Plot size: 546x242

%%
close all
filename = '1-to-8-lossy-';
kb_1toN_hfss_plot(folder, filename);

filename = '1-to-16-lossy-';
kb_1toN_hfss_plot(folder, filename);

filename = '1-to-32-lossy-';
kb_1toN_hfss_plot(folder, filename);

filename = '1-to-64-lossy-';
kb_1toN_hfss_plot(folder, filename);

filename = '1-to-128-lossy-';
kb_1toN_hfss_plot(folder, filename);

hFig = figure(1); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(5:-1:1), {'1-to-8', '1-to-16', '1-to-32', '1-to-64', '1-to-128'}, 'NumColumns', 5);

hFig = figure(3); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(5:-1:1), {'1-to-8', '1-to-16', '1-to-32', '1-to-64', '1-to-128'}, 'NumColumns', 5);

hFig = figure(4); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(5:-1:1), {'1-to-8', '1-to-16', '1-to-32', '1-to-64', '1-to-128'}, 'NumColumns', 5);

%%
close all
filename = '1-to-8-lossy-back2back-';
kb_1toN_hfss_plot(folder, filename);

filename = '1-to-16-lossy-back2back-';
kb_1toN_hfss_plot(folder, filename);

hFig = figure(1); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(2:-1:1), {'1-to-8-to-1', '1-to-16-to-1'}, 'NumColumns', 5);

hFig = figure(2); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(2:-1:1), {'1-to-8-to-1', '1-to-16-to-1'}, 'NumColumns', 5);

hFig = figure(3); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(2:-1:1), {'1-to-8-to-1', '1-to-16-to-1'}, 'NumColumns', 5);
