SetupPath;
clear;
close all;

folder = 'h:\Server\Exports\';

% Plot size: 546x242

%%
close all
filename = '1-to-8-lossy-back2back-';
kb_1toN_hfss_plot(folder, filename);

filename = '1-to-16-lossy-back2back-';
kb_1toN_hfss_plot(folder, filename);

filename = '1-to-32-lossy-back2back-';
kb_1toN_hfss_plot(folder, filename);

hFig = figure(1); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(3:-1:1), {'1-to-8-to-1', '1-to-16-to-1', '1-to-32-to-1'}, 'NumColumns', 1);
delete(hAx.Children(1));
movelegend(hAx, 'n');

hFig = figure(2); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(3:-1:1), {'1-to-8-to-1', '1-to-16-to-1', '1-to-32-to-1'}, 'NumColumns', 1);
delete(hAx.Children(1));
movelegend(hAx, 's');

hFig = figure(3); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(3:-1:1), {'1-to-8-to-1', '1-to-16-to-1', '1-to-32-to-1'}, 'NumColumns', 1);
delete(hAx.Children(1));
movelegend(hAx, 's');
