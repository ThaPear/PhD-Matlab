SetupPath;
clear;
close all;

folder = 'h:\Server\Exports\';

% Plot size: 546x242

%%
close all
filename = '1-to-8-lossy-';
kb_1toN_hfss_plot(folder, filename);

[hAx, hFig] = figureex(1);
addlegendentry(hAx, '1-to-8');
drawnow; movelegend(hAx, 'n');

filename = '1-to-16-lossy-';
kb_1toN_hfss_plot(folder, filename);
addlegendentry(hAx, '1-to-16');
drawnow; movelegend(hAx, 'n');


filename = '1-to-32-lossy-';
kb_1toN_hfss_plot(folder, filename);
addlegendentry(hAx, '1-to-32');
drawnow; movelegend(hAx, 'n');


filename = '1-to-64-lossy-';
kb_1toN_hfss_plot(folder, filename);
addlegendentry(hAx, '1-to-64');
drawnow; movelegend(hAx, 'n');


filename = '1-to-128-lossy-';
kb_1toN_hfss_plot(folder, filename);
addlegendentry(hAx, '1-to-128');
drawnow; movelegend(hAx, 'n');


% hFig = figure(1); hAx = hFig.CurrentAxes;
% legend(hAx, hAx.Children(5:-1:1), {'1-to-8', '1-to-16', '1-to-32', '1-to-64', '1-to-128'}, 'NumColumns', 5);

hFig = figure(3); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(5:-1:1), {'1-to-8', '1-to-16', '1-to-32', '1-to-64', '1-to-128'}, 'NumColumns', 2);
legendlinelength(hAx, 15);
ylabel('Losses [dB]');

hFig = figure(4); hAx = hFig.CurrentAxes;
legend(hAx, hAx.Children(5:-1:1), {'1-to-8', '1-to-16', '1-to-32', '1-to-64', '1-to-128'}, 'NumColumns', 2);
legendlinelength(hAx, 15);

l = [14 20 30 38 55];
L = [-0.35 -0.54 -0.72 -0.93 -1.32];
interp1(l, L, 140, 'linear', 'extrap')
figureex;
    plot(l, L, 'o-');
    % 1-to-1024 length is 140mm
