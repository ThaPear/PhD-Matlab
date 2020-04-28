SetupPath;
clear;
close all;

folder = 'h:\Server\Exports\';

% Plot size: 546x242

%%
close all
filename = 'transition-';
kb_1toN_hfss_plot(folder, filename);

alignplot(figureex(1), 8, 4, 1, [], 2);
ylim([-37.5 0]);
