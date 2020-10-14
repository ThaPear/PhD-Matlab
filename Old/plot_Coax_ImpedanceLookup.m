close all;
clear;
SetupPath;
clc;

% parameters = readtable('Temp/result_navigator.csv', 'Delimiter', ';');
% feed_shield_totalangle = [cellfun(@str2double, parameters.feed_shield_totalangle)];
% feed_shield_Nvias = [cellfun(@str2double, parameters.feed_shield_Nvias)];
% feed_shield_distance = [cellfun(@str2double, parameters.feed_shield_distance)];

% Zs(viadistance < 0.3) = nan;

load('Coax_ImpedanceLookup.mat');


for(n = 2:6)
    runIDs = [cellfun(@str2double, parameters.x3DRunID)];
    runIDs = runIDs(feed_shield_Nvias == n);

    Zsmat = [];
    xs = [];
    ys = [];
    for(i = runIDs.')
        xs(feed_shield_totalangle(i)/5-3) = feed_shield_totalangle(i);
        ys(floor(feed_shield_distance(i)*20)-4) = feed_shield_distance(i);
        Zsmat(feed_shield_totalangle(i)/5-3, floor(feed_shield_distance(i)*20)-4) = Zs(i);
    end
    figureex; hold off; imagesc(ys, xs, Zsmat);
    title(['Nvias = ', num2str(n)]);
    colorbar;
    ax = gca;    map = ax.Colormap;    map(1,:) = 1;    ax.Colormap = map; % Make 0 (or nan) white.
    alignplot(gcf, 4, 2, n-1, [], 2);
end