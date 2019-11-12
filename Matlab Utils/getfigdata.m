function [X, Y] = getfigdata(fig)
% ---- GetFigData.m

X = [];
Y = [];
ax = fig.CurrentAxes;
for(i = 1:length(ax.Children))
    if(isa(ax.Children(i), 'matlab.graphics.chart.primitive.Line'))
        X = [X; ax.Children(i).XData];
        Y = [Y; ax.Children(i).YData];
    end
end

% close all;
% clear variables;
% clc;

% uiopen('E:\Zooi\Stack\Schoolwerk\TU\Master Year 1\EE4580 - Quasi Optical Systems\Lab Sessions\Lab 2\figures\feed_farfield.fig',1)
% h = gcf;
% axesObjs = h.CurrentAxes;
% dataObjs = axesObjs.Children;
% xdata1 = get(dataObjs(1), 'XData');
% ydata1 = get(dataObjs(1), 'YData');
% xdata2 = get(dataObjs(2), 'XData');
% ydata2 = get(dataObjs(2), 'YData');