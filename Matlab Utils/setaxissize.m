function setaxissize(hAxes, szx, szy)
    % If it's a figure handle
    if(isa(hAxes, 'matlab.ui.Figure'))
        hFig = hAxes;
        hAxes = hFig.CurrentAxes;
    elseif(isa(hAxes, 'matlab.graphics.axis.Axes'))
        hFig = hAxes.Parent;
    end
    
    hAxes.Units = 'pixels';
    hFig.Units = 'pixels';
    
    figSz = hFig.Position(3:4);
    axSz = hAxes.Position(3:4);
    
    padding = figSz - axSz;
    hAxes.Position(3:4) = [szx, szy];
    hFig.Position(3:4) = [szx, szy] + padding;
end