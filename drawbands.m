function drawbands(hAx)
    if(~isa(hAx, 'matlab.graphics.axis.Axes'))
        if(isa(hAx, 'matlab.ui.Figure'))
            hAx = hAx.CurrentAxes;
        else
            error('drawbands: Expected axes or figure handle, received %s.\n', class(hAx));
        end
    end

    ku = 0;
    ka = 0;
    for(ichild = 1:length(hAx.Children))
        child = hAx.Children(ichild);
        if(isa(child, 'matlab.graphics.primitive.Patch'))
            % Check if the patch coordinates match those of either band.
            if(length(child.XData) == 4 && all(child.XData == [13.75 14.5 14.5 13.75].'))
                ku = 1;
            end
            if(length(child.XData) == 4 && all(child.XData == [28 31 31 28].'))
                ka = 1;
            end
        end
    end
    if(~ku)
        patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    if(~ka)
        patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
end