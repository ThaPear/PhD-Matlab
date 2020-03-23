% ---- movelegend.m

function movelegend(hAx, pos)
%     hFig = gcf;
    switch(nargin)
        case 0
            pos = 'sw';
            hAx = gca;
        case 1
            pos = hAx;
            hAx = gca;
    end
    if(~isa(hAx, 'matlab.graphics.axis.Axes'))
        if(isa(hAx, 'matlab.ui.Figure'))
            hAx = hAx.CurrentAxes;
        else
            error('Invalid hAx type %s.\n', class(hAx));
        end
    end
    
    hLeg = legend(hAx);
    hLeg.Units = hAx.Units;
    switch(pos)
        case 'n'
            hLeg.Position(1) = hAx.Position(1) + hAx.Position(3) / 2 - hLeg.Position(3) / 2;
            hLeg.Position(2) = hAx.Position(2) + hAx.Position(4) - hLeg.Position(4);
        case 'ne'
            hLeg.Position(1) = hAx.Position(1) + hAx.Position(3) - hLeg.Position(3);
            hLeg.Position(2) = hAx.Position(2) + hAx.Position(4) - hLeg.Position(4);
        case 'e'
            hLeg.Position(1) = hAx.Position(1) + hAx.Position(3) - hLeg.Position(3);
            hLeg.Position(2) = hAx.Position(2) + hAx.Position(4) / 2 - hLeg.Position(4) / 2;
        case 'se'
            hLeg.Position(1) = hAx.Position(1) + hAx.Position(3) - hLeg.Position(3);
            hLeg.Position(2) = hAx.Position(2);
        case 's'
            hLeg.Position(1) = hAx.Position(1) + hAx.Position(3) / 2 - hLeg.Position(3) / 2;
            hLeg.Position(2) = hAx.Position(2);
        case 'sw'
            hLeg.Position(1) = hAx.Position(1);
            hLeg.Position(2) = hAx.Position(2);
        case 'w'
            hLeg.Position(1) = hAx.Position(1);
            hLeg.Position(2) = hAx.Position(2) + hAx.Position(4) / 2 - hLeg.Position(4) / 2;
        case 'nw'
            hLeg.Position(1) = hAx.Position(1);
            hLeg.Position(2) = hAx.Position(2) + hAx.Position(4) - hLeg.Position(4);
        otherwise
            disp('movelegend: Only ''n'', ''ne'', ''e'', ''se'', ''s'', ''sw'', ''w'', ''nw'' are supported');
    end
end