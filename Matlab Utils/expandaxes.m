% ---- expandaxes.m

function expandaxes(hFig)
    if(nargin < 1)
        hFig = gcf;
    end
                          % Find axes     % Filter suptitle
    hAxes = findobj(hFig, 'Type', 'Axes', '-not', 'Tag', 'suptitle');
    suptitle = findobj(hFig, 'Type', 'Axes', '-and', 'Tag', 'suptitle');
    hasSuptitle = numel(suptitle);

    Naxes = numel(hAxes);
    
    Inset = hAxes(1).TightInset;
    if(Naxes == 1)
        hAxes.Position = [Inset(1:2), 1-Inset(1)-Inset(3), 1-Inset(2)-Inset(4)];
    elseif(Naxes == 2)
        i1 = 1;
        i2 = 2;
        % Ensure i1 points to the left subfigure.
        if(hAxes(i1).Position(1) > hAxes(i2).Position(2))
            [i1, i2] = deal(i2, i1);
        end
        hAxes(i1).Position = [Inset(1:2), 0.5-Inset(1)-Inset(3), 1-Inset(2)-Inset(4)];
        hAxes(i2).Position = [[0.5,0]+Inset(1:2), 0.5-Inset(1)-Inset(3), 1-Inset(2)-Inset(4)];
        
        %% Add spacing between the top and the axes if a suptitle is present.
        if(hasSuptitle)
            suptitle.Position(2) = 1 - suptitle.TightInset(4);
            suptitle.Position(4) = suptitle.TightInset(4);
            hAxes(i1).Position(4) = hAxes(i1).Position(4) - 2*suptitle.Position(4);
            hAxes(i2).Position(4) = hAxes(i2).Position(4) - 2*suptitle.Position(4);
        end
    else
        error('expandaxes(): Only supports 1 or 2 subplots.');
    end
%         for n = 1:Naxes
%             pos1(n) = h.Children(n).Position(1);
%             pos2(n) = h.Children(n).Position(2);
%         end
%         Ncols = numel(unique(pos1));
%         Nrows= numel(unique(pos2));
end