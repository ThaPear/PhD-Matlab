function [hFig, hAx] = a_fig(number)
    linewidth = 1;
    axlinewidth = 1;
    
    if(nargin > 0 && ~isempty(number))
        [hFig, hAx] = figureex(number);
    else
        [hFig, hAx] = figureex;
    end
    
    if(length(hAx.Children) < 2)
        patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    
%     alignplot(hFig, 4, 3, hFig.Number, [], 1);
    alignplot(hFig, 8, 4, hFig.Number, [], 2);
    hAx.LineWidth = axlinewidth;
    hAx.FontSize = 9;
    set(hAx, 'DefaultLineLineWidth', linewidth);
    
    xlim(hAx, [12 32]);
    hAx.XLimMode = 'manual';
end