function [hFig, hAx] = a_figVSWR(number)
    if(nargin < 1)
        number = [];
    end
    [hFig, hAx] = a_fig(number);
    
    xlabel(hAx, 'Frequency [GHz]');
    clrname = caller('name');
    if(clrname(1) == 'a' || clrname(1) == 'b' || clrname(1) == 'c')
        ylabel(hAx, 'VSWR');
    else
        ylabel(hAx, 'Active VSWR');
    end
    
    ylim(hAx, [1 8]);
    hAx.YLimMode = 'manual';
end