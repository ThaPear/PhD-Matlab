function [hFig, hAx] = a_figS(number)
    if(nargin < 1)
        number = [];
    end
    [hFig, hAx] = a_fig(number);
    
    xlabel(hAx, 'Frequency [GHz]');
    ylabel(hAx, '|\Gamma| [dB]');
    
    ylim(hAx, [-30 0]);
    hAx.YLimMode = 'manual';
end