function [hFig, hAx] = a_figZ(number, fs, zfeed)
    if(nargin < 1)
        number = [];
    end
    [hFig, hAx] = a_fig(number);
    
    xlabel(hAx, 'Frequency [GHz]');
    clrname = caller('name');
    if(clrname(1) == 'a' || clrname(1) == 'b' || clrname(1) == 'c')
        ylabel(hAx, 'Input Impedance [\Omega]');
    else
        ylabel(hAx, 'Active Input Impedance [\Omega]');
    end
    
    hAx.ColorOrder = reshape(repmat(lines(7), 1, 2).', [], 14).';
    if(nargin > 2)
        plot(hAx, [min(fs), max(fs)]./1e9, [zfeed, zfeed], 'k--');
    else
        warning('No feed impedance provided.\n');
    end
end