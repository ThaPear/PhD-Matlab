clear;
close all;
SetupPath;


linewidth = 1;
axlinewidth = 1;
nplotx = 6;
nploty = 4;
xlims = [12 32];

roughness = 0.0;

% Connector
folder = 'h:\Server\Exports\';
[parameters, SCST] = CST.LoadData([folder, sprintf('Connector_TR_Deembed_%gum.s2p', roughness)]);
fCST = parameters.frequencies.';

roughness = 1.3;

clrs = lines(7);
clrs(4,:) = [];

file = '16.s2p';
[parameters, data] = CST.LoadData(file);
fs = parameters.frequencies/1e9;

s11 = squeeze(data(1,1,:));
s12 = squeeze(data(1,2,:));
s21 = squeeze(data(2,1,:));
s22 = squeeze(data(2,2,:));
[hFig, hAx] = figureex;
    if(length(hAx.Children) < 2)
        patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    hAx.ColorOrder = clrs;
    alignplot(hFig, nplotx, nploty, hFig.Number, [], 1);
%         hAx.Title.String = strrep(file, '.s2p', '');
    plot(hAx, fs, 20*log10(abs(s11)), 'LineWidth', linewidth);
    addlegendentry(hAx, 'S11');
    plot(hAx, fs, 20*log10(abs(s12)), 'LineWidth', linewidth);
    addlegendentry(hAx, 'S12');
%     plot(hAx, fs, 20*log10(abs(s21)), 'LineWidth', linewidth);
%     addlegendentry(hAx, 'S21');
    plot(hAx, fs, 20*log10(abs(s22)), 'LineWidth', linewidth);
    addlegendentry(hAx, 'S22');

    xlabel(hAx, 'Frequency [GHz]');
    ylabel(hAx, 'S_{ij} [dB]');
    xlim(hAx, xlims);
    ylim(hAx, [-30 0]);
    movelegend(hAx, 'se');
    hAx.LineWidth = axlinewidth;
    
folder = 'h:\Server\Exports\';
filename = sprintf('1-to-16-lossy-back2back-rough%g-', roughness);
[f, ~, S11HFSS] = HFSS.LoadData([folder, filename, 'S11.csv']);
[f, ~, S12HFSS] = HFSS.LoadData([folder, filename, 'S12.csv']);
SHFSS = struct('s11', 10.^(S11HFSS./20), 's12', 10.^(S12HFSS./20), 's21', 10.^(S12HFSS./20), 's22', 10.^(S11HFSS./20));
ABCDHFSS = S2ABCD(SHFSS, 80, 80);
S11i = interp1(fCST./1e9, squeeze(SCST(1,1,:)), f);
S12i = interp1(fCST./1e9, squeeze(SCST(1,2,:)), f);
S21i = interp1(fCST./1e9, squeeze(SCST(2,1,:)), f);
S22i = interp1(fCST./1e9, squeeze(SCST(2,2,:)), f);

SCST = struct('s11', S11i, 's12', S12i, 's21', S21i, 's22', S22i);
Z1 = 80;%50.974587;
Z2 = 80;%77.987955;
ABCDCST = S2ABCD(SCST, Z1, Z2);

ABCDtotal = ABCDCST.mul(ABCDHFSS.mul(ABCDCST.reverse()));

Stotal = ABCD2S(ABCDtotal, 50, 50);
% Stotal = SHFSS;
    
[hFig2, hAx2] = figureex;
    alignplot(hFig2, nplotx, nploty, hFig2.Number, [], 1);
%     title(hAx2, '16');
    plot(hAx2, fs, 10*log10(abs(s11).^2+abs(s12).^2), '-', 'Color', clrs(1,:), 'LineWidth', linewidth);
    addlegendentry(hAx2, 'Measurement');

%     plot(hAx, f, S11HFSS, 'o', 'Color', clrs(1,:), 'LineWidth', linewidth);
    plot(hAx, f, 20*log10(abs(Stotal.s11)), 'k:', 'LineWidth', linewidth*2);
    addlegendentry(hAx, 'Simulation');
    plot(hAx, f, 20*log10(abs(Stotal.s12)), 'k:', 'LineWidth', linewidth*2);
%     plot(hAx, f, S12HFSS, 'x', 'Color', clrs(2,:), 'LineWidth', linewidth);
%     plot(hAx2, f, 10*log10(10.^(S12HFSS./20).^2 + 10.^(S11HFSS./20).^2), 'x', 'Color', clrs(2,:), 'LineWidth', linewidth);
%     addlegendentry(hAx2, '0.25\mum');
%     plot(hAx2, f, 10*log10(10.^(S12HFSS./20).^2 + 10.^(S11HFSS./20).^2) + 2*10*log10(abs(S11i).^2+abs(S12i).^2), 'o', 'Color', clrs(2,:), 'LineWidth', linewidth);
    plot(hAx2, f, 10*log10(abs(Stotal.s11).^2 + abs(Stotal.s12).^2), '-', 'Color', clrs(2,:), 'LineWidth', linewidth);
    addlegendentry(hAx2, 'Simulation');% + connector');
        
    xlim(hAx2, xlims);
    ylim(hAx2, [-3 0]);
    xlabel(hAx2, 'Frequency [GHz]');
    ylabel(hAx2, 'Efficiency');
    
    legend(hAx, 'NumColumns', 5)
    legendlinelength(hAx, 15);
    drawnow;
    movelegend(hAx, 's');
        
        
    dispex('%s: %.2f at 12, %.2f at 32.\n', file, 20*log10(s12(fs == 12)), 20*log10(s12(fs == 32)));
    
%     figureex; plot(fCST./1e9, 20*log10(abs(SCST.s11)));
%     plot(fCST./1e9, 20*log10(abs(SCST.s12)));


for(ifig = 1:100)
    if(ishandle(ifig))
        hFig = figure(ifig);
        hAx = hFig.CurrentAxes;
        cropfigure(hFig);
        legendlinelength(hFig, 15);
        drawnow;
        legend(hAx, 'NumColumns', 5)
        movelegend(hFig, 'se');
    end
end

