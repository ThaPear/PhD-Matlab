clear;
close all;
SetupPath;


linewidth = 1;
axlinewidth = 1;
nplotx = 6;
nploty = 4;
xlims = [12 32];

% Connector
folder = 'h:\Server\Exports\';
[parameters, SCST] = CST.LoadData([folder, 'Connector_TR_Deembed_0um.s2p']);
fCST = parameters.frequencies.';

files = {'Line.s2p', '16.s2p', '64.s2p', '256.s2p', '1024.s2p'};
clrs = lines(7);
clrs(4,:) = [];
for(ifile = 1:length(files))
    file = files{ifile};
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
%         plot(hAx, fs, 20*log10(abs(s21)), 'LineWidth', linewidth);
%         addlegendentry(hAx, 'S21');
        plot(hAx, fs, 20*log10(abs(s22)), 'LineWidth', linewidth);
        addlegendentry(hAx, 'S22');
        
    %% 1-to-1
    if(ifile == find(strcmp(files, 'Line.s2p')))
        [hFig2, hAx2] = figureex;
        alignplot(hFig2, nplotx, nploty, hFig2.Number, [], 1);
        
        folder = 'h:\Server\Exports\';
        filename = '1-to-1-lossy_rough0um-';
        [f, ~, S11] = HFSS.LoadData([folder, filename, 'S11.csv']);
%             plot(hAx, f, S11, 'x', 'Color', clrs(1,:), 'LineWidth', linewidth);
%             addlegendentry(hAx, '0.0\mum');
        [f, ~, S12] = HFSS.LoadData([folder, filename, 'S12.csv']);
% %             plot(hAx, f, S12, 'x', 'Color', clrs(2,:), 'LineWidth', linewidth);
%             plot(hAx2, f, 10*log10(10.^(S12./20).^2 + 10.^(S11./20).^2), 'x', 'Color', clrs(1,:), 'LineWidth', linewidth);
%             addlegendentry(hAx2, '0.0\mum');
        S11i = interp1(fCST./1e9, squeeze(SCST(1,1,:)), f);
        S12i = interp1(fCST./1e9, squeeze(SCST(1,2,:)), f);
            plot(hAx2, f, 10*log10(10.^(S12./20).^2 + 10.^(S11./20).^2) + 2*10*log10(abs(S11i).^2+abs(S12i).^2), 'o', 'Color', clrs(1,:), 'LineWidth', linewidth);
            addlegendentry(hAx2, 'Simulation');% + connector');
        
        
        
        
        plot(hAx2, fs, 10*log10(abs(s11).^2+abs(s12).^2), '-', 'Color', clrs(1,:), 'LineWidth', linewidth);
        addlegendentry(hAx2, 'Measured');
        
        xlim(hAx2, xlims);
        xlabel(hAx2, 'Frequency [GHz]');
        ylabel(hAx2, '|S_{11}|^2+|S_{12}|^2');
    end
        
    
    %% 1-to-16
    if(ifile == find(strcmp(files, '16.s2p')))
        [hFig2, hAx2] = figureex;
        alignplot(hFig2, nplotx, nploty, hFig2.Number, [], 1);
%         title(hAx2, '16');
        
        
        folder = 'h:\Server\Exports\';
        filename = '1-to-16-lossy-back2back-';
        [f, ~, S11] = HFSS.LoadData([folder, filename, 'S11.csv']);
%             plot(hAx, f, S11, 'o', 'Color', clrs(1,:), 'LineWidth', linewidth);
%             addlegendentry(hAx, 'Simulation');
        [f, ~, S12] = HFSS.LoadData([folder, filename, 'S12.csv']);
% %             plot(hAx, f, S12, 'x', 'Color', clrs(2,:), 'LineWidth', linewidth);
%             plot(hAx2, f, 10*log10(10.^(S12./20).^2 + 10.^(S11./20).^2), 'x', 'Color', clrs(1,:), 'LineWidth', linewidth);
%             addlegendentry(hAx2, 'Simulation');
        S11i = interp1(fCST./1e9, squeeze(SCST(1,1,:)), f);
        S12i = interp1(fCST./1e9, squeeze(SCST(1,2,:)), f);
            plot(hAx2, f, 10*log10(10.^(S12./20).^2 + 10.^(S11./20).^2) + 2*10*log10(abs(S11i).^2+abs(S12i).^2), 'o', 'Color', clrs(1,:), 'LineWidth', linewidth);
            addlegendentry(hAx2, 'Simulation');
        
        
        plot(hAx2, fs, 10*log10(abs(s11).^2+abs(s12).^2), '-', 'Color', clrs(1,:), 'LineWidth', linewidth);
        addlegendentry(hAx2, 'Measured');
        xlim(hAx2, xlims);
        xlabel(hAx2, 'Frequency [GHz]');
        ylabel(hAx2, '|S_{11}|^2+|S_{12}|^2');
    end
    
    if(ifile == find(strcmp(files, '64.s2p')))
        [hFig2, hAx2] = figureex;
        plot(hAx2, fs, 10*log10(abs(s11).^2+abs(s12).^2), '-', 'Color', clrs(1,:), 'LineWidth', linewidth);
        addlegendentry(hAx2, 'Measured');
    end
        
        
        xlabel(hAx, 'Frequency [GHz]');
        ylabel(hAx, 'S_{ij} [dB]');
        xlim(hAx, xlims);
        ylim(hAx, [-30 0]);
        movelegend(hAx, 'se');
        hAx.LineWidth = axlinewidth;
    dispex('%s: %.2f at 12, %.2f at 32.\n', file, 20*log10(s12(fs == 12)), 20*log10(s12(fs == 32)));
end

for(ifig = 1:100)
    if(ishandle(ifig))
        hFig = figure(ifig);
        hAx = hFig.CurrentAxes;
        cropfigure(hFig);
        legendlinelength(hFig, 15);
        legend(hAx, 'NumColumns', 5);
        drawnow;
        movelegend(hFig, 'se');
    end
end

