function kb_1to8_hfss_plot(folder, filenamebase)

extension = '.csv';
nplotsX = 6;
nplotsY = 4;
linewidth = 1.5;
axlinewidth = 1;

figureex(1);
figureex(2);
figureex(3);
figureex(4);


%% S11
filename = [filenamebase, 'S11', extension];
[f, ~, S11] = HFSS.LoadData([folder, filename]);
[hFig, hAx] = figureex(1);
    if(length(hAx.Children) < 1)
        patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    
    plot(hAx, f, S11, 'LineWidth', linewidth);

    alignplot(hFig, nplotsX, nplotsY, hFig.Number, [], 1);
    xlabel(hAx, 'Frequency [GHz]');
    ylabel(hAx, '|S_{11}| [dB]');
    hFig.Name = filename;
    hAx.LineWidth = axlinewidth;
    ylim([-30 0]);

%% S12 for back2back
if(contains(lower(filenamebase), 'back2back'))
    filename = [filenamebase, 'S12', extension];
    [f, ~, S12] = HFSS.LoadData([folder, filename]);
    [hFig, hAx] = figureex(2);
        if(length(hAx.Children) < 1)
            patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end

        plot(hAx, f, S12, 'LineWidth', linewidth);

        alignplot(hFig, nplotsX, nplotsY, hFig.Number, [], 1);
        xlabel(hAx, 'Frequency [GHz]');
        ylabel(hAx, '|S_{12}| [dB]');
        hFig.Name = filename;
        hAx.LineWidth = axlinewidth;
        ylim([-1.5 0]);
        
        
    val = (10.^(S11./20)).^2 + (10.^(S12./20)).^2;
    [hFig, hAx] = figureex(3);
        if(length(hAx.Children) < 1)
            patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
        
        plot(hAx, f, 10*log10(val), 'LineWidth', linewidth);
        
        alignplot(hFig, nplotsX, nplotsY, hFig.Number, [], 1);
        xlabel(hAx, 'Frequency [GHz]');
        ylabel(hAx, '|S_{11}|^2 + |S_{12}|^2 [dB]');
        hFig.Name = filename;
        hAx.LineWidth = axlinewidth;
        ylim([-1.5 0]);
        
    return;
end

%% Sx1
filename = [filenamebase, 'Sx1', extension];
[f, ~, Sx1] = HFSS.LoadData([folder, filename]);
[hFig, hAx] = figureex;
    if(length(hAx.Children) < 1)
        patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    
    plot(hAx, f, Sx1, 'LineWidth', linewidth);

    alignplot(hFig, nplotsX, nplotsY, hFig.Number, [], 1);
    xlabel(hAx, 'Frequency [GHz]');
    ylabel(hAx, '|S_{11}| [dB]');
    hFig.Name = filename;
    hAx.LineWidth = axlinewidth;
    
    mx = max(Sx1(:));
    mn = min(Sx1(:));
    
    ylim([floor(mn) ceil(mx)]);
    ylim([ceil(mx)-2 ceil(mx)]);

%% S11 + Sx1
Sx12 = (10.^(Sx1./20)).^2;
S112 = (10.^(S11./20)).^2;
val = sum(Sx12, 2) + S112;
[hFig, hAx] = figureex(3);
    if(length(hAx.Children) < 1)
        patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end

    plot(hAx, f, 10*log10(val), 'LineWidth', linewidth);

    alignplot(hFig, nplotsX, nplotsY, hFig.Number, [], 1);
    xlabel(hAx, 'Frequency [GHz]');
    ylabel(hAx, '\Sigma_{i}|S_{1i}|^2 [dB]');
    hFig.Name = filename;
    hAx.LineWidth = axlinewidth;
    ylim([-1.5 0]);


%% Sx1 phase
filename = [filenamebase, 'Sx1-Phase', extension];
[f, ~, phase] = HFSS.LoadData([folder, filename]);

for(i = 1:size(phase, 2))
    for(fi = 1:size(phase, 1))
        if((phase(fi, i) - phase(fi, 1)) > 300)
            phase(fi, i) = phase(fi, i) - 360;
        elseif((phase(fi, i) - phase(fi, 1)) < -300)
            phase(fi, i) = phase(fi, i) + 360;
        end
    end
end

phase = max(phase, [], 2) - min(phase, [], 2);

[hFig, hAx] = figureex(4);
    if(length(hAx.Children) < 1)
        patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    
    plot(hAx, f, phase, 'LineWidth', linewidth);
    drawnow;

    alignplot(hFig, nplotsX, nplotsY, hFig.Number, [], 1);
    xlabel(hAx, 'Frequency [GHz]');
    ylabel(hAx, 'Max phase difference [\circ]');
    hFig.Name = filename;
    hAx.LineWidth = axlinewidth;
    ylim([0 15]);
%     ylim([0 5*ceil(max(values)/5)]);
    
%     mx = max(values(:));
%     mn = min(values(:));
%     
%     ylim([floor(mn), ceil(mx)]);
