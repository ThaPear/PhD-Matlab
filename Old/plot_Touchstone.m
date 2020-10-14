SetupPath;
close all;
clear;

path = GetFullPath('Temp/');

% filebase = '0.127_feed_schem';
% fileext = '.s2p';
% runIDs = 1:2;

filebase = '80_dual_below';
fileext = '.s4p';
runIDs = 0;
C = inf;

filebase = '0.127_feed_below';
fileext = '.s2p';
runIDs = 1:3;
C = inf;

% path = 'e:\ Simulations\ Unit Cell Design Path\Export 01\';
% filebase = '01 - 0.127';
% fileext = '.s2p';
% % runIDs = [38 20:25]; % Plane
% % runIDs = [38 26:31]; % Plane
% % runIDs = [38 32:37]; % D-plane
% runIDs = [38 37 25]; % Broadside 60 60
% C = 0.2e-12;

% path = 'e:\ Simulations\ Unit Cell Design Path\Export 06\';
% filebase = '06 - 0.127_feed';
% fileext = '.s2p';
% runIDs = [1 2 3];
% C = inf;

path = 'e:\ Simulations\ Unit Cell Design Path\Export 08\';
filebase = '08 - 0.127_feed_below';
fileext = '.s2p';
runIDs = [1 2 3];
C = inf;

linewidth = 1;
axlinewidth = 1;

for(id = runIDs)
    if(id == 0)
        filename = [path, filebase, fileext];
    else
        filename = [path, filebase, '_', num2str(id), fileext];
    end
    [parameters, S] = CST.LoadData(filename);
    f = parameters.frequencies;
    
    S = struct('s11', squeeze(S(1,1,:)).', ...
               's12', squeeze(S(1,2,:)).', ...
               's21', squeeze(S(2,1,:)).', ...
               's22', squeeze(S(2,2,:)).');
           
           
    % Add capacitance.
    Zin1 = (1+S.s11)./(1-S.s11).*parameters.slot_impedance;
    Zin2 = (1+S.s22)./(1-S.s22).*parameters.slot_impedance;
           
    Zcap = 1 ./ (1j .* 2.*pi.*f .* C);
    Zin1 = Zin1 + Zcap;
    Zin2 = Zin2 + Zcap;
    S.s11 = (Zin1 - parameters.slot_impedance) ./ (Zin1 + parameters.slot_impedance);
    S.s22 = (Zin2 - parameters.slot_impedance) ./ (Zin2 + parameters.slot_impedance);

	VSWR1 = S2VSWR(S.s11);
	VSWR2 = S2VSWR(S.s22);    

    switch(fileext)
        case '.s2p'
            [figS, axS] = figureex(1);
                axS.LineWidth = axlinewidth;
                plot(axS, f/1e9, 20*log10(abs(S.s11)), 'LineWidth', linewidth);
            [figS, axS] = figureex(2);
                axS.LineWidth = axlinewidth;
                plot(axS, f/1e9, 20*log10(abs(S.s22)), 'LineWidth', linewidth);
            [figVSWR, axVSWR] = figureex(3);
                axVSWR.LineWidth = axlinewidth;
                plot(axVSWR, f/1e9, VSWR1, 'LineWidth', linewidth);
            [figVSWR, axVSWR] = figureex(4);
                axVSWR.LineWidth = axlinewidth;
                plot(axVSWR, f/1e9, VSWR2, 'LineWidth', linewidth);
        case '.s4p'
            [figS, axS] = figureex;
                axS.LineWidth = axlinewidth;
                plot(axS, f/1e9, 20*log10(abs(S.s11)), 'LineWidth', linewidth);
                plot(axS, f/1e9, 20*log10(abs(S.s22)), 'LineWidth', linewidth);
                plot(axS, f/1e9, 20*log10(abs(S.s33)), 'LineWidth', linewidth);
                plot(axS, f/1e9, 20*log10(abs(S.s44)), 'LineWidth', linewidth);
    end
    
    
    
    
end
for(i = 1:100)
    if(~ishandle(i))
        break;
    end
    figS = figureex(i); axS = figS.CurrentAxes;
    xlabel('Frequency [GHz]');
    xlim([12 32]);
    if(i < 3)
        ylabel('|\Gamma| [dB]');
        ylim([-30 0]);
    else
        ylabel('VSWR');
        ylim([1 8]);
    end
    patch(axS, [28 31 31 28], [-30 -30 10 10], [0 0 0], ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none');
    patch(axS, [13.75 14.5 14.5 13.75], [-30 -30 10 10], [0 0 0], ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    alignplot(figS, 6, 4, figS.Number, [], 1);
    
    if(i == 1 || i == 3)
        title('X-oriented Slot');
        legend({'Broadside', 'H-plane 60\circ', 'E-plane 60\circ'});
        movelegend('se');
        drawnow;
        cropfigure;
        drawnow;
        movelegend('se');
    elseif(i == 2 || i == 4)
        title('Y-oriented Slot');
        legend({'Broadside', 'E-plane 60\circ', 'H-plane 60\circ'});
        movelegend('se');
        drawnow;
        cropfigure;
        drawnow;
        movelegend('se');
    end
end