%% Simulation setup.
df = 0.2;
fs = 1e9 * (df:df:40);
ths = [eps 60 60] * pi/180;
phs = [0   0  90] * pi/180;

% ths = [eps 10:10:60] * pi/180;
% phs = [0   0 0 0 0 0 0] * pi/180;
% phs = [0   90 90 90 90 90 90] * pi/180;
% phs = [0   45 45 45 45 45 45 45] * pi/180;

%% Perform simulation.
% Impedance of feed capacitance.    
Zcap = 1 ./ (1j .* 2.*pi.*fs .* C);

fig = figureex(1);
alignplot(fig, 1, 1, 1, [], 2);

array = InfiniteArray(slot, tlineup, tlinedown);

linewidth = 1;
axlinewidth = 1;

styles = {'k-', 'r--', 'b:'};
styles = {'-', '-', '-', '-', '-', '-', '-', '-', '-'};
for(iangle = 1:length(ths))
    th = ths(iangle);
    ph = phs(iangle);
    
    % Open figures.
    figGamma = figureex(1);%+3*(iangle-1));
        axGamma = figGamma.CurrentAxes;
        alignplot(figGamma, 6, 4, figGamma.Number, [], 1);
    figReal = figureex(2);%+3*(iangle-1));
        axReal = figReal.CurrentAxes;
        alignplot(figReal, 6, 4, figReal.Number, [], 1);
        axReal.ColorOrder = reshape(repmat(lines(7), 1, 2).', [], 14).';
    figVSWR = figureex(3);%+3*(iangle-1));
        axVSWR = figVSWR.CurrentAxes;
        alignplot(figVSWR, 6, 4, figVSWR.Number, [], 1);
%     figImag = figureex(3+3*(iangle-1));
%         axImag = figImag.CurrentAxes;
%         alignplot(figImag, 6, 4, figImag.Number, [], 1);
        
%     axGamma = subplot(3, 3, 1+3*(iangle-1));
%         hold on; box on; grid on;
%         title(axGamma, sprintf('\\Gamma, \\theta = %.0f, \\phi = %.0f', th * 180/pi, ph * 180/pi));
%     axReal = subplot(3, 3, 2+3*(iangle-1));
%         hold on; box on; grid on;
%         title(axReal, sprintf('Re(Z), \\theta = %.0f, \\phi = %.0f', th * 180/pi, ph * 180/pi));
%     axImag = subplot(3, 3, 3+3*(iangle-1));
%         hold on; box on; grid on;
%         title(axImag, sprintf('Im(Z), \\theta = %.0f, \\phi = %.0f', th * 180/pi, ph * 180/pi));

    axImag = axReal;
    
    % If the axes are empty, put in the frequency range indicator.
    if(length(axGamma.Children) < 1)
%         plot(axGamma, [12 12 nan 31 31], [-1e9 1e9 nan -1e9 1e9], 'k');
%         patch(axGamma, [13.75 14.5 14.5 13.75], [-10+3*(iangle>1) -10+3*(iangle>1) 0 0], [0 0 0], ...
%             'FaceAlpha', 0.1, 'EdgeColor', 'none');
%         patch(axGamma, [28 31 31 28], [-10+3*(iangle>1) -10+3*(iangle>1) 0 0], [0 0 0], ...
%             'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(axGamma, [13.75 14.5 14.5 13.75], [-30 -30 0 0], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(axGamma, [28 31 31 28], [-30 -30 0 0], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        patch(axReal, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(axReal, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        
        patch(axVSWR, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(axVSWR, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
%         patch(axImag, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
%             'FaceAlpha', 0.1, 'EdgeColor', 'none');
%         patch(axImag, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
%             'FaceAlpha', 0.1, 'EdgeColor', 'none');

        xlabel(axGamma, 'Frequency [GHz]');
        xlabel(axReal, 'Frequency [GHz]');
        xlabel(axVSWR, 'Frequency [GHz]');
        ylabel(axGamma, '|\Gamma| [dB]');
        ylabel(axReal, 'Active input impedance [\Omega]');
        ylabel(axVSWR, 'VSWR');
        axGamma.LineWidth = axlinewidth;
        axReal.LineWidth = axlinewidth;
        axVSWR.LineWidth = axlinewidth;
    end
    
    % Calculate input impedance.
    tc = tic;
    Zas = array.GetInputImpedance(fs, th, ph);
%     dispex('Calculated input impedance in %.3fs.\n', toc(tc));
    ZasC = Zas + Zcap; % With series capacitance.
    
    % Calculate reflection coefficient.
    Gamma = (ZasC - zfeed) ./ (ZasC + zfeed);
    VSWR = (1 + abs(Gamma)) ./ (1 - abs(Gamma));
    
    % Output worst values.
%     ind = find(20*log10(abs(Gamma)) == max(20*log10(abs(Gamma(fs >= 13e9 & fs <= 31e9)))));
    ind = find(20*log10(abs(Gamma)) == max(20*log10(abs(Gamma(fs >= 13.75e9 & fs <= 14.5e9 | fs >= 28e9 & fs <= 31e9)))));
    dispex('Worst at %02.0f,%02.0f is %.2f, f = %.1fGHz, Z = %.2f + %.2fj.\n', ...
        th*180/pi, ph*180/pi, 20*log10(abs(Gamma(ind))), fs(ind)/1e9, real(ZasC(ind)), imag(ZasC(ind)));
    
    % Plot results.
    % Reflection coefficient.
    plot(axGamma, fs./1e9, 10*log10(abs(Gamma).^2), styles{iangle}, 'LineWidth', linewidth);
        ylim(axGamma, [-30 0]);
    % Input impedance.
    plot(axReal, fs./1e9, real(ZasC), styles{iangle}, 'LineWidth', linewidth);
        ylim(axReal, [0 200]);
    plot(axImag, fs./1e9, imag(ZasC), '--', 'LineWidth', linewidth);%styles{iangle});
        ylim(axImag, [-100 200]);
    % VSWR
    plot(axVSWR, fs./1e9, VSWR, 'LineWidth', linewidth);%styles{iangle});
        ylim(axVSWR, [1 8]);
        addlegendentry([num2str(round(th * 180/pi)), '\circ']);
%         legend(axVSWR, 'NumColumns', 7);
        
    drawnow;
    
%     figureex(iangle+10);
%         plot(fs/1e9, real(Zas), 'k');
%         plot(fs/1e9, imag(Zas), 'k--');
end
% legend('B', 'H', 'E');
