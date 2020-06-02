clearvars -except Eth Eph;
close all;

path = [cd, '\PhD-Matlab\ Reports\WP_2100&2400\ra_circular\'];

ARs = zeros(13, 25);

fs = 9:33;
ports = 1:2;
f0s = [14 16 30];
th0s = [25 25 15];
% th0s = [55 55 40];
ths = 0:5:80;
phs = [0 45]; %0:45:90; % 90 is identical to 0.
ifig = 0;
for(f0 = f0s)
    ifig = ifig + 3;
    th0 = th0s(f0s == f0);
    for(ph = phs)
        if(ph == 45 && f0 == 30)
            th0 = th0 + 25;
        end
        %% Get the weight for a single angle.
        load([path, '0.127_ff_', num2str(th0), '_', num2str(ph), '.mat']);
        EthX = Eth(1,:);
        EphX = Eph(1,:);
        EthY = Eth(2,:);
        EphY = Eph(2,:);
        
        EthX0 = interp1(fs, EthX, f0);
        EphX0 = interp1(fs, EphX, f0);
        EthY0 = interp1(fs, EthY, f0);
        EphY0 = interp1(fs, EphY, f0);

        ax = 1;
        % ay = 1j;
        ay = -ax * (EthX0 - 1j * EphX0) / (EthY0 - 1j * EphY0);
        
%             [hFig, hAx] = figureex(ifig);
%                 if(length(hAx.Children) < 2)
%                     patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
%                         'FaceAlpha', 0.1, 'EdgeColor', 'none');
%                     patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
%                         'FaceAlpha', 0.1, 'EdgeColor', 'none');
%                 end
%                 plot(hAx, fs, 20*log10(abs(AR)), 'o-');
%                 if(f0 < 20)
%                     itxt = find(20*log10(AR) < 6, 1);
%                     text(fs(itxt), 20*log10(AR(itxt)), [num2str(th), '  '], 'HorizontalAlignment', 'right');
%                 else
%                     itxt = find(20*log10(AR) < 6, 1, 'last');
%                     text(fs(itxt), 20*log10(AR(itxt)), ['  ', num2str(th)], 'HorizontalAlignment', 'left');
%                 end
%     %             ylim(hAx, [0 max(20*log10(abs(AR)))]);
%                 ylim(hAx, [0 6]);
%                 hFig.Name = sprintf('Phi = %i, f0 = %g, th0 = %i', ph, f0, th0);
%                 alignplot(hFig, 6, 4, hFig.Number, [], 1);
%             [hFig, hAx] = figureex(ifig-1);
%                 if(length(hAx.Children) < 2)
%                     patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
%                         'FaceAlpha', 0.1, 'EdgeColor', 'none');
%                     patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
%                         'FaceAlpha', 0.1, 'EdgeColor', 'none');
%                 end
%                 plot(hAx, fs, tau*180/pi, 'o-');
%                 ylim(hAx, [-180 180]);
%                 hFig.Name = sprintf('Phi = %i, f0 = %g, th0 = %i', ph, f0, th0);
%                 alignplot(hFig, 6, 4, hFig.Number, [], 1);
        
        ifig = ifig + 1;
        for(th = ths)
            load([path, '0.127_ff_', num2str(th), '_', num2str(ph), '.mat']);

            EthX = Eth(1,:);
            EphX = Eph(1,:);
            EthY = Eth(2,:);
            EphY = Eph(2,:);

            Er = ((ax .* EthX + ay .* EthY) + 1j .* (ax * EphX + ay * EphY)) ./ sqrt(2);
            El = ((ax .* EthX + ay .* EthY) - 1j .* (ax * EphX + ay * EphY)) ./ sqrt(2);

            tau = angle(Er ./ El);

            AR = (abs(Er) + abs(El)) ./ (abs(Er) - abs(El));
            ARs(th/5+1,:) = AR;
        end
        [hFig, hAx] = figureex(ifig);
            contours = 0:1:6;
            [hAxImage, hAxContour, hColorbar] = imagecontour(hFig, fs, ths, 20*log10(abs(ARs)), contours, 'pcolor', 'colorbar');
            hold(hAxImage, 'on'); plot(hAxImage, f0, th0, 'w+');
            plot(hAxImage, [13.75 13.75 nan 14.5 14.5 nan 28 28 nan 31 31], ...
                           [-1e3 1e3 nan -1e3 1e3 nan -1e3 1e3 nan -1e3 1e3], 'Color', [0.8 0.8 0.8]);
            caxis(hAxImage, [0 6]);
            shading(hAxImage, 'interp');
            grid(hAxContour, 'on');
            hAxContour.Color = 'none';
            colormap(hAxImage, 'jet');
            xlabel(hAxImage, 'Frequency [GHz]');
            ylabel(hAxImage, '\theta [\circ]');
            hFig.Name = sprintf('Phi = %i, f0 = %g, th0 = %i', ph, f0, th0);
            ylabel(hColorbar, 'Axial Ratio [dB]');
            
            alignplot(hFig, 8, 4, hFig.Number, [], 2);
            cropfigure(hFig, 10);
            drawnow;
            cropfigure(hFig);
%             pcolor(fs, ths, 20*log10(abs(ARs)));
%             alignplot(hFig, 6, 3, hFig.Number, [], 1);
%             colormap jet;
%             caxis([0 6]);
%             shading(hAx, 'interp');
    end
end

