clearvars -except Eth Eph;
close all;

path = sprintf('%s/ Reports/WP_2100&2400/%s/', resultdir_matlab, mfilename);

fs = 9:33;
ports = 1:2;
ifig = 0;
for(f0 = [14 30])
    for(ph = 0:45:90)
        ifig = ifig + 2;
        ths = [0:10:60];
%         ths = [0:5:80 84 89]
        for(th = ths)
            if(~exist([path, '0.127_ff_', num2str(th), '_', num2str(ph), '.mat'], 'file'))
                project = CST.Application.Active3D();
                fdsolver = project.FDSolver();
                
                project.StoreParameter('aa_theta', th);
                project.StoreParameter('aa_phi', ph);
                project.Rebuild();
                if(~fdsolver.Start()); dispex('Simulation failed.\n'); return; end
                [Eth, Eph] = CST.GetFarfieldResult(project, fs, ports, 'aa_theta', 'aa_phi');
                save([path, '0.127_ff_', num2str(th), '_', num2str(ph), '.mat'], 'Eth', 'Eph');
            else
                load([path, '0.127_ff_', num2str(th), '_', num2str(ph), '.mat']);
            end


            EthX = Eth(1,:);
            EphX = Eph(1,:);
            EthY = Eth(2,:);
            EphY = Eph(2,:);

    %         f0 = 30;
            EthX0 = interp1(fs, EthX, f0);
            EphX0 = interp1(fs, EphX, f0);
            EthY0 = interp1(fs, EthY, f0);
            EphY0 = interp1(fs, EphY, f0);

            ax = 1;
            % ay = 1j;
            ay = -ax * (EthX0 - 1j * EphX0) / (EthY0 - 1j * EphY0);

            Er = ((ax .* EthX + ay .* EthY) + 1j .* (ax * EphX + ay * EphY)) ./ sqrt(2);
            El = ((ax .* EthX + ay .* EthY) - 1j .* (ax * EphX + ay * EphY)) ./ sqrt(2);

            tau = angle(Er ./ El);

            AR = (abs(Er) + abs(El)) ./ (abs(Er) - abs(El));

            [hFig, hAx] = figureex(ifig);
                if(length(hAx.Children) < 2)
                    patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                        'FaceAlpha', 0.1, 'EdgeColor', 'none');
                    patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                        'FaceAlpha', 0.1, 'EdgeColor', 'none');
                end
                plot(hAx, fs, 20*log10(abs(AR)), 'LineWidth', 1);
                hAx.LineWidth = 1;
%                 itxt = find(20*log10(AR) < 6, 1);
%                 text(fs(itxt), 20*log10(AR(itxt)), [num2str(th), '  '], 'HorizontalAlignment', 'right');
    %             ylim(hAx, [0 max(20*log10(abs(AR)))]);
                ylim(hAx, [0 6]);
                hFig.Name = sprintf('Phi = %i, f0 = %g', ph, f0);
                alignplot(hFig, 8, 4, hFig.Number, [], 2);
                xlabel(hAx, 'Frequency [GHz]');
                ylabel(hAx, 'Axial Ratio [dB]');
            [hFig, hAx] = figureex(ifig-1);
                if(length(hAx.Children) < 2)
                    patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                        'FaceAlpha', 0.1, 'EdgeColor', 'none');
                    patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                        'FaceAlpha', 0.1, 'EdgeColor', 'none');
                end
                plot(hAx, fs, tau*180/pi, 'LineWidth', 1);
                hAx.LineWidth = 1;
                ylim(hAx, [-180 180]);
                hFig.Name = sprintf('Phi = %i, f0 = %g', ph, f0);
                alignplot(hFig, 8, 4, hFig.Number, [], 2);
                xlabel(hAx, 'Frequency [GHz]');
                ylabel(hAx, 'Tau [\circ]');
        end
        [hFig, hAx] = figureex(ifig);
            hLeg = legend(hAx, hAx.Children(length(ths):-1:1), strcat({'\theta = '}, strsplit(num2str(ths)), '\circ'), ...
                        'NumColumns', 1);%, 'Location', 'northoutside');
            legendlinelength(hAx, 10);
            movelegend(hAx, 'eneo');
            xlim(hAx, [12 32]);
    end
end

%% EuCAP2021
for(ifig = 1:100)
    if(ishandle(ifig))
        [hFig, hAx] = figureex(ifig);
        if(mod(ifig, 2) == 1 || ifig == 6)
            close(ifig);
            continue;
        elseif(ifig == 12)
            hLeg = legend(hAx, 'NumColumns', 4);
            legendlinelength(hAx, 15);
            drawnow; movelegend(hAx, 'no');
            hFig.Units = 'pixels';
            hAx.Units = 'pixels';
            hLeg.Units = 'pixels';
            drawnow;
            hFig.Position(3:4) = hLeg.Position(3:4) + [1 1];
            drawnow;
            hLeg.Position(1:2) = [1 1];
            continue;
        end
        alignplot(hFig, 8, 4, hFig.Number/2, [], 2);
        legend(hAx, 'off');
    end
end
        


