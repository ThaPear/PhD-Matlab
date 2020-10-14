clearvars -except Eth Eph;
close all;

path = sprintf('%s/ Reports/WP_2100&2400/ra_circular/', resultdir_matlab);

ARs = zeros(13, 25);

horizontal = 1;
if(horizontal)
    hAxes = TightSubplot(4, 2, [0.005 0.03], [0.05 0.01], [0.12 0.06]);
    hAxes = hAxes([1 5 2 6 3 7 4 8]);
    hFig = gcf;
    alignplot(hFig, 4, 3.8, 1, [], 2);
    
    titles = {'Broadside', 'H-plane 60\circ', 'D-plane 60\circ', 'E-plane 60\circ'};
else
    hAxes = TightSubplot(2, 4, [0.01 0.02], [0.05 0.01], [0.06 0.01]);
    hFig = gcf;
    alignplot(hFig, 4, 2, 1, [], 2);
end


fs = 9:33;
ports = 1:2;
f0s = [14 30];
ths = 0:60:60;
phs = 0:45:90;
errors = -20:10:20;[0:5:20];%[0:-5:-10 5:5:10];
ifig = 0;
for(ph = phs)
    for(th = ths)
        for(f0 = f0s)
            % Only do broadside once per f0.
            if(th == 0 && ph > 0)
                continue;
            end
            ifig = ifig + 1;
            for(error = errors)
                load([path, '0.127_ff_', num2str(th), '_', num2str(ph), '.mat']);

                EthX = Eth(1,:);
                EphX = Eph(1,:);
                EthY = Eth(2,:);
                EphY = Eph(2,:);

                EthX0 = interp1(fs, EthX, f0);
                EphX0 = interp1(fs, EphX, f0);
                EthY0 = interp1(fs, EthY, f0);
                EphY0 = interp1(fs, EphY, f0);

                ax = 1;
                ay = -ax * (EthX0 - 1j * EphX0) / (EthY0 - 1j * EphY0);
                ay = ay * exp(1j * error * pi/180);

                Er = ((ax .* EthX + ay .* EthY) + 1j .* (ax * EphX + ay * EphY)) ./ sqrt(2);
                El = ((ax .* EthX + ay .* EthY) - 1j .* (ax * EphX + ay * EphY)) ./ sqrt(2);

                tau = angle(Er ./ El);

                AR = (abs(Er) + abs(El)) ./ (abs(Er) - abs(El));
                ARs(th/5+1,:) = AR;
                %% Plot AR for this error.
                hAx = hAxes(ifig);
                    grid(hAx, 'on'); box(hAx, 'on'); hold(hAx, 'on');
                    hAx.LineWidth = 1;
                    set(hAx, 'DefaultLineLineWidth', 1);
                    if(length(hAx.Children) < 2)
                        patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                            'FaceAlpha', 0.1, 'EdgeColor', 'none');
                        patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                            'FaceAlpha', 0.1, 'EdgeColor', 'none');
                    end
%                 [hFig, hAx] = a_fig(ifig);
                    plot(hAx, fs, 20*log10(abs(AR)), 'LineWidth', 1);
                    xlim(hAx, [12 32]);
                    ylim(hAx, [0 6]);
                    
                    if(horizontal && mod(ifig, 2) ~= 0)
                        title(hAx, titles{ceil(ifig/2)})
                    end
                    if(ifig > 6 && ~horizontal || mod(ifig, 2) == 0 && horizontal)
                        xlabel(hAx, 'Frequency [GHz]');
                        xticks(hAx, 5:5:35);
                        xticklabels(hAx, 'auto');
                    end
                    if(mod(ifig, 2) ~= 0 && ~horizontal || ifig < 3 && horizontal)
                        ylabel(hAx, 'Axial Ratio [dB]');
                        yticks(hAx, 0:2:6);
                        yticklabels(hAx, 'auto');
                    end
                    if(ifig == 8 && ~horizontal || ifig == 7 && horizontal)
                        addlegendentry(hAx, sprintf('%+03.0f\\circ', error));
                        drawnow;
                        if(horizontal)
                            movelegend(hAx, 'se');
                        else
                            movelegend(hAx, 'nw');
                        end
                    end
            end
        end
    end
end

drawnow;
print -clipboard -dmeta;
dispex('Figure copied to clipboard.\n');

