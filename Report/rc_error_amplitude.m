clearvars -except Eth Eph;
close all;

path = [cd, '\PhD-Matlab\Report\ra_circular\'];

ARs = zeros(13, 25);

fs = 9:33;
ports = 1:2;
f0s = [14 30];
ths = 0:60:60;
phs = 0:45:90;
errors = [0:-0.1:-0.3 0.1:0.1:0.4];
ifig = 0;
for(f0 = f0s)
    for(ph = phs)
        for(th = ths)
            ifig = ifig + 2;
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
                ay = ay * (1+error);

                Er = ((ax .* EthX + ay .* EthY) + 1j .* (ax * EphX + ay * EphY)) ./ sqrt(2);
                El = ((ax .* EthX + ay .* EthY) - 1j .* (ax * EphX + ay * EphY)) ./ sqrt(2);

                tau = angle(Er ./ El);

                AR = (abs(Er) + abs(El)) ./ (abs(Er) - abs(El));
                ARs(th/5+1,:) = AR;
                %% Plot Positive error
                if(error+eps >= 0)
                [hFig, hAx] = a_fig(ifig);
                    hFig.NumberTitle = 0;
                    plot(hAx, fs, 20*log10(abs(AR)), 'LineWidth', 1);
                    addlegendentry(sprintf('%+03.0f%%', error*100));
%                     if(f0 < 20)
%                         itxt = find(20*log10(AR) < 6, 1);
%                         text(fs(itxt), 20*log10(AR(itxt)), [num2str(error*100, '%+02.0f'), '%  '], 'HorizontalAlignment', 'right');
%                     else
%                         itxt = find(20*log10(AR) < 6, 1, 'last');
%                         text(fs(itxt), 20*log10(AR(itxt)), ['  ', num2str(error*100, '%+02.0f'), '%'], 'HorizontalAlignment', 'left');
%                     end
        %             ylim(hAx, [0 max(20*log10(abs(AR)))]);
                    ylim(hAx, [0 6]);
                    hFig.Name = sprintf('Phi = %i, f0 = %g, th = %i', ph, f0, th);
                    ylabel(hAx, 'Axial Ratio [dB]');
                end
                %% Plot Negative error
                if(error-eps <= 0)
                [hFig, hAx] = a_fig(ifig-1);
                    hFig.NumberTitle = 0;
                    plot(hAx, fs, 20*log10(abs(AR)), 'LineWidth', 1);
                    addlegendentry(sprintf('%+03.0f%%', error*100));
%                     if(f0 < 20)
%                         itxt = find(20*log10(AR) < 6, 1);
%                         text(fs(itxt), 20*log10(AR(itxt)), [num2str(error*100, '%+02.0f'), '%  '], 'HorizontalAlignment', 'right');
%                     else
%                         itxt = find(20*log10(AR) < 6, 1, 'last');
%                         text(fs(itxt), 20*log10(AR(itxt)), ['  ', num2str(error*100, '%+02.0f'), '%'], 'HorizontalAlignment', 'left');
%                     end
        %             ylim(hAx, [0 max(20*log10(abs(AR)))]);
                    ylim(hAx, [0 6]);
                    hFig.Name = sprintf('Phi = %i, f0 = %g, th = %i', ph, f0, th);
                    ylabel(hAx, 'Axial Ratio [dB]');
                end
            end
        end
    end
end

