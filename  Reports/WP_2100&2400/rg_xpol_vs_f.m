clearvars -except Eth Eph;
close all;

path = sprintf('%s/ Reports/WP_2100&2400/ra_circular/', resultdir_matlab);

taus = pi/180 * -45;%(0:10:90);

fs = 9:33;
f0s = [nan 30];
% ths = 0:30:60;
phs = 45;%0:45:90;
ths = 60;%[0 60];
% phs = 0;
ifig = 0;
for(tau = taus)
    for(f0 = f0s)
        fi = find(f0s == f0);
        for(ph = phs)
            phi = find(phs == ph);
            for(th = ths)
                thi = find(ths == th);
                load([path, '0.127_ff_', num2str(th), '_', num2str(ph), '.mat']);

                % Fields for both slots.
                EthX = Eth(1,:);
                EphX = Eph(1,:);
                EthY = Eth(2,:);
                EphY = Eph(2,:);

                % Get the value at f0.
                EthX0 = interp1(fs, EthX, f0);
                EphX0 = interp1(fs, EphX, f0);
                EthY0 = interp1(fs, EthY, f0);
                EphY0 = interp1(fs, EphY, f0);

    %             ax = 1;
    %             
    %             % Eo->0
    %             ay = -ax .* (EphX0 .* cos(tau) - EthX0 .* sin(tau)) ./ (EphY0 .* cos(tau) - EthY0 .* sin(tau));
    %             
    %             % Ep->0
    % %             ay = -ax .* (EthX0 .* cos(tau) + EphX0 .* sin(tau)) ./ (EthY0 .* cos(tau) + EphY0 .* sin(tau));

                ay = 1;
                % Eo->0
                if(isnan(f0))
                    ax = 0;
                else
                    ax = -ay .* (EphY0 .* cos(tau) - EthY0 .* sin(tau)) ./ (EphX0 .* cos(tau) - EthX0 .* sin(tau));
                end

                if(abs(ax) > 1)
                    ay = ay ./ ax;
                    ax = ax ./ ax;
                end

                Ep = (ax .* EthX + ay .* EthY) .* cos(tau) + (ax .* EphX + ay .* EphY) .* sin(tau);
                Eo = (ax .* EthX + ay .* EthY) .* sin(tau) - (ax .* EphX + ay .* EphY) .* cos(tau);

                [hFig, hAx] = a_fig(1);
    %                 plot(fs, abs(Ep));
    %                 plot(fs, abs(Eo));
                    plot(hAx, fs, abs(Eo ./ Ep));
                    xlabel(hAx, 'Frequency [GHz]');
                    ylabel(hAx, '|Eo/Ep|');
%                     legend({'Ep', 'Eo'});
                    hFig.Name = sprintf('theta = %f, phi = %f, f0 = %f', th, ph, f0);
                    xlim(hAx, [12 32]);
                    ylim(hAx, [0 1]);
                [hFig, hAx] = a_fig(2);
                    plot(hAx, fs, 20*log10(abs(Eo ./ Ep)));
                    xlabel(hAx, 'Frequency [GHz]');
%                     ylabel(hAx, '|Eo/Ep| [dB]');
                    ylabel(hAx, 'Xpol [dB]');
                    xlim(hAx, [12 32]);
                    ylim(hAx, [-30 0]);
            end
        end
    end
end

