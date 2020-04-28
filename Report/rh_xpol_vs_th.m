clearvars -except Eth Eph;
close all;

path = [cd, '\PhD-Matlab\Report\ra_circular\'];

taus = pi/180 * -45;%(0:10:90);

fs = 9:33;
f0s = nan;
% ths = 0:30:60;
phs = 45;%0:45:90;
ths = [0:5:80 84 85 89];
% phs = 0;
ifig = 0;

Eps = zeros(length(f0s), length(ths), length(phs), length(fs));
Eos = zeros(length(f0s), length(ths), length(phs), length(fs));

for(tau = taus)
    for(if0 = 1:length(f0s))
        f0 = f0s(if0);
        for(iph = 1:length(phs))
            ph = phs(iph);
            for(ith = 1:length(ths))
                th = ths(ith);
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
                
                Eps(if0, ith, iph, :) = Ep;
                Eos(if0, ith, iph, :) = Eo;
            end
        end
    end
end
f0s = [29];
for(f0 = f0s)
    [hFig, hAx] = figureex(1);
        plot(hAx, ths, abs(Eos(1, :, 1, find(fs == f0)) ./ Eps(1, :, 1, find(fs == f0))), 'LineWidth', 1);
        hAx.LineWidth = 1;
        xlabel(hAx, 'Theta [\circ]');
        ylabel(hAx, '|Eo/Ep|');
        hFig.Name = sprintf('theta = %f, phi = %f, f0 = %f', th, ph, f0);
        xlim(hAx, [0 90]);
        ylim(hAx, [0 1]);
        alignplot(hFig, 8, 4, hFig.Number, [], 2);
    [hFig, hAx] = figureex(2);
        plot(hAx, ths, 20*log10(abs(Eos(1, :, 1, find(fs == f0)) ./ Eps(1, :, 1, find(fs == f0)))), 'LineWidth', 1);
        hAx.LineWidth = 1;
        xlabel(hAx, 'Theta [\circ]');
    %     ylabel(hAx, '|Eo/Ep| [dB]');
        ylabel(hAx, 'Xpol [dB]');
        xlim(hAx, [0 90]);
        alignplot(hFig, 8, 4, hFig.Number, [], 2);
        ylim(hAx, [-30 0]);
        alignplot(hFig, 8, 4, hFig.Number, [], 2);
end