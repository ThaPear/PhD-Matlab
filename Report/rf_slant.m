clearvars -except Eth Eph;
close all;

path = [cd, '\PhD-Matlab\Report\ra_circular\'];

taus = pi/180 * (0:15:90);

fs = 9:33;
ports = 1:2;
f0s = [14 30];
% ths = 0:30:60;
phs = 45;%0:45:90;
ths = 60;%[0 60];
% phs = 0;
ifig = 0;
for(itau = 1:length(taus))
    tau = taus(itau);
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
                ax = -ay .* (EphY0 .* cos(tau) - EthY0 .* sin(tau)) ./ (EphX0 .* cos(tau) - EthX0 .* sin(tau));

                if(abs(ax) > 1)
                    ay = ay ./ ax;
                    ax = ax ./ ax;
                end

                Ep = (ax .* EthX + ay .* EthY) .* cos(tau+eps) + (ax .* EphX + ay .* EphY) .* sin(tau+eps);
                Eo = (ax .* EthX + ay .* EthY) .* sin(tau+eps) - (ax .* EphX + ay .* EphY) .* cos(tau+eps);

                [hFig, hAx] = a_fig(if0*2-1);
    %                 plot(fs, abs(Ep));
    %                 plot(fs, abs(Eo));
                    plot(fs, abs(Eo ./ Ep));
                    xlabel('Frequency [GHz]');
                    ylabel('|E|');
%                     legend({'Ep', 'Eo'});
                    hFig.Name = sprintf('theta = %g, phi = %g, f0 = %g', th, ph, f0);
                    alignplot(hFig, 8, 4, hFig.Number, [], 2);
                [hFig, hAx] = a_fig(if0*2);
                    plot(fs, 20*log10(abs(Eo ./ Ep)+eps));
                    alignplot(hFig, 8, 4, hFig.Number, [], 2);
                    xlabel('Frequency [GHz]');
                    ylabel('Xpol [dB]')
                    xlim([12 32]);
                    ylim([-30 0]);
            end
        end
    end
end
[hFig, hAx] = figureex(2);
legend(hAx, hAx.Children(length(taus):-1:1), strcat({'\alpha = '}, strsplit(num2str(round(taus*180/pi))), '\circ'), ...
    'NumColumns', 2);
legendlinelength(hAx, 25);
drawnow; movelegend(hAx, 'n');

[hFig, hAx] = figureex(4);
legend(hAx, hAx.Children(length(taus):-1:1), strcat({'\alpha = '}, strsplit(num2str(round(taus*180/pi))), '\circ'), ...
    'NumColumns', 2);
legendlinelength(hAx, 25);
drawnow; movelegend(hAx, 'n');