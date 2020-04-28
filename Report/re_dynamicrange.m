clearvars -except Eth Eph;
close all;

project = CST.Application.Active3D();
fdsolver = project.FDSolver();

path = [cd, '\PhD-Matlab\Report\ra_circular\'];

ARs = zeros(13, 25);

fs = 9:33;
ports = 1:2;
f0s = [14 30];
ths = 0:5:60;
phs = 0:45:90;
ifig = 0;
for(f0 = f0s)
    fi = find(f0s == f0);
    for(ph = phs)
        phi = find(phs == ph);
        for(th = ths)
            thi = find(ths == th);
            for(error = -20:10:20)
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
                
                ays(fi, phi, thi) = ay;
            end
        end
    end
end

for(fi = 1:length(f0s))
    [hFig, hAx] = figureex(3+fi);
        plot(hAx, ths, 10*log10(abs(squeeze(ays(fi,1,:)))), '-');
        plot(hAx, ths, 10*log10(abs(squeeze(ays(fi,2,:)))), '-');
        plot(hAx, ths, 10*log10(abs(squeeze(ays(fi,3,:)))), '-');
        legend(hAx, strcat({'\phi = '}, strsplit(num2str(phs))))
        xlabel(hAx, '\theta [\circ]');
        ylabel(hAx, '|a_{y}| [dB]');
        hFig.Name = ['f0 = ', num2str(f0s(fi))];
        alignplot(hFig, 6, 4, hFig.Number, [], 1);
    [hFig, hAx] = figureex(9+fi);
        plot(hAx, ths, abs(squeeze(ays(fi,1,:))), '-');
        plot(hAx, ths, abs(squeeze(ays(fi,2,:))), '-');
        plot(hAx, ths, abs(squeeze(ays(fi,3,:))), '-');
        legend(hAx, strcat({'\phi = '}, strsplit(num2str(phs))))
        xlabel(hAx, '\theta [\circ]');
        ylabel(hAx, '|a_{y}|');
        hFig.Name = ['f0 = ', num2str(f0s(fi))];
        alignplot(hFig, 6, 4, hFig.Number, [], 1);
end

