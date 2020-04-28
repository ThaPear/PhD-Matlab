SetupPath;
clear;
close all;

c0 = Constants.c0;
z0 = Constants.z0;
f0 = 31e9;          % 13.75 to 31 GHz.

tlineup = FreeSpace();
tlinedown = ShortedLine(1, c0/f0/4);
% tlinedown = FreeSpace();

%% Simulation
fs = (12:0.1:32) * 1e9;
fs = f0;
ths = pi/180 * [eps, (1:90)];
phs = pi/180 * 45;

EL3cos = zeros(length(ths), length(phs), length(fs));
EL3xs = zeros(length(ths), length(phs), length(fs));

for(ith = 1:length(ths))
    th = ths(ith);
    for(iph = 1:length(phs))
        ph = phs(iph);

        % Calculate propagation constants.
        [k0, kx0, ky0, kz0] = k(fs, 1, th, ph);
        % [kd, ~, ~, kzd] = k(fs, (z0/z1)^2, th, ph);
        kr = sqrt(kx0.^2 + ky0.^2);

        % Calculate impedances.
        [z0, z0te, z0tm] = z(1, k0, kz0);

        % Calculate currents on slot plane.
        zupte = tlineup.GetInputImpedance(1, fs, k0, kr);
        zuptm = tlineup.GetInputImpedance(0, fs, k0, kr);

        zdownte = tlinedown.GetInputImpedance(1, fs, k0, kr);
        zdowntm = tlinedown.GetInputImpedance(0, fs, k0, kr);

        zte = 1 ./ (1 ./ zupte + 1 ./ zdownte);
        ztm = 1 ./ (1 ./ zuptm + 1 ./ zdowntm);

        ite = sin(ph);
        itm = cos(ph);

        vte = ite .* zte;
        vtm = itm .* ztm;
        
%         iteup = vte ./ zupte;
%         itmup = vtm ./ zuptm;
        
        iteup = ite .* zdownte ./ (zupte + zdownte);
        itmup = itm .* zdowntm ./ (zuptm + zdowntm);

        % Transform currents through ABCD to top of ADL.
        ABCDte = ABCDMatrix.identity();
        ABCDtm = ABCDMatrix.identity();

        % vtetop = (  ABCDte.D .* vte - ABCDte.B .* ite) ./ (ABCDte.A .* ABCDte.D - ABCDte.B .* ABCDte.C);
        % itetop = (- ABCDte.C .* vte + ABCDte.A .* ite) ./ (ABCDte.A .* ABCDte.D - ABCDte.B .* ABCDte.C);
        % 
        % vtmtop = (  ABCDtm.D .* vtm - ABCDtm.B .* itm) ./ (ABCDtm.A .* ABCDtm.D - ABCDtm.B .* ABCDtm.C);
        % itmtop = (- ABCDtm.C .* vtm + ABCDtm.A .* itm) ./ (ABCDtm.A .* ABCDtm.D - ABCDtm.B .* ABCDtm.C);

        vtetop = ABCDte.A .* vte + ABCDte.B .* iteup;
        itetop = ABCDte.C .* vte + ABCDte.D .* iteup;

        vtmtop = ABCDtm.A .* vtm + ABCDtm.B .* itmup;
        itmtop = ABCDtm.C .* vtm + ABCDtm.D .* itmup;

        % Green's function
        [Gej] = SpectralGF.ej(z0, k0, kx0, ky0, vtmtop, vtetop, itmtop, itetop);
        Gxx = -(vtm.*cos(ph).^2+vte.*sin(ph).^2);
        Gyx = (vte-vtm).*sin(ph).*cos(ph);
        Gzx = (z0.*vtm)./ztm.*sin(th).*cos(ph);
        Gzx2 = z0.*itmup.*sin(th).*cos(ph);
        if(0)
            Gej.xx - Gxx
            Gej.yx - Gyx
            Gej.zx - Gzx
            Gej.zx - Gzx2
        end
        
        % Fields in xyz
        Ex = Gej.xx;
        Ey = Gej.yx;
        Ez = Gej.zx;
%         Ex = Gxx;
%         Ey = Gyx;
%         Ez = Gzx;

        % Fields in theta-phi
        Er  = sin(th) .* cos(ph) .* Ex + sin(th) .* sin(ph) .* Ey + cos(th) .* Ez;
        Eth = cos(th) .* cos(ph) .* Ex + cos(th) .* sin(ph) .* Ey - sin(th) .* Ez;
        Eph =           -sin(ph) .* Ex +            cos(ph) .* Ey;

        EL3co = cos(ph) .* Eth - sin(ph) .* Eph;
        EL3x  = sin(ph) .* Eth + cos(ph) .* Eph;
        
%         EL3co = -vtm./(2*cos(th)) - vte/2;
%         EL3x  = -vtm./(2*cos(th)) + vte/2;
        
%         EL3co = -z0./4.*exp(-1j.*kz0.*z0).*(1+sec(th));
%         EL3x  = -z0./4.*exp(-1j.*kz0.*z0).*(1-sec(th));
        
        EL3cos(ith, iph, :) = EL3co;
        EL3xs(ith, iph, :) = EL3x;
    end
end

[hFig, hAx] = figureex;
    plot(hAx, ths * 180/pi, 20*log10(abs(EL3cos(:, 1, 1))));
    ylabel('EL3co [dB]');
[hFig, hAx] = figureex;
    plot(hAx, ths * 180/pi, 20*log10(abs(EL3xs(:, 1, 1))));
    ylabel('EL3x [dB]');
[hFig, hAx] = figureex;
    plot(hAx, ths * 180/pi, 20*log10(abs(EL3xs(:, 1, 1) ./ EL3cos(:, 1, 1))), 'LineWidth', 1);
    hAx.LineWidth = 1;
    xlim([0 90]);
    ylim([-40 0]);
    xlabel('Theta [\circ]');
    ylabel('Xpol [dB]');
    alignplot(hFig, 8, 4, hFig.Number, [], 2);
% [hFig, hAx] = figureex;
%     plot(hAx, fs ./ 1e9, 20*log10(abs(EL3x)));
%     ylabel('EL3x [dB]')
% [hFig, hAx] = figureex;
%     plot(hAx, fs ./ 1e9, 20*log10(abs(Ex)));
%     plot(hAx, fs ./ 1e9, 20*log10(abs(Ey)));
%     legend(hAx, {'Ex', 'Ey'});
% [hFig, hAx] = figureex;
%     plot(hAx, fs ./ 1e9, 20*log10(abs(Eth)));
%     plot(hAx, fs ./ 1e9, 20*log10(abs(Eph)));
%     legend(hAx, {'Eth', 'Eph'});










