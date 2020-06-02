SetupPath;
clear;
close all;

global test;
test = 0;

c0 = Constants.c0;
z0 = Constants.z0;
f0 = 31e9;          % 13.75 to 31 GHz.

z1 = 80;
z2 = z0;

dx = 4.35e-3; %4.35e-3; %4.354838709677e-3; %c0/f0 * 0.45;
dy = dx;
erback = 2.2;
% hback = 1.9e-3;%1.9496e-3; % = (c0/f0)/sqrt(erback*0.7)/4
hback = (c0/f0)/sqrt(erback*0.7)/4;

p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 21e9;
f0design = 29e9;
slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);

tlineup = TerminatedTLine(slab, FreeSpace());
tlinedown = ShortedLine(erback * 0.7, hback, erback);

% test = 1;

% tlineup = FreeSpace();
% tlinedown = FreeSpace();

%% Simulation
% fs = f0;
ths = pi/180 * [eps, (1:90)];
phs = pi/180 * 45;
fs = (9:1:33) * 1e9;

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

        %% Calculate currents on slot plane.
        % Impedances of up- and down T-lines
        zupte = tlineup.GetInputImpedance(1, fs, k0, kr);       zuptm = tlineup.GetInputImpedance(0, fs, k0, kr);
        zdownte = tlinedown.GetInputImpedance(1, fs, k0, kr);   zdowntm = tlinedown.GetInputImpedance(0, fs, k0, kr);
        % Parallel impedance
        zte = 1 ./ (1 ./ zupte + 1 ./ zdownte);                 ztm = 1 ./ (1 ./ zuptm + 1 ./ zdowntm);
        % Current from Norton equivalent circuit
        ite = sin(ph);                                          itm = cos(ph);
        vte = ite .* zte;                                       vtm = itm .* ztm;
        % Current in upwards T-line
        iteup = vte ./ zupte;                                   itmup = vtm ./ zuptm;

        %% Transform currents through ABCD to top of ADL.
        ABCDte = slab.GetABCD(1, fs, k0, kr);                   ABCDtm = slab.GetABCD(0, fs, k0, kr);

        % Invert ABCD matrix
        invABCDte = ABCDte.invert();                            invABCDtm = ABCDtm.invert();
        % Get voltage and current on top of stratification
        [vtetop, itetop] = invABCDte.mul(vte, iteup);           [vtmtop, itmtop] = invABCDtm.mul(vtm, itmup);

        % Green's function
        [Gej] = SpectralGF.ej(z0, k0, kx0, ky0, vtmtop, vtetop, itmtop, itetop);
        
        % Fields in xyz
        Ex = Gej.xx;
        Ey = Gej.yx;
        Ez = Gej.zx;

        % Fields in theta-phi
        Er  = sin(th) .* cos(ph) .* Ex + sin(th) .* sin(ph) .* Ey + cos(th) .* Ez;
        Eth = cos(th) .* cos(ph) .* Ex + cos(th) .* sin(ph) .* Ey - sin(th) .* Ez;
        Eph =           -sin(ph) .* Ex +            cos(ph) .* Ey;

        EL3co = cos(ph) .* Eth - sin(ph) .* Eph;
        EL3x  = sin(ph) .* Eth + cos(ph) .* Eph;
        
        EL3cos(ith, iph, :) = EL3co;
        EL3xs(ith, iph, :) = EL3x;
    end
end

% [hFig, hAx] = figureex;
%     plot(hAx, ths * 180/pi, 20*log10(abs(EL3cos(:, 1, 1))));
%     xlabel('\theta [\circ]');
%     ylabel('EL3co [dB]');
%     alignplot(hFig, 8, 4, hFig.Number, [], 2);
% [hFig, hAx] = figureex;
%     plot(hAx, ths * 180/pi, 20*log10(abs(EL3xs(:, 1, 1))));
%     xlabel('\theta [\circ]');
%     ylabel('EL3x [dB]');
%     alignplot(hFig, 8, 4, hFig.Number, [], 2);
[hFig, hAx] = figureex;
    plot(hAx, ths * 180/pi, 20*log10(abs(EL3xs(:, 1, 21) ./ EL3cos(:, 1, 21))), 'LineWidth', 1); % Index 21 = 29 GHz.
    hAx.LineWidth = 1;
    xlim([0 90]);
    ylim([-40 0]);
    xlabel('Theta [\circ]');
    ylabel('Xpol [dB]');
    alignplot(hFig, 8, 4, hFig.Number, [], 2);
    
    
% [hFig, hAx] = figureex;
%     plot(hAx, ths * 180/pi, 20*log10(abs(EL3cos(:, 1, 1))));
%     xlabel('\theta [\circ]');
%     ylabel('EL3co [dB]');
%     alignplot(hFig, 8, 4, hFig.Number, [], 2);
% [hFig, hAx] = figureex;
%     plot(hAx, ths * 180/pi, 20*log10(abs(EL3xs(:, 1, 1))));
%     xlabel('\theta [\circ]');
%     ylabel('EL3x [dB]');
%     alignplot(hFig, 8, 4, hFig.Number, [], 2);

[hFig, hAx] = figureex;
    plot(hAx, fs ./ 1e9, squeeze(20*log10(abs(EL3xs(61, 1, :) ./ EL3cos(61, 1, :)))));
    xlim([-inf inf]);
    ylim([-40 0]);
    xlabel('Frequency [GHz]');
    ylabel('Xpol [dB]');
    alignplot(hFig, 8, 4, hFig.Number, [], 2);
    
% [hFig, hAx] = figureex;
%     plot(hAx, fs ./ 1e9, 20*log10(abs(Eth ./ Eph)));
%     ylabel('Eth ./ Eph [dB]');
%     alignplot(hFig, 8, 4, hFig.Number, [], 2);
% [hFig, hAx] = figureex;
%     plot(hAx, fs ./ 1e9, 20*log10(abs(Ex)));
%     plot(hAx, fs ./ 1e9, 20*log10(abs(Ey)));
%     legend(hAx, {'Ex', 'Ey'});
%     alignplot(hFig, 8, 4, hFig.Number, [], 2);
% [hFig, hAx] = figureex;
%     plot(hAx, fs ./ 1e9, 20*log10(abs(Eth)));
%     plot(hAx, fs ./ 1e9, 20*log10(abs(Eph)));
%     legend(hAx, {'Eth', 'Eph'});
%     alignplot(hFig, 8, 4, hFig.Number, [], 2);










