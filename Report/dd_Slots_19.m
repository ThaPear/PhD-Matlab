SetupPath;
clear;
close all;

c0 = Constants.c0;
z0 = Constants.z0;
f0 = 31e9;          % 13.75 to 31 GHz.

z1 = 80;
z2 = z0;
zfeed = 80;

dx = 4.35e-3; %4.35e-3; %4.354838709677e-3; %c0/f0 * 0.45;
dy = dx;
erback = 2.2;
hback = 1.9e-3;%1.9496e-3; % = (c0/f0)/sqrt(erback*0.7)/4
wslot = 1.4e-3;
dslot = 2e-3;
walled = 1;

C = inf;%0.2e-12;

p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 19e9;
f0design = 29e9;
slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);

tlineup = TerminatedTLine(slab, FreeSpace());
tlinedown = ShortedLine(erback, hback, erback);

% The slot.
slot = Slot(dx, dy, wslot, dslot, tlineup, tlinedown, walled);
% slot = Slot_Dualpol_Bowtie(dx, dy, wslot, dslot, tlineup, tlinedown);

%% Perform simulation.
fs = 1e9 * (0.1:0.1:40);
% ths = [eps 20 30 40 50 60] * pi/180;
% phs = [0   0  0  0  0  0] * pi/180;
ths = [eps 60 60] * pi/180;
phs = [0   0  90] * pi/180;

Zcap = 1 ./ (1j .* 2.*pi.*fs .* C);

clrmap = lines(7);
for(iangle = 1:length(ths))
    th = ths(iangle);
    ph = phs(iangle);    
    
    tc = tic;
    Zas = slot.GetInputImpedance(fs, th, ph);
%     dispex('Calculated input impedance in %.3fs.\n', toc(tc));
    ZasC = Zas + Zcap; % With series capacitance.
    
    % Calculate reflection coefficient.
    Gamma = (ZasC - zfeed) ./ (ZasC + zfeed);
    VSWR = (1 + abs(Gamma)) ./ (1 - abs(Gamma));
    
    [figS, axS] = a_figS(6);
    [figZ, axZ] = a_figZ(7, fs, zfeed);
    [figVSWR, axVSWR] = a_figVSWR(8);
    
    plot(axS, fs/1e9, 20*log10(abs(Gamma)));
    plot(axVSWR, fs/1e9, VSWR);

    plot(axZ, fs/1e9, real(ZasC));
    plot(axZ, fs/1e9, imag(ZasC), '--');

    ind = find(20*log10(abs(Gamma)) == max(20*log10(abs(Gamma(fs >= 13.75e9 & fs <= 14.5e9 | fs >= 28e9 & fs <= 31e9)))));
    dispex('Worst at %02.0f,%02.0f is %.2f, f = %.1fGHz, Z = %.2f + %.2fj.\n', ...
        th*180/pi, ph*180/pi, 20*log10(abs(Gamma(ind))), fs(ind)/1e9, real(ZasC(ind)), imag(ZasC(ind)));
    
    xlim(axS, [12 35]);
    ylim(axS, [-30 0]);
    xlim(axZ, [12 32]);
    ylim(axZ, [-100 200]);
    xlim(axVSWR, [12 32]);
    ylim(axVSWR, [1 8]);
end

legend(axVSWR, axVSWR.Children(3:-1:1), {'Broadside', 'H-plane 60\circ', 'E-plane 60\circ'});
movelegend(axVSWR, 'n');






















