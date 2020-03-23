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
walled = 0;

C = inf;%0.2e-12;

p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 21e9;
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
ths = [eps] * pi/180;
phs = [0  ] * pi/180;

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
    
    [figS, axS] = figureex(7);
    axS.ColorOrder = clrmap;
%     if(iangle > 1)
%         axS.ColorOrder = reshape(repmat(clrmap, 1, 2).', [], 14).';
%     end
    alignplot(figS, 8, 4, figS.Number, [], 1);
    if(length(axS.Children) < 2)
        patch(axS, [28 31 31 28], [-30 -30 0 0], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(axS, [13.75 14.5 14.5 13.75], [-30 -30 0 0], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    [figZ, axZ] = figureex(8);
        alignplot(figZ, 8, 4, figZ.Number, [], 1);
        axZ.ColorOrder = reshape(repmat(clrmap, 1, 2).', [], 14).';
        if(length(axZ.Children) < 2)
            patch(axZ, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(axZ, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            plot(axZ, [min(fs), max(fs)]./1e9, [zfeed, zfeed], 'k--');
        end
    
    plot(axS, fs/1e9, 20*log10(abs(Gamma)));

    plot(axZ, fs/1e9, real(ZasC));
    plot(axZ, fs/1e9, imag(ZasC), '--');

    ind = find(20*log10(abs(Gamma)) == max(20*log10(abs(Gamma(fs >= 13.75e9 & fs <= 14.5e9 | fs >= 28e9 & fs <= 31e9)))));
    dispex('Worst at %02.0f,%02.0f is %.2f, f = %.1fGHz, Z = %.2f + %.2fj.\n', ...
        th*180/pi, ph*180/pi, 20*log10(abs(Gamma(ind))), fs(ind)/1e9, real(ZasC(ind)), imag(ZasC(ind)));
    
%     xlim(axS, [12 35]);
    ylim(axS, [-30 0]);
    xlim(axZ, [12 32]);
    ylim(axZ, [-100 200]);
    xlabel(axS, 'Frequency [GHz]');
    ylabel(axS, '|\Gamma| [dB]');
    xlabel(axZ, 'Frequency [GHz]');
    ylabel(axZ, 'Input Impedance [\Omega]');
    
    linewidth = 1;
    axS.LineWidth = linewidth;
    axZ.LineWidth = linewidth;
    for(i = 1:length(axS.Children))
        axS.Children(i).LineWidth = linewidth;
    end
    for(i = 1:length(axZ.Children))
        axZ.Children(i).LineWidth = linewidth;
    end
end

























