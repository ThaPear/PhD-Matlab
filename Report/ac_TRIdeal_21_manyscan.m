SetupPath;
clear;
close all;

c0 = Constants.c0;
z0 = Constants.z0;
f0 = 31e9;          % 13.75 to 31 GHz.

z1 = 80;
z2 = z0;
zfeed = 80;

dx = 4.35e-3;
dy = dx;

p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 21e9;
% f0design = 29e9;

Zs = Chebyshev.GetImpedances(gamma, z1, z2, N);

tlines = {};
for(i = 1:N)
    er = (z0 / Zs(i))^2;
    L = (c0/f0match)/4 / sqrt(er);
    tlines = [tlines, {Line(er, L)}];
end

slab = TLine(tlines);

%% Perform simulation.
fs = 1e9 * (0.1:0.1:40);
ths = [eps 10 20 30 40 50 60] * pi/180;
phs = [0   0  0  0  0  0  0] * pi/180;
% ths = [eps 60] * pi/180;
% phs = [0   0] * pi/180;

clrmap = lines(7);
for(iangle = 1:length(ths))
    th = ths(iangle);
    ph = phs(iangle);

    % Calculate propagation constants.
    [k0, kx0, ky0, kz0] = k(fs, 1, th, ph);
    [kd, ~, ~, kzd] = k(fs, (z0/z1)^2, th, ph);
    kr = sqrt(kx0.^2 + ky0.^2);
    % Calculate impedances.
    [z0, z0te, z0tm] = z(1, k0, kz0);

    % Calculate ABCD matrix - TE only since only H-plane is important.
    isTE = 1;
    ABCDte = slab.GetABCD(isTE, fs, k0, kr);
    Ste = ABCD2S(ABCDte, zfeed, z0te);
    Zte = ABCD2Z(ABCDte);
    Zinte = Zte.z11 - Zte.z12.*Zte.z21 ./ (Zte.z22+z0te);
    Zinte(fs == f0match) = (Zinte(circshift(fs == f0match, -1)) + Zinte(circshift(fs == f0match, 1)))/2;
    
    [figSte, axSte] = figureex(1);
        axSte.ColorOrder = clrmap;
        alignplot(figSte, 6, 4, figSte.Number, [], 1);
        if(length(axSte.Children) < 2)
            patch(axSte, [28 31 31 28], [-30 -30 0 0], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(axSte, [13.75 14.5 14.5 13.75], [-30 -30 0 0], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
    [figStm, axStm] = figureex(2);
        axStm.ColorOrder = clrmap;
        alignplot(figStm, 6, 4, figStm.Number, [], 1);
        if(length(axStm.Children) < 2)
            patch(axStm, [28 31 31 28], [-30 -30 0 0], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(axStm, [13.75 14.5 14.5 13.75], [-30 -30 0 0], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
    [figZte, axZte] = figureex(3);
        alignplot(figZte, 6, 4, figZte.Number, [], 1);
        axZte.ColorOrder = reshape(repmat(clrmap, 1, 2).', [], 14).';
        if(length(axZte.Children) < 2)
            patch(axZte, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(axZte, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            plot(axZte, [min(fs), max(fs)]./1e9, [zfeed, zfeed], 'k--');
        end
    [figZtm, axZtm] = figureex(4);
        alignplot(figZtm, 6, 4, figZtm.Number, [], 1);
        axZtm.ColorOrder = reshape(repmat(clrmap, 1, 2).', [], 14).';
        if(length(axZtm.Children) < 2)
            patch(axZtm, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(axZtm, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            plot(axZtm, [min(fs), max(fs)]./1e9, [zfeed, zfeed], 'k--');
        end
    
    plot(axSte, fs/1e9, 20*log10(abs(Ste.s11)));
    addlegendentry(axSte, sprintf('\\theta=%2.0f\\circ', th*180/pi));

    plot(axZte, fs/1e9, real(Zinte));
    addlegendentry(axZte, sprintf('\\theta=%2.0f\\circ', th*180/pi));
    plot(axZte, fs/1e9, imag(Zinte), '--');
    
    isTE = 0;
    ABCDtm = slab.GetABCD(isTE, fs, k0, kr);
    Stm = ABCD2S(ABCDtm, zfeed, z0tm);
    Ztm = ABCD2Z(ABCDtm);
    Zintm = Ztm.z11 - Ztm.z12.*Ztm.z21 ./ (Ztm.z22+z0tm);
    Zintm(fs == f0match) = (Zintm(circshift(fs == f0match, -1)) + Zintm(circshift(fs == f0match, 1)))/2;

    plot(axStm, fs/1e9, 20*log10(abs(Stm.s11)));
    addlegendentry(axStm, sprintf('\\theta=%2.0f\\circ', th*180/pi));

    plot(axZtm, fs/1e9, real(Zintm));
    addlegendentry(axZtm, sprintf('\\theta=%2.0f\\circ', th*180/pi));
    plot(axZtm, fs/1e9, imag(Zintm), '--');
    
    
end
%     xlim(axS, [12 32]);
    ylim(axSte, [-30 0]);
    xlim(axZte, [12 32]);
    ylim(axZte, [-100 200]);
    xlabel(axSte, 'Frequency [GHz]');
    ylabel(axSte, '|\Gamma| [dB]');
    xlabel(axZte, 'Frequency [GHz]');
    ylabel(axZte, 'Input Impedance [\Omega]');
    ylim(axStm, [-30 0]);
    xlim(axZtm, [12 32]);
    ylim(axZtm, [-100 200]);
    xlabel(axStm, 'Frequency [GHz]');
    ylabel(axStm, '|\Gamma| [dB]');
    xlabel(axZtm, 'Frequency [GHz]');
    ylabel(axZtm, 'Input Impedance [\Omega]');

    legend(axZte, 'Orientation', 'horizontal', 'NumColumns', 4);
    legend(axZtm, 'Orientation', 'horizontal', 'NumColumns', 4);
    drawnow;
    
    movelegend(axSte, 'sw');
    movelegend(axStm, 'sw');
    movelegend(axZte, 's');
    movelegend(axZtm, 's');




















