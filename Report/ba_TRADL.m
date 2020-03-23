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
f0design = 29e9;

linewidth = 1;
axlinewidth = 1;

for(iopt = 1:2)
switch(iopt)
    case 1
        Zs = Chebyshev.GetImpedances(gamma, z1, z2, N);

        tlines = {};
        for(i = 1:N)
            er = (z0 / Zs(i))^2;
            L = (c0/f0match)/4 / sqrt(er);
            tlines = [tlines, {Line(er, L)}]; %#ok<AGROW>
        end

        slab = TLine(tlines);
    case 2
        slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);
end


%% Perform simulation.
fs = 1e9 * (0.1:0.1:40);
ths = [eps 20 30 40 50 60] * pi/180;
phs = [0   0  0  0  0  0] * pi/180;
ths = [eps 60] * pi/180;
phs = [0   0] * pi/180;

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
    VSWRte = S2VSWR(Ste.s11);
    
    if(iangle == 1 && iopt == 1)
        Zinte(fs == 22e9) = (Zinte(circshift(fs == 22e9, -1)) + Zinte(circshift(fs == 22e9, 1)))/2;
    end
    
    if(iangle == 1 && iopt == 2)
        [figS, axS] = figureex(1);
        [figVSWR, axVSWR] = figureex(2);
    else
        [figS, axS] = figureex;
        axS.ColorOrder = clrmap;
        if(iangle > 1)
            axS.ColorOrder = reshape(repmat(clrmap, 1, 2).', [], 14).';
        end
        alignplot(figS, 6, 4, figS.Number, [], 1);
        if(length(axS.Children) < 2)
            patch(axS, [28 31 31 28], [-30 -30 0 0], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(axS, [13.75 14.5 14.5 13.75], [-30 -30 0 0], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
        
        [figVSWR, axVSWR] = figureex;
        alignplot(figVSWR, 6, 4, figVSWR.Number, [], 1);
        if(length(axVSWR.Children) < 2)
            patch(axVSWR, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(axVSWR, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
    end
    axS.LineWidth = axlinewidth;
    [figZ, axZ] = figureex;
        alignplot(figZ, 6, 4, figZ.Number, [], 1);
        axZ.ColorOrder = reshape(repmat(clrmap, 1, 2).', [], 14).';
        if(length(axZ.Children) < 2)
            patch(axZ, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(axZ, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            plot(axZ, [min(fs), max(fs)]./1e9, [zfeed, zfeed], 'k--', 'LineWidth', linewidth);
        end
        axZ.LineWidth = axlinewidth;
    
    plot(axS, fs/1e9, 20*log10(abs(Ste.s11)), 'LineWidth', linewidth);
    addlegendentry(axS, 'TE');
    plot(axVSWR, fs/1e9, VSWRte, 'LineWidth', linewidth);
    addlegendentry(axVSWR, 'TE');

    plot(axZ, fs/1e9, real(Zinte), 'LineWidth', linewidth);
    addlegendentry(axZ, 'TE');
    plot(axZ, fs/1e9, imag(Zinte), '--', 'LineWidth', linewidth);
    
    if(iangle > 1)
        isTE = 0;
        ABCDtm = slab.GetABCD(isTE, fs, k0, kr);
        Stm = ABCD2S(ABCDtm, zfeed, z0tm);
        Ztm = ABCD2Z(ABCDtm);
        Zintm = Ztm.z11 - Ztm.z12.*Ztm.z21 ./ (Ztm.z22+z0tm);
        VSWRtm = S2VSWR(Stm.s11);

        plot(axS, fs/1e9, 20*log10(abs(Stm.s11)), '--', 'LineWidth', linewidth);
        addlegendentry(axS, 'TM');
        plot(axVSWR, fs/1e9, VSWRtm, 'LineWidth', linewidth);
        addlegendentry(axVSWR, 'TE');

        plot(axZ, fs/1e9, real(Zintm), 'LineWidth', linewidth);
        addlegendentry(axZ, 'TM');
        plot(axZ, fs/1e9, imag(Zintm), '--', 'LineWidth', linewidth);
    end
%     xlim(axS, [12 32]);
    ylim(axS, [-30 0]);
    xlim(axZ, [12 35]);
    ylim(axZ, [-100 200]);
    xlim(axVSWR, [12 35]);
    ylim(axVSWR, [1 4]);
    xlabel(axS, 'Frequency [GHz]');
    ylabel(axS, '|\Gamma| [dB]');
    xlabel(axZ, 'Frequency [GHz]');
    ylabel(axZ, 'Input Impedance [\Omega]');
    xlabel(axVSWR, 'Frequency [GHz]');
    ylabel(axVSWR, 'VSWR');
    
    linewidth = 1;
    axS.LineWidth = linewidth;
    axZ.LineWidth = linewidth;
    axVSWR.LineWidth = linewidth;
    for(i = 1:length(axS.Children))
        axS.Children(i).LineWidth = linewidth;
    end
    for(i = 1:length(axZ.Children))
        axZ.Children(i).LineWidth = linewidth;
    end
    for(i = 1:length(axVSWR.Children))
        axVSWR.Children(i).LineWidth = linewidth;
    end
end
end

























