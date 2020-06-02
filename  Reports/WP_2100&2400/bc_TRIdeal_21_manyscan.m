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

linewidth = 1;
axlinewidth = 1;

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
    
    [figS, axS] = a_figS;
        axS.ColorOrder = clrmap;
        if(iangle > 1)
            axS.ColorOrder = reshape(repmat(clrmap, 1, 2).', [], 14).';
        end
    [figZ, axZ] = a_figZ([], fs, zfeed);
    
    plot(axS, fs/1e9, 20*log10(abs(Ste.s11)));
    addlegendentry(axS, 'TE');

    plot(axZ, fs/1e9, real(Zinte));
    addlegendentry(axZ, 'Real');
    plot(axZ, fs/1e9, imag(Zinte), '--');
    addlegendentry(axZ, 'Imaginary');
    legend(axZ, 'NumColumns', 2);
    
    if(iangle > 1)
        isTE = 0;
        ABCDtm = slab.GetABCD(isTE, fs, k0, kr);
        Stm = ABCD2S(ABCDtm, zfeed, z0tm);
        Ztm = ABCD2Z(ABCDtm);
        Zintm = Ztm.z11 - Ztm.z12.*Ztm.z21 ./ (Ztm.z22+z0tm);

        plot(axS, fs/1e9, 20*log10(abs(Stm.s11)), '--');
        addlegendentry(axS, 'TM');

        plot(axZ, fs/1e9, real(Zintm));
        addlegendentry(axZ, 'TM');
        plot(axZ, fs/1e9, imag(Zintm), '--');
    end
%     xlim(axS, [12 32]);
    ylim(axS, [-30 0]);
    xlim(axZ, [12 32]);
    ylim(axZ, [-100 200]);
end





















