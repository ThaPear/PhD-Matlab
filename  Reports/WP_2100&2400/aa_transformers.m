SetupPath;
clear;
close all;

for(iopt = 1:4)

    c0 = Constants.c0;
    z0 = Constants.z0;
    f0 = 31e9;          % 13.75 to 31 GHz.

    z1 = 80;
    z2 = z0;
    zfeed = z1;

    dx = 4.35e-3;
    dy = dx;

    p = dx / 2;
    gamma = 10^(-10/20);%0.2;
    N = 2;
    f0match = 22.35e9;
    % f0design = 29e9;
    L = N * 0.25 * (c0 / f0match);
    
    switch(iopt)
        case 1
            leg = 'Quarter-wave';
            Zs = sqrt(z1 * z2);
%         case 2
%             leg = 'Multi-Qwave';
%             for(n = 1:N)
%                 Zs(n) = z1^((N+1-n)/(N+1)) * (z2^(n/(N+1)));
%             end
        case 2
            leg = 'Exponential';
            
            L =  N * 0.25 * (c0 / f0match);
            
            a = 1 / L * log(z2 / z1);
            
            Z = @(z) z1 * exp(a * z);
            len = L/N;
            for(n = 1:N)
                Zs(n) = Z(n*len - len/2);
            end
            
            Gamma = @(f) log(z2/z1)/2 .* exp(-1j .* 2.*pi.*f/c0 .* L) .* sin(2.*pi.*f/c0 .* L) ./ (2.*pi.*f/c0 .* L);
        case 3
            leg = 'Chebyshev';
            Zs = Chebyshev.GetImpedances(gamma, z1, z2, N);
        case 4
            leg = 'Binomial';
            Zs = Binomial.GetImpedances(z1, z2, N);
    end

    tlines = {};
    for(i = 1:length(Zs))
        er = (z0 / Zs(i))^2;
        l = L/N / sqrt(er);
        tlines = [tlines, {Line(l, er)}];
    end

    slab = TLine(tlines);

    %% Perform simulation.
    fs = 1e9 * (0.1:0.1:40);
    ths = [eps] * pi/180;
    phs = [0  ] * pi/180;

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

        [figS, axS] = a_figS(1);
        [figZ, axZ] = a_figZ(2, nan, nan);

        if(iopt == 2)
            plot(axS, fs/1e9, 20*log10(abs(Gamma(fs))));
        else
            plot(axS, fs/1e9, 20*log10(abs(Ste.s11)));
        end
    %     addlegendentry(axS, 'TE');

        plot(axZ, fs/1e9, real(Zinte));

        xlim(axS, [12 32]);
        ylim(axS, [-30 0]);
        xlim(axZ, [12 32]);
        ylim(axZ, [-100 200]);
    end
    
    addlegendentry(axS, leg);
end

alignplot(figS, 4, 4, figS.Number, [], 2);
figureex(1);
legend(axS, 'NumColumns', 2);
movelegend(axS, 'n');



















