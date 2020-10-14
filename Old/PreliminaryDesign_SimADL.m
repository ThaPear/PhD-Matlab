%% Simulation setup.
fs = 1e9 * (0.1:0.1:40);
ths = [eps 20 30 40 50 60] * pi/180;
phs = [0   0  0  0  0  0] * pi/180;
ths = [eps 60] * pi/180;
phs = [0   0] * pi/180;

%% Perform simulation.
styles = {'-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'};
for(iangle = 1:length(ths))
    th = ths(iangle);
    ph = phs(iangle);

    % Calculate propagation constants.
    [k0, kx0, ky0, kz0] = k(fs, 1, th, ph);
%     [kd, ~, ~, kzd] = k(fs, (z0/z1)^2, th, ph);
    kr = sqrt(kx0.^2 + ky0.^2);
    % Calculate impedances.
    [z0, z0te, z0tm] = z(1, k0, kz0);
%     [zd, zdte, zdtm] = z((z0/z1)^2, kd, kzd);

    % Calculate ABCD matrix - TE only since only H-plane is important.
    isTE = 1;
    ABCDte = slab.GetABCD(isTE, fs, k0, kr);
    Ste = ABCD2S(ABCDte, zfeed, z0te);
    Zte = ABCD2Z(ABCDte);
    Zinte = Zte.z11 - Zte.z12.*Zte.z21 ./ (Zte.z22+z0te);

    % If the axes are empty, put in the frequency range indicator.
    % Plot results.
    figS = figureex(iangle); axS = figS.CurrentAxes;
        alignplot(figS, 4, 3, figS.Number, [], 2);
        if(length(axS.Children) < 2)
%             patch(ax, [29 31 31 29], [-10+3*(iangle>1) -10+3*(iangle>1) 0 0], [0 0 0], ...
%                 'FaceAlpha', 0.1, 'EdgeColor', 'none');
%             patch(ax, [13.75 14 14 13.75], [-10+3*(iangle>1) -10+3*(iangle>1) 0 0], [0 0 0], ...
%                 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(axS, [28 31 31 28], [-30 -30 0 0], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
            patch(axS, [13.75 14.5 14.5 13.75], [-30 -30 0 0], [0 0 0], ...
                'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
        ylim([-30 0]);
        xlim([-inf inf]);
        xlim([12 32]);
        xlabel('Frequency [GHz]');
        ylabel('|\Gamma| [dB]');

        plot(axS, fs/1e9, 20*log10(abs(Ste.s11)), styles{iangle});
        addlegendentry(axS, 'TE');
        
    figZ = figureex(iangle+length(ths)); axZ = figZ.CurrentAxes;
        alignplot(figZ, 4, 3, figZ.Number, [], 2);
        axZ.ColorOrder = reshape(repmat(lines(7), 1, 2).', [], 14).';
        plot(axZ, fs/1e9, real(Zinte));
        addlegendentry(axZ, 'TE');
        plot(axZ, fs/1e9, imag(Zinte), '--');
        
    if(iangle > 1)
        isTE = 0;
        ABCDtm = slab.GetABCD(isTE, fs, k0, kr);
        Stm = ABCD2S(ABCDtm, zfeed, z0tm);
        Ztm = ABCD2Z(ABCDtm);
        Zintm = Ztm.z11 - Ztm.z12.*Ztm.z21 ./ (Ztm.z22+z0tm);

        plot(axS, fs/1e9, 20*log10(abs(Stm.s11)), styles{iangle});
        addlegendentry(axS, 'TM');

        plot(axZ, fs/1e9, real(Zintm));
        addlegendentry(axZ, 'TM');
        plot(axZ, fs/1e9, imag(Zintm), '--');
    end
end