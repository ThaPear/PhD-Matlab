SetupPath;
clear;
close all;

c0 = Constants.c0;
z0 = Constants.z0;
f0 = 31e9;          % 13.75 to 31 GHz.

z1 = 90;
z2 = z0;
zfeed = 80;

dx = 4.35e-3; %4.35e-3; %4.354838709677e-3; %c0/f0 * 0.45;
dy = dx;
erback = 2.2;
hback = 0.9e-3;
% wslot = 0.75e-3; dslot = 0.25e-3; wfeed = 0.4e-3;                               % - 9.51 & -4.94dB
% wslot = 0.65e-3; dslot = 1.0e-3; wfeed = wslot;                                 % -10.62 & -5.11dB
% wslot = 0.5e-3; dslot = 1.0e-3; wfeed = wslot;                                  % -12.85 & -5.25dB
% wslot = 0.5e-3; dslot = 0.7e-3; wfeed = wslot;                                  % -13.76 & -5.70dB
% wslot = 0.6e-3; dslot = 0.7e-3; wfeed = 0.6e-3;                                 % -11.31 & -5.40dB
% wslot = 0.4e-3; dslot = 0.7e-3; wfeed = 0.6e-3;                                 % -15.80 & -6.27dB
% wslot = 0.3e-3; dslot = 0.7e-3; wfeed = 0.6e-3;                                 % -17.04 & -6.23dB
% wslot = 0.35e-3; dslot = 0.7e-3; wfeed = 0.6e-3;                                % -16.54 & -6.37dB
wslot = 0.35e-3; dslot = 0.6e-3; wfeed = 0.6e-3;                                % -15.76 & -6.74dB
% wslot = 0.3e-3; dslot = 0.6e-3; wfeed = 0.6e-3;                                 % -16.34 & -6.63dB
% wslot = 0.35e-3; dslot = 0.5e-3; wfeed = 0.6e-3;                                % -15.31 & -6.06dB
% wslot = 0.35e-3; dslot = 0.6e-3; wfeed = 0.55e-3;                               % -16.80 & -6.50dB

walled = 1;

C = 0.2e-12;
C = inf;

p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 19e9;
f0design = 29e9;
dualpol = 1;
cavity = 1;
vias = 1;
slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);

tlineup = TerminatedTLine(slab, FreeSpace());
tlinedown = ShortedLine(erback, hback, erback);

% The slot.
if(~dualpol)
    slot = Slot(dx, dy, wslot, dslot, walled);
else
    slot = Slot_Dualpol_Bowtie(dx, dy, wslot, dslot);
end

array = InfiniteArray(slot, tlineup, tlinedown);

params = struct('z1', z1, ...
                'z2', z2, ...
                'zfeed', zfeed, ...
                'hback', hback, ...
                'wslot', wslot, ...
                'dslot', dslot, ...
                'wfeed', wfeed, ...
                'dualpol', dualpol, ...
                'th', 0, ...
                'ph', 0, ...
                'cavity', cavity, ...
                'vias', vias);

%% Perform simulation.
% ths = [eps 20 30 40 50 60] * pi/180;
% phs = [0   0  0  0  0  0] * pi/180;
% ths = [eps 60 60] * pi/180;
% phs = [0   0  90] * pi/180;
ths = [eps 60] * pi/180;
phs = [0   0 ] * pi/180;

folder = sprintf('%s/ Reports/WP_2100&2400/%s/', resultdir_cst, mfilename);
touchstonefolder = sprintf('%s/ Reports/WP_2100&2400/%s/', resultdir_matlab, mfilename);
[~, ~] = mkdir(folder);
[~, ~] = mkdir(touchstonefolder);

clrmap = lines(7);
for(iangle = 1:length(ths))
    params.th = ths(iangle);
    params.ph = phs(iangle);
    
    filename = struct2string(params);
    
    cstfile = [folder, filename, '.cst'];
    if(~dualpol)
        touchstonefile = [touchstonefolder, filename, '.s1p'];
    else
        touchstonefile = [touchstonefolder, filename, '.s2p'];
    end
    
    if(~exist(cstfile, 'file'))
        dispex('Building simulation...'); tc = tic;
        ea_CST_Build;
        
        project.StoreParameter(CST.Defaults.ThetaName, params.th * 180/pi);
        project.StoreParameter(CST.Defaults.PhiName, params.ph * 180/pi);
        if(params.ph > 1)
            project.StoreParameter('nsamplesperGHz', 1);
        end
        
        project.SaveAs(cstfile, 0);
        fprintf(' Done. Took %.1fs.\n', toc(tc));
        pause(5);
        project.Quit();
        pause(10);
    end
    if(~exist(touchstonefile, 'file'))
        dispex('Simulating...'); tc = tic;
        project = CST.Application.OpenFile(cstfile);
        dsproject = CST.Application.ActiveDS();
        
        project.Rebuild();
        
        fdsolver = project.FDSolver();
        if(~fdsolver.Start()); disp('Simulation failed.'); return; end
        
        touchstone = project.TOUCHSTONE();
        touchstone.Reset();
        touchstone.Impedance(zfeed);
        touchstone.Renormalize(1);
        touchstone.FileName(touchstonefile(1:end-4));
        touchstone.Write();
    
        project.Save();
        project.Quit();
        fprintf(' Done. Took %.1fs.\n', toc(tc));
    end
    
    [parameters, S] = CST.LoadData(touchstonefile);
    fs = parameters.frequencies;
    
    Z = (1+S)./(1-S).*parameters.slot_impedance;
    
    tc = tic;
    Zas11 = squeeze(Z(1,1,:)).';
    Zas22 = squeeze(Z(2,2,:)).';
    Zcap = 1 ./ (1j .* 2.*pi.*fs .* C);
    ZasC11 = Zas11 + Zcap; % With series capacitance.
    ZasC22 = Zas22 + Zcap; % With series capacitance.
    
    % Calculate reflection coefficient.
    Gamma11 = (ZasC11 - zfeed) ./ (ZasC11 + zfeed);
    Gamma22 = (ZasC22 - zfeed) ./ (ZasC22 + zfeed);
    VSWR11 = (1 + abs(Gamma11)) ./ (1 - abs(Gamma11));
    VSWR22 = (1 + abs(Gamma22)) ./ (1 - abs(Gamma22));
    
    [figS, axS] = a_figS(6);
    [figZ, axZ] = a_figZ(7, fs, zfeed);
    [figVSWR, axVSWR] = a_figVSWR(8);
    
    plot(axS, fs/1e9, 20*log10(abs(Gamma11)));

    plot(axZ, fs/1e9, real(ZasC11));
    plot(axZ, fs/1e9, imag(ZasC11), '--');
    
    plot(axVSWR, fs/1e9, VSWR11);
    
    if(iangle > 1)
        plot(axS, fs/1e9, 20*log10(abs(Gamma22)));
        plot(axZ, fs/1e9, real(ZasC22));
        plot(axZ, fs/1e9, imag(ZasC22), '--');
        plot(axVSWR, fs/1e9, VSWR22);
    end

    if(iangle > 1)
        ind = find(20*log10(abs(Gamma11)) == max(20*log10(abs(Gamma11(fs >= 13.75e9 & fs <= 14.5e9 | fs >= 28e9 & fs <= 31e9)))), 1);
        ind2 = find(20*log10(abs(Gamma22)) == max(20*log10(abs(Gamma22(fs >= 13.75e9 & fs <= 14.5e9 | fs >= 28e9 & fs <= 31e9)))), 1);
        if(20*log10(abs(Gamma11(ind))) < 20*log10(abs(Gamma22(ind2))))
            dispex('Worst at %02.0f,%02.0f is %.2f, f = %.1fGHz, Z = %.2f + %.2fj.\n', ...
                params.th*180/pi, (params.ph*180/pi+90), 20*log10(abs(Gamma22(ind2))), fs(ind2)/1e9, real(ZasC22(ind2)), imag(ZasC22(ind2)));
        else
            dispex('Worst at %02.0f,%02.0f is %.2f, f = %.1fGHz, Z = %.2f + %.2fj.\n', ...
                params.th*180/pi, params.ph*180/pi, 20*log10(abs(Gamma11(ind))), fs(ind)/1e9, real(ZasC11(ind)), imag(ZasC11(ind)));
        end
    else
        ind = find(20*log10(abs(Gamma11)) == max(20*log10(abs(Gamma11(fs >= 13.75e9 & fs <= 14.5e9 | fs >= 28e9 & fs <= 31e9)))));
        dispex('Worst at %02.0f,%02.0f is %.2f, f = %.1fGHz, Z = %.2f + %.2fj.\n', ...
            params.th*180/pi, params.ph*180/pi, 20*log10(abs(Gamma11(ind))), fs(ind)/1e9, real(ZasC11(ind)), imag(ZasC11(ind)));
    end
    
    xlim(axS, [12 inf]);
    ylim(axS, [-30 0]);
    xlim(axZ, [12 32]);
    ylim(axZ, [-100 200]);
    xlim(axVSWR, [12 32]);
    ylim(axVSWR, [1 8]);
end

legend(axVSWR, axVSWR.Children(3:-1:1), {'Broadside', 'H-plane 60\circ', 'E-plane 60\circ'});
movelegend(axVSWR, 'n');

























