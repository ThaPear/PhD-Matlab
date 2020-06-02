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
dualpol = 0;
cavity = 1;
vias = 0;
slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);

tlineup = TerminatedTLine(slab, FreeSpace());
tlinedown = ShortedLine(erback, hback, erback);

% The slot.
if(~dualpol)
    slot = Slot(dx, dy, wslot, dslot, tlineup, tlinedown, walled);
else
    slot = Slot_Dualpol_Bowtie(dx, dy, wslot, dslot, tlineup, tlinedown);
end

params = struct('z1', z1, ...
                'z2', z2, ...
                'zfeed', zfeed, ...
                'dx', dx, ...
                'dy', dy, ...
                'erback', erback, ...
                'hback', hback, ...
                'wslot', wslot, ...
                'dslot', dslot, ...
                'walled', walled, ...
                'p', p, ...
                'gamma', gamma, ...
                'N', N, ...
                'f0match', f0match, ...
                'f0design', f0design, ...
                'dualpol', dualpol, ...
                'th', 0, ...
                'ph', 0, ...
                'cavity', cavity, ...
                'vias', vias);

%% Perform simulation.
% ths = [eps 20 30 40 50 60] * pi/180;
% phs = [0   0  0  0  0  0] * pi/180;
ths = [eps 60 60] * pi/180;
phs = [0   0  90] * pi/180;


folder = ['e:\ Reports\WP_2100&2400\', mfilename, '\'];
[~, ~] = mkdir(folder);

clrmap = lines(7);
for(iangle = 1:length(ths))
    params.th = ths(iangle);
    params.ph = phs(iangle);
    
    filename = struct2string(params);
    
    filepath = [folder, filename, '.cst'];
    if(~exist(filepath, 'file'))
        dispex('Building simulation...'); tc = tic;
        ea_CST_Build;
        
        project.StoreParameter(CST.Defaults.ThetaName, params.th * 180/pi);
        project.StoreParameter(CST.Defaults.PhiName, params.ph * 180/pi);
        if(params.ph > 1)
            project.StoreParameter('nsamplesperGHz', 1);
        end
    
        project.SaveAs(filepath, 0);
        fprintf(' Done. Took %.1fs.\n', toc(tc));
        pause(5);
        project.Quit();
        pause(10);
    end
    if(~dualpol)
        touchstonefile = [filepath(1:end-4), '.s1p'];
    else
        touchstonefile = [filepath(1:end-4), '.s2p'];
    end
    if(~exist(touchstonefile, 'file'))
        dispex('Simulating...'); tc = tic;
        project = CST.Application.OpenFile(filepath);
        dsproject = CST.Application.ActiveDS();
        
        project.Rebuild();
        
        fdsolver = project.FDSolver();
        if(~fdsolver.Start()); disp('Simulation failed.'); return; end
        
        touchstone = project.Touchstone();
        touchstone.Reset();
        touchstone.Impedance(zfeed);
        touchstone.Renormalize(1);
        touchstone.FileName(filepath(1:end-4));
        touchstone.Write();
    
        project.Save();
        project.Quit();
        fprintf(' Done. Took %.1fs.\n', toc(tc));
    end
    
    [parameters, S] = CST.LoadData(touchstonefile);
    fs = parameters.frequencies;
    
    Z = (1+S)./(1-S).*parameters.slot_impedance;
    
    tc = tic;
    Zas = squeeze(Z).';
    Zcap = 1 ./ (1j .* 2.*pi.*fs .* C);
    ZasC = Zas + Zcap; % With series capacitance.
    
    % Calculate reflection coefficient.
    Gamma = (ZasC - zfeed) ./ (ZasC + zfeed);
    VSWR = (1 + abs(Gamma)) ./ (1 - abs(Gamma));
    
    [figS, axS] = a_figS(6);
    [figZ, axZ] = a_figZ(7, fs, zfeed);
    [figVSWR, axVSWR] = a_figVSWR(8);
    
    plot(axS, fs/1e9, 20*log10(abs(Gamma)));

    plot(axZ, fs/1e9, real(ZasC));
    plot(axZ, fs/1e9, imag(ZasC), '--');
    
    plot(axVSWR, fs/1e9, VSWR);

    ind = find(20*log10(abs(Gamma)) == max(20*log10(abs(Gamma(fs >= 13.75e9 & fs <= 14.5e9 | fs >= 28e9 & fs <= 31e9)))), 1);
    dispex('Worst at %02.0f,%02.0f is %.2f, f = %.1fGHz, Z = %.2f + %.2fj.\n', ...
        params.th*180/pi, params.ph*180/pi, 20*log10(abs(Gamma(ind))), fs(ind)/1e9, real(ZasC(ind)), imag(ZasC(ind)));
    
    xlim(axS, [12 inf]);
    ylim(axS, [-30 0]);
    xlim(axZ, [12 32]);
    ylim(axZ, [-100 200]);
    xlim(axVSWR, [12 32]);
    ylim(axVSWR, [1 8]);
end

legend(axVSWR, axVSWR.Children(3:-1:1), {'Broadside', 'H-plane 60\circ', 'E-plane 60\circ'});
movelegend(axVSWR, 'n');






















