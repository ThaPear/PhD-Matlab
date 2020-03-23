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
wslot = 1.1e-3;
dslot = 2.5e-3;
walled = 1;

C = inf;%0.2e-12;

p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 19e9;
f0design = 29e9;
dualpol = 0;
cavity = 1;
vias = 1;
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


folder = [cd, '\Report\', mfilename, '\'];
mkdir(folder);

clrmap = lines(7);
for(iangle = 1:length(ths))
    params.th = ths(iangle);
    params.ph = phs(iangle);
    
    filename = struct2string(params);
    
    filepath = [folder, filename, '.cst'];
    if(~exist(filepath, 'file'))
        ea_CST_Build;
        
        project.StoreParameter(CST.Defaults.ThetaName, params.th * 180/pi);
        project.StoreParameter(CST.Defaults.PhiName, params.ph * 180/pi);
        if(params.ph > 1)
            project.StoreParameter('nsamplesperGHz', 1);
        end
    
        project.SaveAs(filepath, 0);
        project.Quit();
        pause(5);
    end
    if(~dualpol)
        touchstonefile = [filepath(1:end-4), '.s1p'];
    else
        touchstonefile = [filepath(1:end-4), '.s2p'];
    end
    if(~exist(touchstonefile, 'file'))
        project = CST.Application.OpenFile(filepath);
        dsproject = CST.Application.ActiveDS();
        
        fdsolver = project.FDSolver();
        fdsolver.Start();
%         dsproject.UpdateResults();
%         dsproject.Save();
%         dsproject.TouchstoneExport('Tasks\SPara1\S-Parameters\S1,1', [filepath(1:end-3), 's1p'], zfeed);
%         dsproject.SelectTreeItem('Tasks\SPara1\S-Parameters\S1,1');
%         project.SelectTreeItem('1D Results\S-Parameters\S1,1');
        
        touchstone = project.Touchstone();
        touchstone.Reset();
        touchstone.Impedance(zfeed);
        touchstone.Renormalize(1);
        touchstone.FileName(filepath(1:end-4));
        touchstone.Write();
    
        project.Save();
        project.Quit();
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
        params.th*180/pi, params.ph*180/pi, 20*log10(abs(Gamma(ind))), fs(ind)/1e9, real(ZasC(ind)), imag(ZasC(ind)));
    
    xlim(axS, [12 inf]);
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
    drawnow
end

























