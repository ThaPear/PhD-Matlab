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
hgap = -0.1e-3;

shortcore = 1;
patch_area = 1.2;

dslot = 0.3e-3;
wfeed = 0.3e-3;
% wslot = 0.8e-3; wms = 0.13;  hback2 = 1;                            % -15.43dB & -4.63dB
% wslot = 0.9e-3; wms = 0.165;  hback2 = 0.127*4;            % Actual % -17.07dB & -6.29dB & -6.17dB
wslot = 0.9e-3; wms = 0.165;  hback2 = 0.127*4;            % Fake   % -17.48dB & -6.29dB & -6.33dB
%     wfeed = 0.35e-3;                                                % -18.11dB & -5.79dB
%     wfeed = 0.25e-3;                                                % -15.44dB & -6.17dB
%     patch_area = 1.25;                                              % -18.02dB & -6.21dB
%     hback=0.8e-3;                                                   % -13.89dB & -6.25dB
%         patch_area = 1.0;                                           % -12.12dB & -6.05dB
%         patch_area = 0.7;                                           % - 8.39dB & -5.56dB
%         patch_area = 1.4;                                           % - 8.48dB & -2.25dB
%         patch_area = 1.25;                                          % -14.41dB & -6.29dB
%     lms4 = 0.65;                                                    % -16.83dB & -6.25dB & -6.12dB
%     lms4 = 0.67;                                                    % -     dB & -6.22dB & -6.12dB
%     lms4 = 0.7;                                                     % -16.91dB & -6.14dB & -6.30dB
%     lms4 = 0.75;                                                    % -     dB & -6.35dB & -6.01dB
%     lms4 = 0.8;                                                     % -     dB & -6.30dB & -5.97dB
% wslot = 1.0e-3; wms = 0.165;  hback2 = 0.127*4;                     % -17.15dB & -6.16dB
% wslot = 0.85e-3; wms = 0.165;  hback2 = 0.127*4;                    % -16.55dB & -6.25dB


walled = 1;

C = inf;

p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 19e9;
f0design = 29e9;
dualpol = 1;
cavity = 1;
vias = 1;
coax = 1;
slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);

%% Coax feed
dms = 'dms';
lms = 0.5;
wms2 = 'feed_ms_width';
lms2 = 0.4;%'feed_shield_distance-feed_shield_radius-0.1';
core_radius = 0.1125;
core_tophole_radius = 'feed_core_radius';
core_transition_radius = 0.05;
shield_radius = 'feed_core_radius';
shield_distance = 0.53;
shield_startangle = ['90 - feed_shield_totalangle/2 - feed_ms_buried_angle'];
shield_totalangle = 90;
shield_Nvias = 3;

%% Feed below
dms3 = 'hback2/2';
wms3 = 'feed_ms_width';
lms3 = 'feed_shield_distance';
rcore2ms3 = 'feed_core_transition_radius';
rhole3 = 0.25;
dms4 = 0.127;
wms4 = 0.165;
if(~exist('lms4', 'var')); lms4 = 'feed_shield_distance'; end
rcore2ms4 = 'feed_core_transition_radius';
rhole4 = 'feed_groundhole_radius';

%% Patch
patch_angle = 70;
patch_l1 = 1;
patch_capacitance = 0.2e-12;


tlineup = TerminatedTLine(slab, FreeSpace());
tlinedown = ShortedLine(erback, hback, erback);

% The slot.
if(~dualpol)
    slot = Slot(dx, dy, wslot, dslot, walled);
else
    slot = Slot_Dualpol_Bowtie(dx, dy, wslot, dslot);
end

array = InfiniteArray(slot, tlineup, tlinedown);

params = struct('wslot', wslot*1e3, ...
                'dslot', dslot*1e3, ...
                'wfeed', wfeed*1e3, ...
                'th', 0, ...
                'ph', 0, ...
                'coax', coax, ...
                'scor', shortcore, ...
                'Ap', patch_area, ...
                'hb2', hback2, ...
                'wms', wms);
if(hback ~= 0.9e-3)
    params.hb = hback;
end
if(~strcmpi(lms3, 'feed_shield_distance'))
    params.lms3 = lms3;
end
if(~strcmpi(lms4, 'feed_shield_distance'))
    params.lms4 = lms4;
end

%% Perform simulation.
% ths = [eps 20 30 40 50 60] * pi/180;
% phs = [0   0  0  0  0  0] * pi/180;
ths = [eps 60 60] * pi/180;
phs = [0   0  90] * pi/180;
% ths = [eps 60] * pi/180;
% phs = [0   0 ] * pi/180;

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
        gc_below_Build;
        
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
                params.th*180/pi, params.ph*180/pi+90, 20*log10(abs(Gamma22(ind2))), fs(ind2)/1e9, real(ZasC22(ind2)), imag(ZasC22(ind2)));
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
legend(axVSWR, axVSWR.Children([5 4 2 1 3]), {'Broadside', ...
                                              'X-slot H-plane 60\circ', ...
                                              'X-slot E-plane 60\circ', ...
                                              'Y-slot H-plane 60\circ', ...
                                              'Y-slot E-plane 60\circ'});
movelegend(axVSWR, 'n');
%%
hFileOut = fopen([cd, '/PhD-Matlab/ Reports/WP_2100&2400/a_latestCST.m'], 'w');
fprintf(hFileOut, 'winopen(''%s'');', cstfile);
fclose(hFileOut);





















