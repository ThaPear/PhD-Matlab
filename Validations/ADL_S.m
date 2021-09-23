% close all;
clear;
SetupPath;

f0 = 5e9;
c0 = Constants.c0;
l0 = c0 / f0;

z0 = Constants.z0;
zfeed = 180;

%% ADL with given parameters.
%{-
p = 0.0785*l0;
ws = repmat(0.01*l0, 1, 5);
ds = 0.012*l0 * ones(1,4);
% ds = [0, ds, 0];
ds = [ds(1)/2, ds, ds(end)/2];
% ss = repmat(p/2, 1, length(ws));
ss = [0 0.1 0.3 0.4 0] * p;
slab = ADS(p, ds, ss, ws, 1);
filename = sprintf('adltest');
%}
%% ADL with given epsilon.
%{
erdes = 10;
p = 0.2*l0;
L = l0/2/sqrt(erdes);

if(~exist('slab', 'var'))
    slab = DesignADS(f0, p, L, erdes);
end


filename = sprintf('erdes_%g', erdes);
%}

% DrawADS(slab, f0);
% drawnow;

fs = 1e9 * (0.1:0.1:20);
th = 60 * pi/180;
ph = 0 * pi/180;

% Calculate propagation constants.
[k0, kx0, ky0, kz0] = k(fs, 1, th, ph);
kr = sqrt(kx0.^2 + ky0.^2);
% Calculate impedances.
[z0, z0te, z0tm] = z(1, k0, kz0);
%     [zd, zdte, zdtm] = z((z0/z1)^2, kd, kzd);

%% Simulate analytically.
ABCDte = slab.GetABCD(1, fs, k0, kr);                   ABCDtm = slab.GetABCD(0, fs, k0, kr);
% Ste = ABCD2S(ABCDte, zfeed, z0te);                      Stm = ABCD2S(ABCDtm, zfeed, z0tm);
Ste = ABCD2S(ABCDte, z0te, z0te);                      Stm = ABCD2S(ABCDtm, z0tm, z0tm);
Zte = ABCD2Z(ABCDte);                                   Ztm = ABCD2Z(ABCDtm);
Zinte = Zte.z11 - Zte.z12.*Zte.z21 ./ (Zte.z22+z0te);   Zintm = Ztm.z11 - Ztm.z12.*Ztm.z21 ./ (Ztm.z22+z0tm);
% Plot results.
% [figS, axS] = figureex;
%     repeatcolormap(axS, 2);
%     plot(axS, fs/1e9, real(Ste.s11));
%     addlegendentry(axS, 'TE');
%     plot(axS, fs/1e9, imag(Ste.s11), '--');
%     plot(axS, fs/1e9, real(Stm.s11));
%     addlegendentry(axS, 'TM');
%     plot(axS, fs/1e9, imag(Stm.s11), '--');
%     
%     xlabel(axS, 'Frequency [GHz]');
%     ylabel(axS, '\Gamma');
%     ylim(axS, [-1 1]);
%     figS.Name = 'Analytical';
%     drawnow;
% Plot results.
[figAbsS, axAbsS] = figureex;
    repeatcolormap(axAbsS, 2);
    plot(axAbsS, fs/1e9, abs(Ste.s11));
    addlegendentry(axAbsS, 'TE');
    plot(axAbsS, fs/1e9, abs(Ste.s12), '--');
    plot(axAbsS, fs/1e9, abs(Stm.s11));
    addlegendentry(axAbsS, 'TM');
    plot(axAbsS, fs/1e9, abs(Stm.s12), '--');
    
    xlabel(axAbsS, 'Frequency [GHz]');
    ylabel(axAbsS, '|\Gamma|');
    ylim(axAbsS, [0 1]);
    
    figAbsS.Name = 'Analytical';
    drawnow;
    

%% Determine path to place the CST files.
touchstonepath = sprintf('%s/Validations/%s/', resultdir_matlab, mfilename);
path = sprintf('%s/Validations/%s/', resultdir_matlab, mfilename);

cstfilepath = [path, filename];
touchstonefilepath = [touchstonepath, filename];

%% Validate using CST simulation.
if(~exist([touchstonefilepath, '.s4p'], 'file'))
    tc = tic;
    project = CST.InitializeUnitCellProject();
    project.StoreParameter('aa_theta', round(th * 180/pi));
    project.StoreParameter('aa_phi', round(ph * 180/pi));
    project.StoreParameter('fmin', 1);
    project.StoreParameter('fmax', 20);
    project.StoreParameter('fmesh', 20);
    project.StoreParameter('dx', 2*p*1e3);
    project.StoreParameter('dy', 'dx');
    project.StoreParameter('nsamplesperGHZ', 1);

    slab.BuildCST(project);
    project.Rebuild();
    
    % Run simulation.
    fdsolver = project.FDSolver();
    if(~fdsolver.Start())
        warning('CST Simulation failed.');
        return;
    end

    % Export results.
    touchstone = project.TOUCHSTONE();
    touchstone.Reset();
    touchstone.Impedance(z0);
    touchstone.Renormalize(1);
    touchstone.FileName(touchstonefilepath);
    touchstone.Write();
    dispex('CST took %s.\n', fancyduration(toc(tc)));
end

% Import results.
[parameters, S] = CST.LoadData([touchstonefilepath, '.s4p']);
fs = parameters.frequencies;
Ste = [];                                               Stm = [];
Ste.s11 = squeeze(S(1,1,:));                            Stm.s11 = squeeze(S(2,2,:));
Ste.s12 = squeeze(S(1,3,:));                            Stm.s12 = squeeze(S(2,4,:));
Ste.s21 = squeeze(S(3,1,:));                            Stm.s21 = squeeze(S(4,2,:));
Ste.s22 = squeeze(S(3,3,:));                            Stm.s22 = squeeze(S(4,4,:));

% Renormalize
ABCDte = S2ABCD(Ste, z0, z0);                           ABCDtm = S2ABCD(Stm, z0, z0);
% Ste = ABCD2S(ABCDte, zfeed, z0te(1));                   Stm = ABCD2S(ABCDtm, zfeed, z0tm(1));
Ste = ABCD2S(ABCDte, z0te(1), z0te(1));                      Stm = ABCD2S(ABCDtm, z0tm(1), z0tm(1));

% % [figS, axS] = figureex(2);
%     plot(axS, fs/1e9, real(Ste.s11));
%     addlegendentry(axS, 'CSTTE');
%     plot(axS, fs/1e9, imag(Ste.s11), '--');
%     plot(axS, fs/1e9, real(Stm.s11));
%     addlegendentry(axS, 'CSTTM');
%     plot(axS, fs/1e9, imag(Stm.s11), '--');
%     
%     xlabel(axS, 'Frequency [GHz]');
%     ylabel(axS, '|\Gamma| [dB]');
%     ylim(axS, [-1 1]);
%     figS.Name = 'CST';
% Plot results.
% [figAbsS, axAbsS] = figureex;
%     repeatcolormap(axAbsS, 2);
    plot(axAbsS, fs/1e9, abs(Ste.s11));
    addlegendentry(axAbsS, 'TE CST');
    plot(axAbsS, fs/1e9, abs(Ste.s12), '--');
    plot(axAbsS, fs/1e9, abs(Stm.s11));
    addlegendentry(axAbsS, 'TM CST');
    plot(axAbsS, fs/1e9, abs(Stm.s12), '--');
    
    xlabel(axAbsS, 'Frequency [GHz]');
    ylabel(axAbsS, '|\Gamma|');
    ylim(axAbsS, [0 1]);
    
%     figAbsS.Name = 'CST';
%     drawnow;






































