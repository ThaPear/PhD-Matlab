close all;
clear;
SetupPath;

f0 = 10e9;
c0 = Constants.c0;
l0 = c0 / f0;

z0 = Constants.z0;
zfeed = 180;
ery = (z0 / sqrt(z0 * zfeed))^2;
er = [8 ery 1];

slab = Line(l0/sqrt(ery)/4, er);

fs = 1e9 * (0.1:0.1:20);
th = eps * pi/180;
ph = 0 * pi/180;

% Calculate propagation constants.
[k0, kx0, ky0, kz0] = k(fs, 1, th, ph);
% [kd, ~, ~, kzd] = k(fs, zfeed, th, ph);
kr = sqrt(kx0.^2 + ky0.^2);
% Calculate impedances.
[z0, z0te, z0tm] = z(1, k0, kz0);
%     [zd, zdte, zdtm] = z((z0/z1)^2, kd, kzd);

%% Simulate analytically.
ABCDte = slab.GetABCD(1, fs, k0, kr);                   ABCDtm = slab.GetABCD(0, fs, k0, kr);
Ste = ABCD2S(ABCDte, zfeed, z0te);                      Stm = ABCD2S(ABCDtm, zfeed, z0tm);
Zte = ABCD2Z(ABCDte);                                   Ztm = ABCD2Z(ABCDtm);
Zinte = Zte.z11 - Zte.z12.*Zte.z21 ./ (Zte.z22+z0te);   Zintm = Ztm.z11 - Ztm.z12.*Ztm.z21 ./ (Ztm.z22+z0tm);
% Plot results.
[figS, axS] = figureex(1);
    plot(axS, fs/1e9, 20*log10(abs(Ste.s11)));
    plot(axS, fs/1e9, 20*log10(abs(Stm.s11)));
    
    legend(axS, {'TE', 'TM'});
    xlabel(axS, 'Frequency [GHz]');
    ylabel(axS, '|\Gamma| [dB]');
    ylim(axS, [-30 0]);
    figS.Name = 'Analytical';
    
%% Simulate in CST.
touchstonefilebase = ['h:/temp/temp'];
touchstonefile = [touchstonefilebase, '.s4p'];

if(~exist(touchstonefile, 'file'))
    project = CST.InitializeUnitCellProject();
    project.StoreParameter('aa_theta', round(th * 180/pi));
    project.StoreParameter('aa_phi', round(ph * 180/pi));
    project.StoreParameter('fmin', 1);
    project.StoreParameter('fmax', 20);
    project.StoreParameter('fmesh', 20);
    project.StoreParameter('dx', l0/4*1e3);
    project.StoreParameter('dy', 'dx');

    slab.BuildCST(project);
    project.Rebuild();
    
    % Run simulation.
    fdsolver = project.FDSolver();
    fdsolver.Start();

    % Export results.
    touchstone = project.TOUCHSTONE();
    touchstone.Reset();
    touchstone.Impedance(z0);
    touchstone.Renormalize(1);
    touchstone.FileName(touchstonefilebase);
    touchstone.Write();
end

% Import results.
[parameters, S] = CST.LoadData(touchstonefile);
fs = parameters.frequencies;
Ste = [];                                               Stm = [];
Ste.s11 = squeeze(S(1,1,:));                            Stm.s11 = squeeze(S(2,2,:));
Ste.s12 = squeeze(S(1,3,:));                            Stm.s12 = squeeze(S(2,4,:));
Ste.s21 = squeeze(S(3,1,:));                            Stm.s21 = squeeze(S(4,2,:));
Ste.s22 = squeeze(S(3,3,:));                            Stm.s22 = squeeze(S(4,4,:));

% Renormalize
ABCDte = S2ABCD(Ste, z0te(1), z0te(1));                 ABCDtm = S2ABCD(Stm, z0tm(1), z0tm(1));
Ste = ABCD2S(ABCDte, z0te(1), zfeed);                   Stm = ABCD2S(ABCDtm, z0tm(1), zfeed);

[figS, axS] = figureex(2);
    plot(axS, fs/1e9, 20*log10(abs(Ste.s11)));
    plot(axS, fs/1e9, 20*log10(abs(Stm.s11)));
    
    legend(axS, {'TE', 'TM'});
    xlabel(axS, 'Frequency [GHz]');
    ylabel(axS, '|\Gamma| [dB]');
    ylim(axS, [-30 0]);
    figS.Name = 'CST';






































