close all;
clear;
SetupPath;

f0 = 10e9;
c0 = Constants.c0;
l0 = c0 / f0;

z0 = Constants.z0;
zfeed = 180;
ery = (z0 / sqrt(z0 * zfeed))^2;
er = [ery ery 1];
er = [3 3 3];
mu = 1;

% slab = Line(l0/4/sqrt(max(er)), er, mu);
% slab = TLine({Line(l0/4/sqrt(2), [2 2 2], mu), Line(l0/4/sqrt(3), [3 3 3], mu)});

slab = TLine({Line(l0/4/sqrt(2), [2], mu), Line(l0/4/sqrt(3), [3], mu), Line(l0/4/sqrt(4), [4], mu)});

fs = 1e9 * (0.1:0.1:20);
th = 45 * pi/180;
ph = 0 * pi/180;

% Calculate propagation constants.
[k0, kx0, ky0, kz0] = k(fs, 1, th, ph);
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
    
    addlegendentry(axS, {'TE', 'TM'});
    xlabel(axS, 'Frequency [GHz]');
    ylabel(axS, '|\Gamma| [dB]');
    ylim(axS, [-30 0]);
    figS.Name = 'Analytical';
    
%% Simulate in CST.
mu = repmat(mu, 1, (4-length(mu)));
filename = sprintf('%s/Validations/%s/eps_%.1f_%.1f_%.1f_mu_%.1f_%.1f_%.1f_th%i', resultdir_matlab, mfilename, er(1), er(2), er(3), mu(1), mu(2), mu(3), th*180/pi);
touchstonefile = [filename, '.s4p'];
if(~exist(touchstonefile, 'file'))
    project = CST.InitializeUnitCellProject();
    project.StoreParameter('aa_theta', round(th * 180/pi));
    project.StoreParameter('aa_phi', round(ph * 180/pi));
    project.StoreParameter('fmin', 1);
    project.StoreParameter('fmax', 20);
    project.StoreParameter('fmesh', 20);
    project.StoreParameter('dx', l0/4*1e3);
    project.StoreParameter('dy', 'dx');
    project.StoreParameter('nsamplesperGHZ', 1);

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
    touchstone.FileName(filename);
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
ABCDte = S2ABCD(Ste, z0, z0);                           ABCDtm = S2ABCD(Stm, z0, z0);
Ste = ABCD2S(ABCDte, zfeed, z0te(1));                   Stm = ABCD2S(ABCDtm, zfeed, z0tm(1));

% [figS, axS] = figureex(2);
    plot(axS, fs/1e9, 20*log10(abs(Ste.s11)), 'x');
    plot(axS, fs/1e9, 20*log10(abs(Stm.s11)), 'x');
    
    addlegendentry(axS, {'CSTTE', 'CSTTM'});
    xlabel(axS, 'Frequency [GHz]');
    ylabel(axS, '|\Gamma| [dB]');
    ylim(axS, [-30 0]);
    figS.Name = 'CST';






































