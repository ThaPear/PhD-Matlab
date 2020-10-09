close all;
clear;
SetupPath;

f0 = 10e9;
l0 = Constants.c0/f0;
z0 = Constants.z0;

er = [2 4 2];
mu = [1 1 1];

L = 0.25.*l0./sqrt(max(er));

th0 =  eps * pi/180;
th  = 20 * pi/180;
ph = eps * pi/180;

[k0, ~, ~, kzd] = k(f0, 1, th, ph);
[~, z0te, z0tm] = z(1, k0, kzd);

% slab = Line(L, er, mu);

p = 0.2*l0;
ws = repmat(0.01*l0, 1, 2);
ds = 0.02*l0;
% ds = [0, ds, 0];
ds = [ds(1)/2, ds, ds(end)/2];
ss = repmat(p/2, 1, length(ws));
slab = ADS(p, ds, ss, ws, 1);

DrawADS(slab, f0);

% filename = sprintf('h:/Git/PhD-Matlab/Validations/%s/eps_%.1f_%.1f_%.1f_mu_%.1f_%.1f_%.1f', mfilename, er(1), er(2), er(3), mu(1), mu(2), mu(3));
filename = sprintf('h:/Git/PhD-Matlab/Validations/%s/adltest2', mfilename);
%{
if(~exist([filename, '_th0.s4p'], 'file') || ~exist([filename, '_th.s4p'], 'file'))
    project = CST.InitializeUnitCellProject();

    project.StoreParameter('aa_phi', ph * 180/pi);
    project.StoreParameter('f0', f0/1e9);
    project.StoreParameter('fmesh', f0/1e9);
    project.StoreParameter('l0', 'c0/f0/1e9');
    project.StoreParameter('dx', '0.2*l0');
    project.StoreParameter('dy', 'dx');
    project.StoreParameter('openboundary_distance', 'l0/4')

    fdsolver = project.FDSolver();
    fdsolver.ResetSampleIntervals('all');
    fdsolver.AddSampleInterval('f0', 'f0', 1, 'Single', 1);

    slab.BuildCST(project);
    
    if(~exist([filename, '_th0.s4p'], 'file'))
        project.StoreParameter('aa_theta', th0 * 180/pi);
        project.Rebuild();
        fdsolver.Start();
            touchstone = project.TOUCHSTONE();
            touchstone.Reset();
            touchstone.Impedance(z0);
            touchstone.Renormalize(1);
            touchstone.FileName([filename, '_th0']);
            touchstone.Write();
    end

    if(~exist([filename, '_th.s4p'], 'file'))
        project.StoreParameter('aa_theta', th * 180/pi);
        project.Rebuild();
        fdsolver.Start();
            touchstone = project.TOUCHSTONE();
            touchstone.Reset();
            touchstone.Impedance(z0);
            touchstone.Renormalize(1);
            touchstone.FileName([filename, '_th']);
            touchstone.Write();
    end
%     project.Quit();
end

[parameters, data] = CST.LoadData([filename, '_th0.s4p']);
Ste_th0 = []; Ste_th0.s11 = data(1, 1); Ste_th0.s12 = data(1, 3); Ste_th0.s21 = data(3, 1); Ste_th0.s22 = data(3, 3);
Stm_th0 = []; Stm_th0.s11 = data(2, 2); Stm_th0.s12 = data(2, 4); Stm_th0.s21 = data(4, 2); Stm_th0.s22 = data(4, 4);


[parameters, data] = CST.LoadData([filename, '_th.s4p']);
Ste_th = []; Ste_th.s11 = data(1, 1); Ste_th.s12 = data(1, 3); Ste_th.s21 = data(3, 1); Ste_th.s22 = data(3, 3);
Stm_th = []; Stm_th.s11 = data(2, 2); Stm_th.s12 = data(2, 4); Stm_th.s21 = data(4, 2); Stm_th.s22 = data(4, 4);
% Renormalize to z0te and z0tm.
ABCDte = S2ABCD(Ste_th, z0, z0);  Ste_th = ABCD2S(ABCDte, z0te, z0te);
ABCDtm = S2ABCD(Stm_th, z0, z0);  Stm_th = ABCD2S(ABCDtm, z0tm, z0tm);

L = slab.GetHeight();
[er, mu] = EpsilonMu(f0, Ste_th, Stm_th, Ste_th0, Stm_th0, L, th0, th, ph);

dispex('er = [%.2f, %.2f, %.2f], mu = [%.2f, %.2f, %.2f]\n', er.x, er.y, er.z, mu.x, mu.y, mu.z);
%}



[k0, kx, ky, ~] = k(f0, 1, th, ph);
kr = sqrt(kx.^2 + ky.^2);
[k0, ~, ~, kzd] = k(f0, 1, th, ph);
[~, z0te, z0tm] = z(1, k0, kzd);
ABCDte = slab.GetABCD(1, f0, k0, kr);
Ste_th = ABCD2S(ABCDte, z0te, z0te);
ABCDtm = slab.GetABCD(0, f0, k0, kr);
Stm_th = ABCD2S(ABCDtm, z0tm, z0tm);


[k0, kx, ky, ~] = k(f0, 1, eps, ph);
kr = sqrt(kx.^2 + ky.^2);
[k0, ~, ~, kzd] = k(f0, 1, eps, ph);
[~, z0te, z0tm] = z(1, k0, kzd);

ABCDte = slab.GetABCD(1, f0, k0, kr);
Ste_th0 = ABCD2S(ABCDte, z0te, z0te);
ABCDtm = slab.GetABCD(0, f0, k0, kr);
Stm_th0 = ABCD2S(ABCDtm, z0tm, z0tm);

L = slab.GetHeight();


[er, mu] = EpsilonMu(f0, Ste_th, Stm_th, Ste_th0, Stm_th0, L, th0, th, ph);

dispex('er = [%.2f, %.2f, %.2f], mu = [%.2f, %.2f, %.2f]\n', er.x, er.y, er.z, mu.x, mu.y, mu.z);

