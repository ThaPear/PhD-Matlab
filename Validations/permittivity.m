% close all;
clear;
SetupPath;

f0 = 29e9;
l0 = Constants.c0/f0;
z0 = Constants.z0;


th0 =  eps * pi/180;
th  = 20 * pi/180;
ph = eps * pi/180;

[k0, ~, ~, kzd] = k(f0, 1, th, ph);
[~, z0te, z0tm] = z(1, k0, kzd);

style = 1; % Single point.
% style = 2; % Curve

%% Anisotropic slab.
%{
er = [2 4 2];
mu = [1 1 1];
L = 0.25.*l0./sqrt(max(er));
% slab = Line(L, er, mu);
% filename = sprintf('eps_%.1f_%.1f_%.1f_mu_%.1f_%.1f_%.1f', er(1), er(2), er(3), mu(1), mu(2), mu(3));
%}
%% ADL with given parameters.
%{
p = 0.2*l0;
ws = repmat(0.01*l0, 1, 2);
ds = 0.02*l0;
% ds = [0, ds, 0];
ds = [ds(1)/2, ds, ds(end)/2];
ss = repmat(p/2, 1, length(ws));
slab = ADS(p, ds, ss, ws, 1);
filename = sprintf('adltest');
%}
%% ADL with given epsilon.
%{
erdes = 8.4;
p = 4.35e-3/2;
L = 0.001348936979130;
slab = DesignADS(f0, p, L, erdes);
filename = sprintf('erdes_%g', erdes);
%}

%% Arbitrary stackup
glue = Line(38e-6, 2.32);
substrate = Line(25.4e-6, 3.4);
diel = Line(127e-6, 1.7);
slab = TLine({glue, diel, glue});
filename = sprintf('arbitstack_%i', floor(rand()*10000));


% DrawADS(slab, f0);
drawnow;

%% Determine path to place the CST files.
touchstonepath = sprintf('%s/Validations/%s/', resultdir_matlab, mfilename);
path = sprintf('%s/Validations/%s/', resultdir_matlab, mfilename);

cstfilepath = [path, filename];
touchstonefilepath = [touchstonepath, filename];

%% Validate using CST simulation.
%{*
if(~exist([touchstonefilepath, '_th0.s4p'], 'file') || ~exist([touchstonefilepath, '_th.s4p'], 'file'))
    project = CST.InitializeUnitCellProject();

    project.StoreParameter('aa_phi', ph * 180/pi);
    project.StoreParameter('f0', f0/1e9);
    project.StoreParameter('fmesh', f0/1e9);
    project.StoreParameter('l0', 'c0/f0/1e9');
    project.StoreParameter('dx', '0.4*l0');
    project.StoreParameter('dy', 'dx');
    project.StoreParameter('openboundary_distance', 'l0/4');
    
    meshadaption3d = project.MeshAdaption3D();
    meshadaption3d.MinPasses(6);

    fdsolver = project.FDSolver();
    if(style == 1) % Single point.
        fdsolver.ResetSampleIntervals('all');
        fdsolver.AddSampleInterval('f0', 'f0', 1, 'Single', 1);
    end
    
    wcs = project.WCS();
    wcs.Enable();

    slab.BuildCST(project);
    
    if(~exist([touchstonefilepath, '_th0.s4p'], 'file'))
        project.StoreParameter('aa_theta', th0 * 180/pi);
        project.Rebuild();
        tc = tic;
        fdsolver.Start();
            touchstone = project.TOUCHSTONE();
            touchstone.Reset();
            touchstone.Impedance(z0);
            touchstone.Renormalize(1);
            touchstone.FileName([touchstonefilepath, '_th0']);
            touchstone.Write();
        dispex('CST_th0 took %s.\n', fancyduration(toc(tc)));
    end

    if(~exist([touchstonefilepath, '_th.s4p'], 'file'))
        project.StoreParameter('aa_theta', th * 180/pi);
        project.Rebuild();
        tc = tic;
        fdsolver.Start();
            touchstone = project.TOUCHSTONE();
            touchstone.Reset();
            touchstone.Impedance(z0);
            touchstone.Renormalize(1);
            touchstone.FileName([touchstonefilepath, '_th']);
            touchstone.Write();
        dispex('CST_th took %s.\n', fancyduration(toc(tc)));
    end
%     project.Quit();
end
%% Import results.
[parameters, data] = CST.LoadData([touchstonefilepath, '_th0.s4p']);
Ste_th0 = []; Ste_th0.s11 = squeeze(data(1, 1, :)).'; Ste_th0.s12 = squeeze(data(1, 3, :)).'; Ste_th0.s21 = squeeze(data(3, 1, :)).'; Ste_th0.s22 = squeeze(data(3, 3, :)).';
Stm_th0 = []; Stm_th0.s11 = squeeze(data(2, 2, :)).'; Stm_th0.s12 = squeeze(data(2, 4, :)).'; Stm_th0.s21 = squeeze(data(4, 2, :)).'; Stm_th0.s22 = squeeze(data(4, 4, :)).';
fsCST = parameters.frequencies;


[parameters, data] = CST.LoadData([touchstonefilepath, '_th.s4p']);
Ste_th = []; Ste_th.s11 = squeeze(data(1, 1, :)).'; Ste_th.s12 = squeeze(data(1, 3, :)).'; Ste_th.s21 = squeeze(data(3, 1, :)).'; Ste_th.s22 = squeeze(data(3, 3, :)).';
Stm_th = []; Stm_th.s11 = squeeze(data(2, 2, :)).'; Stm_th.s12 = squeeze(data(2, 4, :)).'; Stm_th.s21 = squeeze(data(4, 2, :)).'; Stm_th.s22 = squeeze(data(4, 4, :)).';
% Renormalize to z0te and z0tm.
ABCDte = S2ABCD(Ste_th, z0, z0);  Ste_th = ABCD2S(ABCDte, z0te, z0te);
ABCDtm = S2ABCD(Stm_th, z0, z0);  Stm_th = ABCD2S(ABCDtm, z0tm, z0tm);

L = slab.GetHeight();
[er, mu] = EpsilonMu(fsCST, Ste_th, Stm_th, Ste_th0, Stm_th0, L, th0, th, ph);

%% Plot results.
if(style == 1) % Single point.
    dispex('erCST = [%.2f, %.2f, %.2f], mu = [%.2f, %.2f, %.2f]\n', er.x, er.y, er.z, mu.x, mu.y, mu.z);
else % Curve.
    [hFig, hAx] = figureex;
        hAx.ColorOrder = lines(3);
        hAx.LineStyleOrder = {'-', '--'};
        plot(fsCST./1e9, er.x);
        plot(fsCST./1e9, er.y);
        plot(fsCST./1e9, er.y);
        plot(fsCST./1e9, mu.x);
        plot(fsCST./1e9, mu.y);
        plot(fsCST./1e9, mu.y);
end
%}


%{*
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

[erML, muML] = EpsilonMu(f0, Ste_th, Stm_th, Ste_th0, Stm_th0, L, th0, th, ph);

dispex('er    = [%.2f, %.2f, %.2f], mu = [%.2f, %.2f, %.2f]\n', erML.x, erML.y, erML.z, muML.x, muML.y, muML.z);
%}
