% close all;
clearvars -except array;
SetupPath;

f0 = 10e9;
f = 10e9;
l0 = 3e8/f0;
[k0, ~, ~, ~] = k(f0, 1, 0, 0);

dx = 0.45*l0;
dy = dx;
wslot = 0.1*l0;
dslot = 0.1*l0;
dedge = 0.25*l0;

zfeed = 80;

ths = linspace(-89, 89, 200) * pi/180;
phs = 0 * pi/180;

Nxs = 5;
for(Nx = Nxs)
Ny = 5;
excitation = ones(Nx, Ny);

tlineup = FreeSpace();
tlinedown = FreeSpace();
% tlinedown = ShortedLine(1, l0/4);

unitcell = Slot(dx, dy, wslot, dslot, 0);

if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
    array = FiniteArray(unitcell, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
end

array.InitializeDs(f);
array.InitializeZMatrix(f);

ff = array.NormalizedFarfield(f, excitation, ths, phs);



[hFig, hAx] = figureex;
hAx.ColorOrder = reshape(repmat(lines(7), 1, 3).', [], 3*7).';

if(length(phs) < 2)
    plot(hAx, ths*180/pi, 20*log10(abs(ff)), ':');
    title(hAx, sprintf('Nx %i, Ny %i, ph %.1f\\circ', Nx, Ny, phs*180/pi));
    ylim(hAx, [-30 0]);
    xlim(hAx, [-inf inf]);
else
    dispex('Plotting over 1 phi not supported.\n');
end

end

return;

project = CST.InitializeBasicProject();
project.StoreParameter('fmin', min(f)/1e9);
project.StoreParameter('fmax', max(f)/1e9);
project.StoreParameter('fmesh', 'fmax');
project.StoreParameter('slot_impedance', zfeed);
project.StoreParameter('nsamplesperGHz', 1);

boundary = project.Boundary();
boundary.Zmin('expanded open');
boundary.Zmax('expanded open');

meshadaption3d = project.MeshAdaption3D();
meshadaption3d.MinPasses(6);

array.BuildCST(project);