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

Nx = 3;
Ny = 3;
excitation = ones(Nx, Ny);

tlineup = FreeSpace();
tlinedown = FreeSpace();
% tlinedown = ShortedLine(1, l0/4);

unitcell = Slot(dx, dy, wslot, dslot);

if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
    array = FiniteArray(unitcell, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
end

array.InitializeDs(f);
array.InitializeZMatrix(f);

kx = linspace(-10*k0, 10*k0, 101);
ky = zeros(size(kx));
ky(2) = inf;

V = array.VoltageSpectrum(f, excitation, kx, ky);



[hFig, hAx] = figureex;
hAx.ColorOrder = reshape(repmat(lines(7), 1, 3).', [], 3*7).';

plot(hAx, kx, real(V(1, 51)));
plot(hAx, kx, imag(V(1, 51)), '--');
plot(hAx, kx, abs(V(1, 51)), ':');

plot(hAx, ky, squeeze(real(V(1, 51, :))));
plot(hAx, ky, squeeze(imag(V(1, 51, :))), '--');
plot(hAx, ky, squeeze(abs(V(1, 51, :))), ':');
