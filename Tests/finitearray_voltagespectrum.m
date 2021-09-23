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

Nx = 1;
Ny = 1;
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

kx = linspace(-10*k0, 10*k0, 10001);
% ky = zeros(size(kx));
% ky(2) = inf;

ny = 0;

Vny = array.VoltageSpectrum(f, excitation, kx, ny);

[hFig, hAx] = figureex;
hAx.ColorOrder = reshape(repmat(lines(7), 1, 3).', [], 3*7).';

plot(hAx, kx, real(Vny));
plot(hAx, kx, imag(Vny), '--');
plot(hAx, kx, abs(Vny), ':');

% plot(hAx, ky, squeeze(real(M(1, 51, :))));
% plot(hAx, ky, squeeze(imag(M(1, 51, :))), '--');
% plot(hAx, ky, squeeze(abs(M(1, 51, :))), ':');
