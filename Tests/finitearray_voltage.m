% close all;
clearvars -except array;
SetupPath;

f0 = 20e9;
f = 15e9;
l0 = 3e8/f0;
[k0, ~, ~, ~] = k(f0, 1, 0, 0);

% dx = 0.45*l0;
% dy = dx;
% wslot = 0.1*l0;
% dslot = 0.1*l0;
% dedge = 0.25*l0;

% zfeed = 80;

dx = 0.45*l0;
dy = dx;
wslot = 0.05*l0;
dslot = 0.05*l0;
walled = 0;
dedge = 0.25*l0;
zfeed = 100;

Nx = 2;
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


% x = linspace(-0.02, 0.02, 2000);
% v = zeros(size(x));
% for(xi = 1:length(x))
%     v(xi) = integral(@(kx) sinc(kx .* dslot ./ (2*pi)) .* exp(-1j .* kx .* x(xi)),-100*k0, 100*k0);
% end

[hFig, hAx] = figureex;
hAx.ColorOrder = reshape(repmat(lines(7), 1, 3).', [], 21).';
for(ny = 0:Ny-1)
    tc = tic;
    [x, v] = array.Voltage(f, excitation, ny);
    dt = toc(tc);
    dispex('Finished voltage on slot %i in %.2fs.\n', ny, dt);
    plot(hAx, x, real(v));
    addlegendentry(sprintf('ny = %i', ny));
    plot(hAx, x, imag(v), '--');
    plot(hAx, x, abs(v), ':');
    
end