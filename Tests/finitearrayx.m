% close all;
closewaitbars;
clear;
SetupPath;
clear global;

f0 = 31e9;
c0 = Constants.c0;
l0 = c0/f0;

dx = 0.45*l0;
dy = dx;
wslot = 0.05*l0;
dslot = 0.05*l0;
walled = 0;
dedge = 0.25*l0;
zfeed = 100;

% Number of unit cells.
Nx = 1;

fs = (12:1:52)*1e9;

ax = ones(1,Nx);

tlineup = FreeSpace();
tlinedown = FreeSpace();

slot = Slot(dx, dy, wslot, dslot, tlineup, tlinedown, walled);

array = FiniteArrayX(slot, Nx, ax, dedge, zfeed);

th = eps * pi/180;
ph = 0 * pi/180;

tic;
Zas = array.GetInputImpedance(fs, th, ph);
toc
%%
[hFig, hAx] = figureex;
    hAx.ColorOrder = lines(min(Nx, 7));%reshape(repmat(lines(7), 1, 2).', [], 14).';
    hAx.LineStyleOrder = {'-', '--', ':', '-.'};
    plot(hAx, fs./1e9, real(Zas));
    addlegendentry(hAx, 'Real');
    plot(hAx, fs./1e9, imag(Zas), '--');
    addlegendentry(hAx, 'Imag');
%%
% global kxpath
% figureex;
%     plot(real(kxpath), imag(kxpath), '.');
%     k0s = 2*pi.*fs./c0;
%     plot([k0s;-k0s], zeros(size(k0s)), 'x');