% close all;
closewaitbars;
clear;
SetupPath;
clear global;

f0 = 20e9;
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
Nx = 3;

fs = 20e9;

ax = ones(1,Nx);

tlineup = FreeSpace();
tlinedown = FreeSpace();

slot = Slot(dx, dy, wslot, dslot, walled);

array = FiniteArrayX(slot, tlineup, tlinedown, Nx, ax, dedge, zfeed);

th = eps * pi/180;
ph = 0 * pi/180;

tic;
Zas = array.GetInputImpedance(fs, th, ph);
toc

[x, v] = array.Voltage(fs, th, ph);

[hFig, hAx] = figureex;
    repeatcolormap(hAx, 3);
    plot(hAx, x, real(v));
    plot(hAx, x, imag(v), '--');
    plot(hAx, x, abs(v), ':');
    
%%
% global kxpath
% figureex;
%     plot(real(kxpath), imag(kxpath), '.');
%     k0s = 2*pi.*fs./c0;
%     plot([k0s;-k0s], zeros(size(k0s)), 'x');