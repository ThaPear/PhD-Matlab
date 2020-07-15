close all;
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

zfeed = 100;

% Number of unit cells.

Ny = 5;
ay = ones(1,Ny);

tlineup = FreeSpace();
tlinedown = FreeSpace();

slot = Slot(dx, dy, wslot, dslot, tlineup, tlinedown, walled);

array = FiniteArrayY(slot, Ny, ay, zfeed);

fs = (5:10)*1e9;

th = eps * pi/180;
ph = 0 * pi/180;

tic;
Zas = array.GetInputImpedance(fs, th, ph);
toc

[hFig, hAx] = figureex;
    hAx.ColorOrder = lines(Ny);%reshape(repmat(lines(7), 1, 2).', [], 14).';
    plot(hAx, fs./1e9, real(Zas));
    plot(hAx, fs./1e9, imag(Zas), '--');