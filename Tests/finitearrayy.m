% close all;
closewaitbars;
clear;
SetupPath;
clear global;

%% Test parameters
%{
f0 = 31e9;
c0 = Constants.c0;
l0 = c0/f0;

dx = 0.45*l0;
dy = dx;
wslot = 0.05*l0;
dslot = 0.05*l0;
walled = 0;

zfeed = 100;

tlineup = FreeSpace();
tlinedown = FreeSpace();

fs = (5:10)*1e9;
%}

%% ESA parameters
%{*
c0 = Constants.c0;
z0 = Constants.z0;
f0 = 31e9;          % 13.75 to 31 GHz.

z1 = 80;
z2 = z0;
zfeed = 80;

dx = 4.35e-3; %4.35e-3; %4.354838709677e-3; %c0/f0 * 0.45;
dy = dx;
erback = 2.2;
hback = 1.9e-3;%1.9496e-3; % = (c0/f0)/sqrt(erback*0.7)/4
wslot = 1.4e-3;
dslot = 2e-3;
walled = 1;

C = inf;%0.2e-12;

p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 19e9;
f0design = 29e9;
slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);

tlineup = TerminatedTLine(slab, FreeSpace());
tlinedown = ShortedLine(erback, hback, erback);

fs = (12:2:32)*1e9;
%}
 
%%
slot = Slot(dx, dy, wslot, dslot, walled);

% Number of unit cells.
Ny = 2;
ay = ones(1,Ny);

array = FiniteArrayY(slot, tlineup, tlinedown, Ny, ay, zfeed);

th = eps * pi/180;
ph = 0 * pi/180;

tic;
Zas = array.GetInputImpedance(fs, th, ph);
toc

[hFig, hAx] = figureex;
    hAx.ColorOrder = lines(Ny);%reshape(repmat(lines(7), 1, 2).', [], 14).';
    plot(hAx, fs./1e9, real(Zas));
    plot(hAx, fs./1e9, imag(Zas), '--');
























