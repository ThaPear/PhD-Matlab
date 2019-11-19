% close all;
clear;
SetupPath;
% clc;

%% Structure setup.
c0 = Constants.c0;
z0 = Constants.z0;
f0 = 31e9;          % 13.75 to 31 GHz.

z1 = 120;
z2 = z0;

% Slot parameters.
dx = 4.4e-3; %4.35e-3; %4.354838709677e-3; %c0/f0 * 0.45;
dy = dx;
erback = Materials.Substrate.permittivity;
hback = 2e-3;
wslot = 1.8e-3;
dslot = 2.5e-3;
walled = 1;

% Feed settings.
C = inf;%1e-12;
zfeed = 120;

% ADL parameters.
p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 20e9;
f0design = 29e9;


% baseheight = slab.elements{1}.elements{1}.L;
global hcut;
hcut = 0;
% for(hcut = (0:0.02:0.2)*1e-3)
% for(wslot = (0.5:0.5:3)*1e-3)
slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);
% Transmission line models for up and down stratifications.
% if(baseheight - hcut < 0)
%     continue;
% end
% slab.elements{1}.elements{1}.L = baseheight - hcut;
tlineup = TerminatedTLine(slab, FreeSpace());
tlinedown = ShortedLine(erback*0.7, hback, erback);

% The slot.
% slot = Slot(dx, dy, wslot, dslot, tlineup, tlinedown);
slot = Slot_Dualpol_Bowtie(dx, dy, wslot, dslot, tlineup, tlinedown);
slot = Slot(dx, dy, wslot, dslot, tlineup, tlinedown);
% return

fprintf('hback = %g, wslot = %g, dslot = %g, z1 = %g, zfeed = %g, f0design = %g, f0match = %g\n', ...
    hback*1e3, wslot*1e3, dslot*1e3, z1, zfeed, f0design/1e9, f0match/1e9);

% PreliminaryDesign_SimADL;
% PreliminaryDesign_SimSlot;
PreliminaryDesign_SimCST;
% end