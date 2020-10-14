close all;



GHz  = 10^9;
f0 = 8*GHz;%7.125*GHz;
fADL = 8.2*GHz;%6.775*GHz;
l0 = physconst('Lightspeed')/f0;
lADL = physconst('Lightspeed')/fADL;

%Array
dx = 0.4*l0;%2*0.190175*lADL; %Spacing between feeding gaps in the x-dimension
dy = dx; %Spacing between slots in the y-dimension
ws = 0.1*l0;%0.21*l0;%0.21*l0; %Width of the slot
deltas = 0.1*l0;%0.24*l0;%0.24*l0; %Length of the gap
h = 0.1*l0; %Distance between the array plane and the backing reflector

%ADL
NLayers = 4; %Number of layers
p =  0.205*lADL;%0.190175*lADL; %Period
w = [0.0265 0.0265 0.089 0.089]*lADL; %Gaps between metal patches
d = [0.0462 (0.046+0.08971)/2 0.08971]*lADL;%0.057542*l0;

hGap = d(1)/2.0;%d(1)/2.45;% d(1)/2; %Distance between the array plane and the first ADL
s = [0.5 0.5 0.5 0.5]*p;%0.5*p*ones(NLayers-1); %Shift between layers
epsr = ones(NLayers+1);

ZLine = 70; %Feeding impedance
C = Inf;%1.5e-12;%1.5e-12; %Capacitor

%Scanning direction (in degrees)
ThetaV = atand(sqrt(2)*tand(45));
PhiV = 45;
Nfm = 20; %Floquet modes are 2*N+1

%%
slab = ADS(p, [hGap, d, 0], s, w, epsr);
br = ShortedLine(1, h);

tlineup = TerminatedTLine(slab, FreeSpace());
tlinedown = br;

slot = Slot(dx, dy, ws, deltas, 0);
array = InfiniteArray(slot, tlineup, tlinedown);

fs = (5:0.1:9)*1e9;

th = ThetaV * pi/180;
ph = PhiV * pi/180;

Zin = array.GetInputImpedance(fs, th, ph);

figureex;
    plot(fs, real(Zin));
    plot(fs, imag(Zin), '--');
    
% DrawADS(slab, fADL)