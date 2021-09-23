
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
erback = 1;
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

% tlineup = FreeSpace();
tlinedown = ShortedLine(erback, hback, erback);
% tlinedown = FreeSpace();

slot = Slot(dx, dy, wslot, dslot, walled);

%%
fs = (12:0.1:32)*1e9;

th = 60*pi/180+eps;
ph = 90*pi/180;

project = CST.InitializePeriodicProject();

array = InfiniteArray(slot, tlineup, tlinedown);
array = InfiniteArray(slot, FreeSpace(), tlinedown);
array.BuildCST(project);

project.StoreParameter('aa_theta', round(th*180/pi));
project.StoreParameter('aa_phi', ph*180/pi);
project.StoreParameter('nsamplesperGHz', 1);

project.Rebuild();

fdsolver = project.FDSolver();
if(~fdsolver.Start())
    error('Simulation failed');
end

touchstonepath = sprintf('%s/Validations/%s/', resultdir_matlab, mfilename);

filename = sprintf('th%i_ph%i', round(th*180/pi), round(ph*180/pi));
touchstonefilepath = [touchstonepath, filename];

touchstone = project.TOUCHSTONE();
touchstone.Reset();
touchstone.Impedance(zfeed);
touchstone.Renormalize(1);
touchstone.FileName(touchstonefilepath);
touchstone.Write();


[parameters, SCST] = CST.LoadData(sprintf('%s.s1p', touchstonefilepath));
fsCST = parameters.frequencies;
ZasCST = squeeze((1+SCST)./(1-SCST).*zfeed);

Zas = array.GetInputImpedance(fs, th, ph);

figureex;
    repeatcolormap(2);
    plot(fs./1e9, real(Zas));
    addlegendentry('MATLAB');
    plot(fs./1e9, imag(Zas), '--');
    
    plot(fsCST./1e9, real(ZasCST));
    addlegendentry('CST');
    plot(fsCST./1e9, imag(ZasCST), '--');
    