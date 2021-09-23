clear;

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
% 
tlineup = TerminatedTLine(slab, FreeSpace());

% tlineup = FreeSpace();
tlinedown = ShortedLine(erback, hback, erback);
% tlinedown = FreeSpace();

slot = Slot(dx, dy, wslot, dslot, walled);

%%
fs = (12:2:32)*1e9;
% fs = 35e9;

th = 0*pi/180+eps;
ph = 0*pi/180;

Ny = 5;
ay = ones(1,Ny);
zfeed = 80;

array = FiniteArrayY(slot, tlineup, tlinedown, Ny, ay, zfeed);

Zas = array.GetInputImpedance(fs, th, ph);

%% Plot CST and Matlab results
% ymax = ceil(max([real(Zas(:)); real(ZasCST(:)); imag(Zas(:)); imag(ZasCST(:))])/100)*100;
% ymin = floor(min([real(Zas(:)); real(ZasCST(:)); imag(Zas(:)); imag(ZasCST(:))])/100)*100;
for(ny = 1:Ny)
    if(mod(ny, 5) == 1)
        [hFig, hAx] = figureex;
            repeatcolormap(hAx, 2);
    end
    hAx = subplot(min(Ny, 5), 1, mod(ny-1, 5)+1);
        hold(hAx, 'on');
        grid(hAx, 'on');
        box(hAx, 'on');
        title(hAx, sprintf('%i', ny));
    plot(hAx, fs./1e9, real(Zas(ny, :)), 'k');
    addlegendentry(hAx, 'MATLAB');
    plot(hAx, fs./1e9, imag(Zas(ny, :)), 'k--');
    
%     ylim(hAx, [ymin ymax]);
    ylim(hAx, [-100 200]);
    xlim(hAx, [-inf inf]);
end
