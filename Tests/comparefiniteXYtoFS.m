close all;
clear;
SetupPath;
% clear global;
f0 = 31e9;
c0 = Constants.c0;
l0 = c0/f0;
z0 = Constants.z0;

f = 12.5 * 1e9;

th = eps * pi/180;
ph = 0 * pi/180;

dx = 0.45*l0;
dy = dx;
wslot = 0.05*l0;
dslot = 0.05*l0;
walled = 0;
dedge = 0.25*l0;
zfeed = 100;


Nx = 1;
Ny = 1;
excitation = ones(Nx, Ny);

tlineup = FreeSpace();
tlinedown = FreeSpace();

slot = Slot(dx, dy, wslot, dslot, walled);

array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);

array.InitializeDs(f);

[k0, ~, ~, ~] = k(f, 1, th, ph);

%% Determine kxs
% Go to +-inf
delta = 0.01.*k0;
lim1 = -100.*k0-1j.*delta;
lim2 = -lim1;
% Deform integration path around branch cuts.
integrationpath = [(-1-1j).*delta, (1+1j).*delta];

% Integration path (dotted)
%               |
%               |       N1
%               |...............
% --------------:---------------
% .............:|\
%      N1       | N2
%               |

N1 = 1000; % Initial number of points on the horizontal sections
N2 = 100;  % Initial number of points on the diagonal section through the origin
kx = [linspace(lim1, integrationpath(1), N1), ...
         linspace(integrationpath(1), integrationpath(2), N2+2), ...
         linspace(integrationpath(2), lim2, N1)];

% Remove the duplicate elements
kx(N1) = [];
kx(end-N1) = [];
%% 

% Equation 2.19 from Daniele's thesis
D = 1./(z0*k0) .* (k0.^2-kx.^2) .* besselj(0,    wslot/4 .* sqrt(k0.^2-kx.^2)) .* ...
                                   besselh(0, 2, wslot/4 .* sqrt(k0.^2-kx.^2));

nyp = 0; % Self of the first (and only) slot
deformedpath = 0; % No deformation
fi = array.D_fs == f;
Ds = array.D_interpolants{deformedpath+1, fi, nyp+1}(real(kx));

[hFig, hAx] = figureex;
    repeatcolormap(hAx, 2);
    plot(real(kx)./k0, real(D));
    addlegendentry('analytical');
    plot(real(kx)./k0, imag(D), '--');

    plot(real(kx)./k0, real(Ds));
    addlegendentry('finitearray');
    plot(real(kx)./k0, imag(Ds), '--');
    
[hFig, hAx] = figureex;
    repeatcolormap(hAx, 2);
    title(hAx, '% Difference');
    plot(real(kx)./k0, 100*abs(real(D)-real(Ds))./abs(real(D)));
    plot(real(kx)./k0, 100*abs(imag(D)-imag(Ds))./abs(imag(D)), '--');
[hFig, hAx] = figureex;
    repeatcolormap(hAx, 2);
    title(hAx, '% Difference');
    title(hAx, 'Abs Difference');
    plot(real(kx)./k0, abs(real(D)-real(Ds)));
    plot(real(kx)./k0, abs(imag(D)-imag(Ds)), '--');
    
    
    
