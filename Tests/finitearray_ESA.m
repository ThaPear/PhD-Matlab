% close all;
closewaitbars;
clearvars -except array;
SetupPath;
clear global;

f0 = 31e9;
c0 = Constants.c0;
l0 = c0/f0;
z0 = Constants.z0;

fs = (12:2:32) * 1e9;

th = eps * pi/180;
ph = 0 * pi/180;

% z1 = 90;
% z2 = z0;
% zfeed = 80;
% 
% dx = 4.35e-3; %4.35e-3; %4.354838709677e-3; %c0/f0 * 0.45;
% dy = dx;
% erback = 2.2;
% hback = 0.9e-3;
% hgap = -0.1e-3;
% 
% dslot = 0.3e-3;
% wfeed = 0.3e-3;
% wslot = 0.8e-3;
% 
dedge = 0.25*l0;
% walled = 1;
% 
% C = 0.2e-12;


% p = dx / 2;
% gamma = 0.2;
% N = 2;
% f0match = 19e9;
% f0design = 29e9;
% slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);
% 
% % Apply the hgap parameter.
% slab.elements{1}.elements{1}.L = slab.elements{1}.elements{1}.L + hgap;


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

% Number of unit cells.
Nx = 5;
Ny = 5;

excitation = ones(Nx, Ny);

% tlineup = TerminatedTLine(slab, FreeSpace());
% tlinedown = ShortedLine(erback*0.7, hback);

slot = Slot(dx, dy, wslot, dslot, walled);

th = eps;
ph = 0;
infarray = InfiniteArray(slot, tlineup, tlinedown);
Zas = infarray.GetInputImpedance(fs, th, ph);
figureex; plot(fs/1e9, real(Zas), fs/1e9, imag(Zas));
% 
% arrayx = FiniteArrayX(slot, tlineup, tlinedown, Nx, excitation(:,1).', dedge, zfeed);
% Zas = arrayx.GetInputImpedance(fs, th, ph);
% figureex; plot(fs, real(Zas), fs, imag(Zas));
% 
% arrayy = FiniteArrayY(slot, tlineup, tlinedown, Ny, excitation(1,:), zfeed);
% Zas = arrayy.GetInputImpedance(fs, th, ph);
% figureex; plot(fs, real(Zas), fs, imag(Zas));


if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
    array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
end

tc = tic;
array.InitializeDs(fs);
save(sprintf('H:\\Git\\PhD-Matlab\\Tests\\%s\\%ix%i_D', mfilename, Nx, Ny), 'array', '-v7');

array.InitializeKyInts(fs);
save(sprintf('H:\\Git\\PhD-Matlab\\Tests\\%s\\%ix%i_Ky', mfilename, Nx, Ny), 'array', '-v7');

array.InitializeZMatrix(fs);
save(sprintf('H:\\Git\\PhD-Matlab\\Tests\\%s\\%ix%i_Zmat', mfilename, Nx, Ny), 'array', '-v7');

Zas = array.GetInputImpedance(fs, excitation);

Zcap = 1 ./ (1j .* 2.*pi.*fs .* C);
for(fi = 1:length(fs))
    ZasC(:,:,fi) = Zas(:,:,fi) + Zcap(fi); % With series capacitance.
end

Sact = (ZasC - zfeed) ./ (ZasC+ zfeed);

dt_ML = toc(tc);
dispex('MATLAB took %.1fs for %ix%i elements.\n', dt_ML, Nx, Ny);

%% 
hFig = figureex;
for(nx = 1:Nx)
    for(ny = 1:Ny)
        hAx = subplot(Nx, Ny, ny+(nx-1)*Ny);
        hold(hAx, 'on');
        grid(hAx, 'on');
        box(hAx, 'on');
        plot(hAx, fs/1e9, squeeze(20*log10(abs(Sact(nx, ny, :)))), 'k');
        ylim(hAx, [-30 0]);
    end
end
% legend(hAx, {'real', 'imag'});

%% Sactive
hFig = figureex;
for(nx = 1:Nx)
    for(ny = 1:Ny)
        hAx = subplot(Nx, Ny, ny+(nx-1)*Ny);
        hold(hAx, 'on');
        grid(hAx, 'on');
        box(hAx, 'on');
        plot(hAx, fs/1e9, squeeze(real(Zas(nx, ny, :))), 'r');
        plot(hAx, fs/1e9, squeeze(imag(Zas(nx, ny, :))), 'r--');
    end
end
legend(hAx, {'real', 'imag'});

[hFig2, hAx2] = figureex;
[hFig3, hAx3] = figureex;
for(nx = 1:Nx)
    for(ny = 1:Ny)
        plot(hAx2, fs/1e9, squeeze(real(Zas(nx, ny, :))), 'r');
        plot(hAx3, fs/1e9, squeeze(imag(Zas(nx, ny, :))), 'r');
    end
end

% Real active input impedance
[hFig, hAx] = figureex;
for(nx = 1:Nx)
    hAx = subplot(6,6,nx);
    hold(hAx, 'on');
    title(hAx, sprintf('nx = %i', nx));
    xlim(hAx, [-inf, inf]);
    ylim(hAx, [0 250]);
    hAx.Colormap = jet(Ny);
    for(ny = 1:Ny)
        plot(hAx, fs/1e9, squeeze(real(Zas(nx, ny, :))));
    end
end
% Imaginary impedance
[hFig, hAx] = figureex;
for(nx = 1:Nx)
    hAx = subplot(6,6,nx);
    hold(hAx, 'on');
    title(hAx, sprintf('nx = %i', nx));
    xlim(hAx, [-inf, inf]);
    ylim(hAx, [-100, 150]);
    hAx.ColorMap = jet(Ny);
    for(ny = 1:Ny)
        plot(hAx, fs/1e9, squeeze(imag(Zas(nx, ny, :))));
    end
end

[hFig4, hAx4] = figureex;
    plot(hAx4, fs/1e9, squeeze(real(Zas(5,5, :))), 'k');
    plot(hAx4, fs/1e9, squeeze(imag(Zas(5,5,:))), 'r');
    
[hFig5, hAx5] = figureex;
for(nx = 1:Nx)
    for(ny = 1:Ny)
        plot(hAx5, fs/1e9, squeeze(20*log10(abs(Sact(nx,ny,:)))), 'r');
    end
end


% savefig(hFig, sprintf('H:\\Git\\PhD-Matlab\\Tests\\%s\\figures\\%ix%i', mfilename, Nx, Ny));

% clear hFig hAx;