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
walled = 0;
dedge = c0/f0/4; % Quarter-wave.

C = inf;%0.2e-12;

% p = dx / 2;
% gamma = 0.2;
% N = 2;
% f0match = 19e9;
% f0design = 29e9;
% slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);
% 
% tlineup = TerminatedTLine(slab, FreeSpace());

tlineup = FreeSpace();
% tlinedown = ShortedLine(erback, hback, erback);
tlinedown = FreeSpace();

slot = Slot(dx, dy, wslot, dslot, walled);

%%
fs = (12:2:32)*1e9;

th = eps*pi/180;
ph = 0*pi/180;

Nx = 5;
ax = ones(1,Nx);
zfeed = 80;

array = FiniteArrayX(slot, tlineup, tlinedown, Nx, ax, dedge, zfeed);

touchstonepath = sprintf('%s/Validations/%s/', resultdir_matlab, mfilename);

filename = sprintf('%i_th%i_ph%i', Nx, round(th*180/pi), round(ph*180/pi));
if(walled)
    filename = [filename, '_walled'];
end
touchstonefilepath = [touchstonepath, filename];

if(~exist(sprintf('%s.s%ip', touchstonefilepath, Nx), 'file'))
    dispex('Building CST simulation.\n');
    project = CST.InitializePeriodicProject();

    boundary = project.Boundary();
    boundary.Xmin('open');
    boundary.Xmax('open');
    boundary.Zmin('expanded open');

    array.BuildCST(project);

    project.StoreParameter('aa_theta', 0);%th*180/pi);
    project.StoreParameter('aa_phi', 0);%ph*180/pi);

    project.Rebuild();

    fdsolver = project.FDSolver();
    if(~fdsolver.Start())
        error('Simulation failed');
    end

    touchstone = project.TOUCHSTONE();
    touchstone.Reset();
    touchstone.Impedance(zfeed);
    touchstone.Renormalize(1);
    touchstone.FileName(touchstonefilepath);
    touchstone.Write();
else
    dispex('CST result already exists.\n');
end


[parameters, SCST] = CST.LoadData(sprintf('%s.s%ip', touchstonefilepath, Nx));
fsCST = parameters.frequencies;

SaCST = zeros(Nx, length(fsCST));
for(fi = 1:length(fsCST))
    f = fsCST(fi);
    [~, ~, ky0, ~] = k(f, 1, th, ph);
    
    ax = ones(1,Nx) .* exp(-1j .* ky0 .* (1:Nx) .* dy);

    for(nx = 1:Nx)
        SaCST(nx, fi) = sum(ax.' .* SCST(:, nx, fi) ./ ax(nx), 1);
    end
end

ZasCST = squeeze((1+SaCST)./(1-SaCST).*zfeed);

Zas = array.GetInputImpedance(fs, th, ph);
%%
ymax = ceil(max([real(Zas(:)); real(ZasCST(:)); imag(Zas(:)); imag(ZasCST(:))])/100)*100;
ymin = floor(min([real(Zas(:)); real(ZasCST(:)); imag(Zas(:)); imag(ZasCST(:))])/100)*100;
for(nx = 1:ceil(Nx/2))
    if(mod(nx, 5) == 1)
        [hFig, hAx] = figureex;
            delete(hFig.CurrentAxes);
            hAxes = TightSubplot(1, ceil(Nx/2), [0.01 0.01], [0.13 0.01], [0.1 0.05]);
            hAxes(end).XTickLabelMode = 'auto';
    end
    hAx = hAxes(mod(nx-1, 5)+1);%subplot(5, 1, mod(nx-1, 5)+1);
        hold(hAx, 'on');
        grid(hAx, 'on');
        box(hAx, 'on');
        hAx.XTickMode = 'auto';
        hAx.YTickLabelMode = 'auto';
%         title(hAx, sprintf('%i', nx));
    plot(hAx, fs./1e9, real(Zas(nx, :)), 'k');
    if(mod(nx, 5) == 0)
        addlegendentry(hAx, 'MATLAB');
    end
    plot(hAx, fs./1e9, imag(Zas(nx, :)), 'k--');

    plot(hAx, fsCST./1e9, real(ZasCST(nx, :)), 'r');
    if(mod(nx, 5) == 0)
        addlegendentry(hAx, 'CST');
    end
    plot(hAx, fsCST./1e9, imag(ZasCST(nx, :)), 'r--');
    xlim(hAx, [12 32]);
    ylim(hAx, [ymin ymax]);
    ylim(hAx, [-100 299]);
    xlabel(hAx, 'Frequency [GHz]');
    if(mod(nx, 2))
        ylabel(hAx, 'Input Impedance [\Omega]');
    end
end