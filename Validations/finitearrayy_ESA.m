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
% wslot = 1.4e-3;
wslot = 0.5e-3;
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
% tlineup = TerminatedTLine(TLine({Line(c0/f0/10, 1), slab}), FreeSpace()); % Extra space between ADL and slot

% tlineup = FreeSpace();
tlinedown = ShortedLine(erback, hback, erback);
% tlinedown = FreeSpace();

slot = Slot(dx, dy, wslot, dslot, walled);

%%
% fs = [0.1, 0.2, 0.3, 0.4, (1:2:32)]*1e9;
fs = (12:2:32)*1e9;
% fs = (24:2:32)*1e9;

th = 0*pi/180+eps;
ph = 0*pi/180;

Ny = 10;
ay = ones(1,Ny);
zfeed = 80;

array = FiniteArrayY(slot, tlineup, tlinedown, Ny, ay, zfeed);

touchstonepath = sprintf('%s/Validations/%s/', resultdir_matlab, mfilename);

filename = sprintf('%i_th%i_ph%i', Ny, round(th*180/pi), round(ph*180/pi));
if(walled)
    filename = [filename, '_walled'];
end
touchstonefilepath = [touchstonepath, filename];
%%
if(~exist(sprintf('%s.s%ip', touchstonefilepath, Ny), 'file'))
    dispex('Building CST simulation.\n');
    project = CST.InitializePeriodicProject();

    boundary = project.Boundary();
    boundary.Ymin('open');
    boundary.Ymax('open');
    boundary.Zmin('expanded open');

    array.BuildCST(project);

    project.StoreParameter('aa_theta', th*180/pi);
    project.StoreParameter('aa_phi', ph*180/pi);
    
    project.StoreParameter('fmin', 12);
    project.StoreParameter('fmax', 32);
    project.StoreParameter('nsamplesperGHz', 0.5);
    
    meshadaption3d = project.MeshAdaption3D();
    meshadaption3d.MinPasses(7);
    
    % Align ADL
    project.StoreParameter('adl_s0', -1.5);

    project.Rebuild();

    fdsolver = project.FDSolver();
    if(~fdsolver.Start())
        error('Simulation failed');
    end
    
    project.SaveAs([touchstonefilepath, '.cst'], 1);

    touchstone = project.TOUCHSTONE();
    touchstone.Reset();
    touchstone.Impedance(zfeed);
    touchstone.Renormalize(1);
    touchstone.FileName(touchstonefilepath);
    touchstone.Write();
else
    dispex('CST result already exists.\n');
end

%%
[parameters, SCST] = CST.LoadData(sprintf('%s.s%ip', touchstonefilepath, Ny));
fsCST = parameters.frequencies;

SaCST = zeros(Ny, length(fsCST));
for(fi = 1:length(fsCST))
    f = fsCST(fi);
    [~, ~, ky0, ~] = k(f, 1, th, ph);
    
    ay = ones(1,Ny) .* exp(-1j .* ky0 .* (1:Ny) .* dy);

    for(ny = 1:Ny)
        SaCST(ny, fi) = sum(ay.' .* SCST(:, ny, fi) ./ ay(ny), 1);
    end
end

ZasCST = squeeze((1+SaCST)./(1-SaCST).*zfeed);

Zas = array.GetInputImpedance(fs, th+eps, ph);

%% Plot CST and Matlab results
ymax = ceil(max([real(Zas(:)); real(ZasCST(:)); imag(Zas(:)); imag(ZasCST(:))])/100)*100;
ymin = floor(min([real(Zas(:)); real(ZasCST(:)); imag(Zas(:)); imag(ZasCST(:))])/100)*100;
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

    plot(hAx, fsCST./1e9, real(ZasCST(ny, :)), 'r');
    addlegendentry(hAx, 'CST');
    plot(hAx, fsCST./1e9, imag(ZasCST(ny, :)), 'r--');
    ylim(hAx, [ymin ymax]);
    ylim(hAx, [-100 200]);
    xlim(hAx, [max([min(fsCST), min(fs)]), min([max(fsCST), max(fs)])]./1e9);
end

return;

%% Overlay infinite array solution
infarray = InfiniteArray(slot, tlineup, tlinedown);
Zasinf = infarray.GetInputImpedance(fs, th+eps, ph);
for(ny = 1:Ny)
    if(mod(ny, 5) == 1)
        [hFig, hAx] = figureex(ceil(ny/5));
        alignplot(hFig, 7, 1, hFig.Number, [], 1);
    end
    hAx = subplot(min(Ny, 5), 1, mod(ny-1, 5)+1);
    plot(hAx, fs./1e9, real(Zasinf), 'b');
    addlegendentry(hAx, 'Infinite');
    plot(hAx, fs./1e9, imag(Zasinf), 'b--');
    drawnow;
    movelegend(hAx, 'se');
end

%% Load Touchstone
[parameters, YCST] = CST.LoadData(sprintf('%s.y%ip', touchstonefilepath, Ny));
fsCST = parameters.frequencies;
figureex;
for(nyp = 1:2)
    hAx = subplot(2, 1, nyp);
    hold(hAx, 'on');
    grid(hAx, 'on');
    box(hAx, 'on');
    repeatcolormap(hAx, 2);
    title(hAx, sprintf('Y_{1,%i}', nyp));
    plot(hAx, fsCST./1e9, real(squeeze(YCST(1,nyp,:))));
    addlegendentry(hAx, 'CST');
    plot(hAx, fsCST./1e9, imag(squeeze(YCST(1,nyp,:))), '--');

    plot(hAx, fs./1e9, real(squeeze(array.Ymat(1,nyp,:))));
    addlegendentry(hAx, 'Matlab');
    plot(hAx, fs./1e9, imag(squeeze(array.Ymat(1,nyp,:))), '--');
    xlim(hAx, [max([min(fsCST), min(fs)]), min([max(fsCST), max(fs)])]./1e9);
end

% %% Set matlab Y matrix to CST
% [parameters, YCST] = CST.LoadData(sprintf('%s_aligned.y%ip', touchstonefilepath, Ny));
% fsCST = parameters.frequencies;
% for(fi = 1:length(fs))
%     f = fs(fi);
%     ficst = find(fsCST == f);
%     
%     array.Ymat(:,:, fi) = YCST(:,:,ficst);
% end
% 
% Zas = array.GetInputImpedance(fs, th, ph);

%% Alternatively use txt
[parameters, data] = CST.LoadData(sprintf('%s_Y.txt', touchstonefilepath));
fsCST = parameters{1}.frequencies.*1e9;
YCST = zeros(Ny, Ny, length(fsCST));
for(ny = 1:Ny)
    for(nyp = 1:Ny)
        YCST(ny, nyp, :) = data{ny + (nyp-1)*Ny};
    end
end
for(nyp = 1:Ny)
    if(mod(nyp, 5) == 1)
        [hFig, hAx] = figureex;
            repeatcolormap(hAx, 2);
%             alignplot(hFig, 7, 1, hFig.Number, [], 1);
    end
    hAx = subplot(min(Ny, 5), 1, mod(nyp-1, 5)+1);
    hold(hAx, 'on');
    grid(hAx, 'on');
    box(hAx, 'on');
    repeatcolormap(hAx, 2);
    title(hAx, sprintf('Y_{1,%i}', nyp));
    plot(hAx, fsCST./1e9, real(squeeze(YCST(1,nyp,:))));
    addlegendentry(hAx, 'CST');
    plot(hAx, fsCST./1e9, imag(squeeze(YCST(1,nyp,:))), '--');

    plot(hAx, fs./1e9, real(squeeze(array.Ymat(1,nyp,:))));
    addlegendentry(hAx, 'Matlab');
    plot(hAx, fs./1e9, imag(squeeze(array.Ymat(1,nyp,:))), '--');
    xlim(hAx, [max([min(fsCST), min(fs)]), min([max(fsCST), max(fs)])]./1e9);
end
