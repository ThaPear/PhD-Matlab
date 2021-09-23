% close all;
clearvars -except array;
SetupPath;

f0 = 31e9;
l0 = 3e8/f0;
[k0, ~, ~, ~] = k(f0, 1, 0, 0);

dx = 0.45*l0;
dy = dx;
wslot = 0.1*l0;
dslot = 0.1*l0;
dedge = 0.25*l0;

% dx = 4.35483870967742e-3;
% dy = dx;
% wslot = 0.483870967741936e-3;
% dslot = 0.483870967741936e-3;
% dedge = 2.41935483870968e-3;

zfeed = 80;

ths = linspace(-89, 89, 500) * pi/180;
phs = [0 90] * pi/180;

thscan = eps + 60 * pi/180;
phscan = eps + 90 * pi/180;

phs = phscan;
for(ph = phs)
    Nx = 32;
    Ny = 32;

    tlineup = FreeSpace();
    tlinedown = FreeSpace();
    % tlinedown = ShortedLine(1, l0/4);

    unitcell = Slot(dx, dy, wslot, dslot, 0);

    if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
        array = FiniteArray(unitcell, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
    else
        dispex('Taking parameters from array variable.\n');
        dx = array.unitcell.dx;
        dy = array.unitcell.dy;
        wslot = array.unitcell.wslot;
        dslot = array.unitcell.dslot;
    end
%     array = FiniteArray_v5_NoPrecomputeKy(unitcell, tlineup, tlinedown, Nx, Ny, dedge, zfeed);

    array.InitializeDs(f0);
    array.InitializeZMatrix(f0);

    excitationphase = generateexcitationphase(Nx, Ny, f0, dx, dy, thscan, phscan, 'full');
%     excitationphase = HFSSexcitation(Nx, Ny, f, dx, dy, thscan, phscan, '1-to-32-v1-lossy.s33p');
    excitation = ones(Nx, Ny) .* excitationphase;
    
    %% HFSS Excitation
    %{-
%         matfile = load('e:\ HFSS\Feed with phase v2\60 degrees\ Results\s1to1024.mat');
%         matfile = load('e:\ HFSS\Feed with phase v2\60 degrees\ Results\s1to1024-lossy-old.mat');
        matfile = load('e:\ HFSS\Feed with phase v2\60 degrees stripline\ Results\s1to1024-lossy.mat');
%         matfile = load('e:\ HFSS\Feed with phase v2\30 degrees\ Results\s1to1024-lossy.mat');
%         matfile = load('e:\ HFSS\Feed with phase v2\30 degrees stripline\ Results\s1to1024-lossy.mat');
        s1to1024 = matfile.s1to1024;
        fsS1x = matfile.fs;
        S1x = squeeze(s1to1024(1,:, :));
        S1xmat = reshape(S1x(2:end,:), 32, 32, length(fsS1x));
        S1xmat = permute(S1xmat, [2 1 3]);

        iif0 = find(fsS1x == f0);
        excitation = squeeze(S1xmat(:,:,iif0));
        
    %}

    ff = array.NormalizedFarfield(f0, excitation, ths, ph);

    [hFig, hAx] = figureex;
        hAx.ColorOrder = reshape(repmat(lines(7), 1, 3).', [], 3*7).';
        plot(hAx, ths*180/pi, 20*log10(abs(ff)), ':');
        hFig.Name = sprintf('Nx %i, Ny %i, f0 %i, ph %.1f\\circ', Nx, Ny, f0/1e9, ph*180/pi);
        xlabel(hAx, 'Observation angle [degrees]');
        ylabel(hAx, 'Normalized direcivity [dB]');
        ylim(hAx, [-30 0]);
        xlim(hAx, [-inf inf]);
        alignplot(hFig, 8, 4, hFig.Number, [], 1);
    

end

return;

project = CST.InitializeBasicProject();
project.StoreParameter('fmin', min(f)/1e9);
project.StoreParameter('fmax', max(f)/1e9);
project.StoreParameter('fmesh', 'fmax');
project.StoreParameter('slot_impedance', zfeed);
project.StoreParameter('nsamplesperGHz', 1);

boundary = project.Boundary();
boundary.Zmin('expanded open');
boundary.Zmax('expanded open');

meshadaption3d = project.MeshAdaption3D();
meshadaption3d.MinPasses(6);

array.BuildCST(project);