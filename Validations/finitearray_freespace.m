% close all;
closewaitbars;
clearvars -except array;
SetupPath;
clear global;

f0 = 31e9;
c0 = Constants.c0;
l0 = c0/f0;

dx = 0.45*l0;
dy = dx;
wslot = 0.05*l0;
dslot = 0.05*l0;
walled = 0;
dedge = 0.25*l0;
zfeed = 100;

% Number of unit cells.
Nx = 5;
Ny = 5;

excitation = ones(Nx, Ny);

tlineup = FreeSpace();
tlinedown = FreeSpace();

slot = Slot(dx, dy, wslot, dslot, walled);


fs = (15:5:55) * 1e9;

th = eps * pi/180;
ph = 0 * pi/180;

array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);

project = CST.InitializeBasicProject();

project.StoreParameter('fmin', 15);
project.StoreParameter('fmax', 55);
project.StoreParameter('fmesh', 'fmax');
project.StoreParameter('slot_impedance', zfeed);

boundary = project.Boundary();
boundary.Zmin('expanded open');
boundary.Zmax('expanded open');

array.BuildCST(project);
%% Determine path to place the CST files
[~, hostname] = system('hostname'); hostname = strsplit(hostname, '\n');
if(strcmp(hostname{1}, 'SRV539'))
    path = sprintf('E:\\data\\Sander\\Finite_Array\\Validations\\%s\\', mfilename);
elseif(strcmp(hostname{1}, 'TUD211735'))
    path = sprintf('E:\\ Simulations\\Finite_Array\\Validations\\%s\\', mfilename);
else
    path = sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\', mfilename);
end
filename = sprintf('%ix%i', Nx, Ny);
filepath = [path, filename];

project.SaveAs([filepath, '.cst'], 0);

% Calculate in CST
tc = tic;
fdsolver = project.FDSolver();
if(~fdsolver.Start()); disp('Simulation failed.'); return; end
dt = toc(tc);
dispex('CST took %.1fs for %ix%i elements.\n', dt, Nx, Ny);
% Calculate in MATLAB
Zas = array.GetInputImpedance(fs, excitation);

touchstone = project.TOUCHSTONE();
touchstone.Reset();
touchstone.Impedance(zfeed);
touchstone.Renormalize(1);
touchstone.FileName(filepath);
touchstone.Write();

[parameters, S] = CST.LoadData(sprintf('%s.s%ip', filepath, Nx*Ny));
fsCST = parameters.frequencies;

%% Define indices in the S matrix for the port at position nx, ny.
portindexing = zeros(Nx,Ny);
portindexing(:,1) = 1:Nx;
for(nx = 1:Nx)
    portindexing(nx, 2:end) = Nx + (nx-1)*(Ny-1) + (1:Ny-1);
end

Sact = zeros(Nx, Ny, size(S, 3));
for(nx = 1:Nx)
    for(ny = 1:Ny)
        Sact(nx, ny, :) = sum(S(:, portindexing(nx, ny), :), 1);
    end
end

Zact = (1+Sact)./(1-Sact).*zfeed;

colors = lines(Nx*Ny);
[figZ, axZ] = figureex;
for(nx = 1:Nx)
    for(ny = 1:Ny)
        clr = colors((nx-1)*Ny+ny, :);
        plot(axZ, fsCST/1e9, squeeze(real(Zact(nx, ny, :))), 'Color', clr);
        plot(axZ, fsCST/1e9, squeeze(imag(Zact(nx, ny, :))), '--', 'Color', clr);
    end
end

for(nx = 1:Nx)
    for(ny = 1:Ny)
        clr = colors((nx-1)*Ny+ny, :);
        plot(axZ, fs/1e9, squeeze(real(Zas(nx, ny, :))), 'x', 'Color', clr);
        plot(axZ, fs/1e9, squeeze(imag(Zas(nx, ny, :))), 'o', 'Color', clr);
    end
end
