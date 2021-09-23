% close all;
closewaitbars;
clearvars -except array*;
SetupPath;
clear global;

f0 = 31e9;
c0 = Constants.c0;
l0 = c0/f0;

fs = (10:2.5:35) * 1e9;

th = eps * pi/180;
ph = 0 * pi/180;

dx = 0.45*l0;
dy = dx;
wslot = 0.05*l0;
dslot = 0.05*l0;
walled = 1;
dedge = 0.25*l0;
zfeed = 100;

% Number of unit cells.
Nx = 3;
Ny = 3;
excitation = ones(Nx, Ny); 

tlineup = FreeSpace();
tlinedown = ShortedLine(1, 0.25*l0);

slot = Slot(dx, dy, wslot, dslot, walled);

if(~exist('array', 'var') || ~contains(class(array), 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny || array.unitcell.walled ~= walled)
    array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
%     array = FiniteArray_NoTermination(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
end

%% Determine path to place the CST files.
touchstonepath = sprintf('%s/Validations/%s/', resultdir_matlab, mfilename);
path = sprintf('%s/Validations/%s/', resultdir_matlab, mfilename);

filename = sprintf('%ix%i', Nx, Ny);
if(walled)
    filename = [filename, '_walled'];
end
cstfilepath = [path, filename];
touchstonefilepath = [touchstonepath, filename];
% touchstonefilepath = [touchstonepath, filename, '_notermination'];

matpath = sprintf('%s/Validations/%s/%ix%i', resultdir_matlab, mfilename, Nx, Ny);
if(~exist(sprintf('%s.mat', matpath), 'file'))
    save(sprintf('%s.mat', matpath), 'array', '-v7');
end

if(~exist(sprintf('%s.s%ip', touchstonefilepath, Nx*Ny), 'file'))
    dispex('Running CST simulation for %ix%i.\n', Nx, Ny);
    tc = tic;
    project = CST.InitializeBasicProject();

    project.StoreParameter('fmin', min(fs)/1e9);
    project.StoreParameter('fmax', max(fs)/1e9);
    project.StoreParameter('fmesh', 'fmax');
    project.StoreParameter('slot_impedance', zfeed);
    project.StoreParameter('nsamplesperGHz', 1e9/(fs(2) - fs(1)));

    boundary = project.Boundary();
    boundary.Zmin('expanded open');
    boundary.Zmax('expanded open');

    array.BuildCST(project);

    project.Rebuild();

    % Calculate in CST
    fdsolver = project.FDSolver();
    if(~fdsolver.Start()); dispex('CST simulation failed.'); return; end

    project.SaveAs([cstfilepath, '.cst'], 1);
    dt_CST = toc(tc);
    dispex('CST took %.1fs for %ix%i elements.\n', dt_CST, Nx, Ny);

    touchstone = project.TOUCHSTONE();
    touchstone.Reset();
    touchstone.Impedance(zfeed);
    touchstone.Renormalize(1);
    touchstone.FileName(touchstonefilepath);
    touchstone.Write();

    project.Quit();
else
    dispex('CST simulation results already exist for %ix%i.\n', Nx, Ny);
end

% Calculate in MATLAB
matContents = who('-file', sprintf('%s.mat', matpath));
if(ismember('array', matContents))
    arrayM = load(sprintf('%s.mat', matpath), 'array');
    if(isa(array, class(arrayM)))
        array = arrayM.array;
    end
end
if(ismember('dt_ML', matContents))
    dt_ML_oldM = load(sprintf('%s.mat', matpath), 'dt_ML');
    dt_ML_old = dt_ML_oldM.dt_ML;
else
    dt_ML_old = [];
end

tc = tic;
% array.InitializeDs(fs);
% dispex('Saving array variable.\n');
% save(matpath, 'array', '-append');
% 
% if(ismethod(array, 'InitializeKyInts'))
%     array.InitializeKyInts(fs);
%     dispex('Saving array variable.\n');
%     save(matpath, 'array', '-append');
% end
% 
array.InitializeZMatrix(fs);
% dispex('Saving array variable.\n');
% save(matpath, 'array', '-append');

Zas = array.GetInputImpedance(fs, excitation);

dt_ML = toc(tc);
dispex('MATLAB took %.1fs for %ix%i elements.\n', dt_ML, Nx, Ny);

[parameters, S] = CST.LoadData(sprintf('%s.s%ip', touchstonefilepath, Nx*Ny));
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

% colors = lines(Nx*Ny);
% [hFig, hAx] = figureex;
% for(nx = 1:Nx)
%     for(ny = 1:Ny)
%         clr = colors((nx-1)*Ny+ny, :);
%         plot(hAx, fsCST/1e9, squeeze(real(Zact(nx, ny, :))), 'Color', clr);
%         plot(hAx, fsCST/1e9, squeeze(imag(Zact(nx, ny, :))), '--', 'Color', clr);
%     end
% end
% 
% for(nx = 1:Nx)
%     for(ny = 1:Ny)
%         clr = colors((nx-1)*Ny+ny, :);
%         plot(hAx, fs/1e9, squeeze(real(Zas(nx, ny, :))), 'x', 'Color', clr);
%         plot(hAx, fs/1e9, squeeze(imag(Zas(nx, ny, :))), 'o', 'Color', clr);
%     end
% end

hFig = figureex;
    alignplot(hFig, 4, 4, hFig.Number, [], 1);
for(nx = 1:Nx)
    for(ny = 1:Ny)
        hAx = subplot(Nx, Ny, ny+(nx-1)*Ny);
        hold(hAx, 'on');
        grid(hAx, 'on');
        box(hAx, 'on');
        plot(hAx, fsCST/1e9, squeeze(real(Zact(nx, ny, :))), 'k');
        plot(hAx, fsCST/1e9, squeeze(imag(Zact(nx, ny, :))), 'k--');
        plot(hAx, fs/1e9, squeeze(real(Zas(nx, ny, :))), 'r');
        plot(hAx, fs/1e9, squeeze(imag(Zas(nx, ny, :))), 'r--');
    end
end
legend(hAx, {'reCST', 'imCST', 'reMATLAB', 'imMATLAB'});

savefig(hFig, matpath);
        
% S = (Zas - zfeed) ./ (Zas + zfeed);
% hFig = figureex;
%     alignplot(hFig, 4, 4, hFig.Number, [], 1);
% for(nx = 1:Nx)
%     for(ny = 1:Ny)
%         hAx = subplot(Nx, Ny, ny+(nx-1)*Ny);
%         hold(hAx, 'on');
%         grid(hAx, 'on');
%         box(hAx, 'on');
%         plot(hAx, fsCST/1e9, squeeze(20*log10(abs(Sact(nx, ny, :)))), 'k');
%         plot(hAx, fs/1e9, squeeze(20*log10(abs(S(nx, ny, :)))), 'r');
%         ylim([-30 0]);
%     end
% end
% legend(hAx, {'CST', 'MATLAB'});

savefig(hFig, [matpath, '_S']);
return;
%% Z-matrix
[parameters, ZCST] = CST.LoadData(sprintf('%s.z%ip', touchstonefilepath, Nx*Ny));
fsCST = parameters.frequencies;
Z = array.ReducedZMatrix(fs);

[hFig, hAx] = figureex;
alignplot(hFig, 4, 4, hFig.Number, [], 1);
nx = 0;
ny = 0;
for(nxp = 0:Nx-1)
    for(nyp = 0:Ny-1)
        hAx = subplot(Nx, Ny, (nxp+1)+Nx*nyp);
            title(hAx, sprintf('%i, %i', nxp, nyp));
            hold(hAx, 'on');
            grid(hAx, 'on');
            box(hAx, 'on');
            repeatcolormap(hAx, 2);
            
        i1 = (nx +1) + Nx*ny;
        i2 = (nxp+1) + Nx*nyp;
        Znxny = squeeze(Z(i1, i2, :)).*2;
        i1 = portindexing(nx+1, ny+1);
        i2 = portindexing(nxp+1, nyp+1);
        ZCSTnxny = squeeze(ZCST(i1,i2,:));
        
        plot(hAx, fs./1e9, real(Znxny./2));
        plot(hAx, fs./1e9, imag(Znxny./2), '--');
        
        plot(hAx, fsCST./1e9, real(ZCSTnxny));
        plot(hAx, fsCST./1e9, imag(ZCSTnxny), '--');
    end
end
legend(hAx, {'reMATLAB', 'imMATLAB', 'reCST', 'imCST'});

