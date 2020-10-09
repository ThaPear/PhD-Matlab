close all;
closewaitbars;
clearvars -except array;
SetupPath;
clear global;

f0 = 31e9;
c0 = Constants.c0;
l0 = c0/f0;
z0 = Constants.z0;

fs = (9:4:33) * 1e9;

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
% C = inf;
% 
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
Nxs = 10;
Nys = 10;
for(Nx = Nxs)
    for(Ny = Nys)

        excitation = ones(Nx, Ny);

        slot = Slot(dx, dy, wslot, dslot, walled);


        if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
            array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
        end

        %% Determine path to place the CST files
        skipCST = 0;
        [~, hostname] = system('hostname'); hostname = strsplit(hostname, '\n');

        touchstonepath = sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\', mfilename);
        if(strcmp(hostname{1}, 'SRV539'))
            path = sprintf('E:\\data\\Sander\\Finite_Array\\Validations\\%s\\', mfilename);
        elseif(strcmp(hostname{1}, 'TUD211735'))
%             skipCST = 1;
            path = sprintf('E:\\ Simulations\\Finite_Array\\Validations\\%s\\', mfilename);
        else
            path = touchstonepath;
        end
        filename = sprintf('%ix%i_esa', Nx, Ny);
        filepath = [path, filename];
        touchstonefilepath = [touchstonepath, filename];
        
        if(~exist(sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\%ix%i.mat', mfilename, Nx, Ny), 'file'))
            save(sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\%ix%i.mat', mfilename, Nx, Ny), 'array', '-v7');
        end

        if(~skipCST && ~exist(sprintf('%s.s%ip', touchstonefilepath, Nx*Ny), 'file'))
            dispex('Running CST simulation for %ix%i.\n', Nx, Ny);
            tc = tic;
            project = CST.InitializeBasicProject();

%             project.StoreParameter('fmin', min(fs)/1e9);
%             project.StoreParameter('fmax', max(fs)/1e9);
            project.StoreParameter('fmesh', 'fmax');
            project.StoreParameter('slot_impedance', zfeed);
            project.StoreParameter('nsamplesperGHz', 1e9/(fs(2) - fs(1)));

            boundary = project.Boundary();
            boundary.Zmin('electric');
            boundary.Zmax('expanded open');

            array.BuildCST(project);

            project.Rebuild();
            project.SaveAs([filepath, '.cst'], 0);

            % Calculate in CST
            fdsolver = project.FDSolver();
            if(~fdsolver.Start()); dispex('CST simulation failed.'); return; end
            dt_CST = toc(tc);
            save(sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\%ix%i.mat', mfilename, Nx, Ny), 'dt_CST', '-v7', '-append');
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
        matpath = sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\%ix%i', mfilename, Nx, Ny);
        savematlabresults = 1;
        matContents = who('-file', sprintf('%s.mat', matpath));
        
        if(ismember('array', matContents))
            arrayM = load(sprintf('%s.mat', matpath), 'array');
            array = arrayM.array;
        end
        if(ismember('dt_ML', matContents))
            dt_ML_oldM = load(sprintf('%s.mat', matpath), 'dt_ML');
            dt_ML_old = dt_ML_oldM.dt_ML;
        else
            dt_ML_old = [];
        end

        fs = (12:2:32)*1e9;
        tc = tic;
        array.InitializeDs(fs);
        dispex('Saving array variable.\n');
        save(matpath, 'array', '-append');
        
        if(ismethod(array, 'InitializeKyInts'))
            array.InitializeKyInts(fs);
            dispex('Saving array variable.\n');
            save(matpath, 'array', '-append');
        end
        
        array.InitializeZMatrix(fs);
        dispex('Saving array variable.\n');
        save(matpath, 'array', '-append');
        
        Zas = array.GetInputImpedance(fs, excitation);
        
        dt_ML = toc(tc);
        dispex('MATLAB took %.1fs for %ix%i elements.\n', dt_ML, Nx, Ny);
%         dt_ML = [dt_ML_old, dt_ML];
%         save(sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\%ix%i', mfilename, Nx, Ny), 'dt_ML', '-v7', '-append');
        

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

%         colors = lines(Nx*Ny);
%         [hFig, hAx] = figureex;
%         for(nx = 1:Nx)
%             for(ny = 1:Ny)
%                 clr = colors((nx-1)*Ny+ny, :);
%                 plot(hAx, fsCST/1e9, squeeze(real(Zact(nx, ny, :))), 'Color', clr);
%                 plot(hAx, fsCST/1e9, squeeze(imag(Zact(nx, ny, :))), '--', 'Color', clr);
%             end
%         end
% 
%         for(nx = 1:Nx)
%             for(ny = 1:Ny)
%                 clr = colors((nx-1)*Ny+ny, :);
%                 plot(hAx, fs/1e9, squeeze(real(Zas(nx, ny, :))), 'x', 'Color', clr);
%                 plot(hAx, fs/1e9, squeeze(imag(Zas(nx, ny, :))), 'o', 'Color', clr);
%             end
%         end

        hFig = figureex;
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
        
        if(savematlabresults)
            savefig(hFig, sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\figures\\%ix%i_Z', mfilename, Nx, Ny));
        end
        
        S = (Zas - zfeed) ./ (Zas + zfeed);
        hFig = figureex;
        for(nx = 1:Nx)
            for(ny = 1:Ny)
                hAx = subplot(Nx, Ny, ny+(nx-1)*Ny);
                hold(hAx, 'on');
                grid(hAx, 'on');
                box(hAx, 'on');
                plot(hAx, fsCST/1e9, squeeze(20*log10(abs(Sact(nx, ny, :)))), 'k');
                plot(hAx, fs/1e9, squeeze(20*log10(abs(S(nx, ny, :)))), 'r');
                ylim([-30 0]);
            end
        end
        legend(hAx, {'CST', 'MATLAB'});
        
        if(savematlabresults)
            savefig(hFig, sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\figures\\%ix%i_S', mfilename, Nx, Ny));
        end
        clear hFig hAx;
    end
end