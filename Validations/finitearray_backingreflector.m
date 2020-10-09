% close all;
closewaitbars;
clearvars -except array;
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
walled = 0;
dedge = 0.25*l0;
zfeed = 100;

% Number of unit cells.
Nxs = 1:5;
Nys = 1:5;
for(Nx = Nxs)
    for(Ny = Nys)

        excitation = ones(Nx, Ny);

        tlineup = FreeSpace();
        tlinedown = ShortedLine(1, 0.25*l0);

        slot = Slot(dx, dy, wslot, dslot, walled);


        if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
            array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
        end

        %% Determine path to place the CST files
        [~, hostname] = system('hostname'); hostname = strsplit(hostname, '\n');

        touchstonepath = sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\', mfilename);
        if(strcmp(hostname{1}, 'SRV539'))
            path = sprintf('E:\\data\\Sander\\Finite_Array\\Validations\\%s\\', mfilename);
        elseif(strcmp(hostname{1}, 'TUD211735'))
            path = sprintf('E:\\ Simulations\\Finite_Array\\Validations\\%s\\', mfilename);
        else
            path = touchstonepath;
        end
        filename = sprintf('%ix%i_br', Nx, Ny);
        filepath = [path, filename];
        touchstonefilepath = [touchstonepath, filename];

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
            boundary.Zmin('electric');
            boundary.Zmax('expanded open');

            array.BuildCST(project);

            project.Rebuild();
            project.SaveAs([filepath, '.cst'], 0);

            % Calculate in CST
            fdsolver = project.FDSolver();
            if(~fdsolver.Start()); dispex('CST simulation failed.'); return; end
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
        matpath = sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\%ix%i', mfilename, Nx, Ny);
        savematlabresults = 1;
        if(exist(sprintf('%s.mat', matpath), 'file'))
            load(sprintf('%s.mat', matpath), 'array');
            savematlabresults = 0;
        end
        tc = tic;
        array.InitializeDs(fs);
        array.InitializeZMatrix(fs);
        Zas = array.GetInputImpedance(fs, excitation);
        if(savematlabresults)
            save(sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\%ix%i', mfilename, Nx, Ny), 'array');
        end
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
            savefig(hFig, sprintf('H:\\Git\\PhD-Matlab\\Validations\\%s\\figures\\%ix%i', mfilename, Nx, Ny));
        end
    end
end