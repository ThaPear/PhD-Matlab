classdef FiniteArray < handle
    properties
        unitcell  % Unit cell in the finite array
        
        Nx % Number of unit cells in x
        Ny % Number of unit cells in y
        dedge % Distance to the slot termination
        zfeed % Impedance of the feeding network
        
        Zmat % Mutual impedance matrix
        Z_fs % Frequencies for impedance matrix
        
        Ds    % Precomputed Ds
        D_kxs % Kx values for precomputed Ds
        D_fs  % Frequencies for precomputed Ds
    end
    methods
        function this = FiniteArray(unitcell, Nx, Ny, dedge, zfeed)
            this.unitcell = unitcell;
            this.Nx = Nx;
            this.Ny = Ny;
            this.dedge = dedge;
            this.zfeed = zfeed;
            
        end
        function InitializeZMatrix(this, fs)
            % Ensure the appropriate D integrals have been calculated.
            if(isempty(this.Ds) || ~isempty(setdiff(fs, this.D_fs)))
                this.InitializeDs(fs);
            end
            tc = tic;
            
            newfs = setdiff(fs, this.Z_fs);
            Nf = length(newfs);
            Nx_ = this.Nx;
            Ny_ = this.Ny;
            dedge_ = this.dedge;
            
            dx = this.unitcell.dx;
            dy = this.unitcell.dy;
            wslot = this.unitcell.wslot;
            dslot = this.unitcell.dslot;
            
            c0 = Constants.c0;
                
            %% Initialize progress bar
            hDataQueue = parallel.pool.DataQueue;
            hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
            afterEach(hDataQueue, @updateWaitbar);
            
%             Zmat_ = zeros((Nx_+2)*Ny_, (Nx_+2)*Ny_, Nf);
%             global fi;
            for(fi = 1:Nf)
                f = newfs(fi);
%                 dispex('%i / %i.\n', fi, length(fs));
                
                k0 = 2*pi*f/c0;
                
                lambda = c0/f;
                g = 5/3 * sqrt(wslot * lambda);
                
                Ffeed = @(kx) sinc(kx .* dslot ./ (2*pi));
                Fright = @(kx) ( besselj(0, kx .* g./2) - 1j .* struve(0, kx.*g./2) ...
                                - 2/pi .* sinc(kx.*g./(4*pi)).*exp(-1j.*kx.*g./4)) .* exp(1j.*kx.*g./2);
                % The basis function of the second termination is the reverse of that of the first
                % termination.
                Fleft = @(kx) Fright(-kx); 
                
                
                
                %% Integration limits & path
                % Go to +-inf
                delta = 0.01.*k0;
                lim1 = -100.*k0-1j.*delta;
                lim2 = -lim1;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta];
                
                %% Set up indexing.
                [nxs, nys, nxps, nyps] = Z_Indexing(Nx_, Ny_);
                
                N = length(nxs);
                progress = -1; send(hDataQueue, nan);
                
                Zvec = zeros(1, N);
                
                xs = [-dedge_, (0:Nx_-1)*dx, (Nx_-1)*dx+dedge_];
                Fs = cell(1,Nx_+2);
                Fs{1} = Fleft;
                Fs(2:Nx_+1) = {Ffeed};
                Fs{end} = Fright;
                
                parfor(in = 1:N)
                    nx = nxs(in);
                    ny = nys(in);
                    nxp = nxps(in);
                    nyp = nyps(in);
%                     nx = 0;
%                     nxp = 0;
%                     ny = 0;
%                     nyp = 0;
                    
                    x = xs(nx+2); %#ok<PFBNS> "The entire array or structure 'xs' is a broadcast variable."
                    xp = xs(nxp+2);
                    F = Fs{nx+2}; %#ok<PFBNS> "The entire array or structure 'Fs' is a broadcast variable."
                    Fp = Fs{nxp+2};
                    
                    
%                     dispex('nx %i, ny %i, nxp %i, nyp %i\n   x %f, xp %f\n', nx, ny, nxp, nyp, x, xp);
                    
                    z = -dy./(2*pi).* ...
                        integral(...
                            @(kx) Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp), ...
                            lim1, lim2, 'Waypoints', integrationpath);
                    Zvec(in) = z;
                    
                    send(hDataQueue, nan); % Update progress bar
                end
                % Reshape the Z vector into the correct matrix, exploiting the block-toeplitz
                % structure of the matrix.
                this.Zmat{end+1} = Z_Matrix(Zvec, Nx_, Ny_);
                this.Z_fs(end+1) = f;
            end
            
            function updateWaitbar(~)
                progress = progress + 1;
                waitbar(((fi-1)*N+progress)/(Nf*N), hWaitbar, ...
                    {sprintf('%.1f%% Calculating Z Matrix...', ((fi-1)*N+progress)/(Nf*N)*100), ...
                     sprintf('%i/%i frequencies done.', fi-1, Nf), ...
                     sprintf('%i/%i impedances done.', progress, N)});
            end
            delete(hWaitbar);
            
            dt = toc(tc);
            dispex('Z matrix took %.1fs for %ix%i elements, %i frequencies.\n', dt, Nx_, Ny_, Nf);
        end
        function InitializeDs(this, fs)
            tc = tic;
            Nf = length(fs);
            Ny_ = this.Ny;
            
            dy = this.unitcell.dy;
            wslot = this.unitcell.wslot;
            tlineup = this.unitcell.tlineup;
            tlinedown = this.unitcell.tlinedown;
            
            z0 = Constants.z0;
            c0 = Constants.c0;
            
            %% Initialize progress bar
            hDataQueue = parallel.pool.DataQueue;
            hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
            afterEach(hDataQueue, @updateWaitbar);
            
            newfs = setdiff(fs, this.D_fs);
            Nf = length(newfs);
            for(fi = 1:Nf)
                f = newfs(fi);
                k0 = 2*pi*f / c0;
                
                % Go to +-inf
                delta = 0.01.*k0;
                lim1 = -100.*k0-1j.*delta;
                lim2 = -lim1;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta];
                
                kx = [linspace(lim1, integrationpath(1), 2950), ...
                      linspace(integrationpath(1), integrationpath(2), 102), ...
                      linspace(integrationpath(2), lim2, 2950)];
                % Remove the duplicate elements
                kx(2951) = [];
                kx(end-2950) = [];
                
                N = length(kx);
                
                progress = -100; send(hDataQueue, nan);
                
                D = zeros(length(kx), Ny_);
                parfor(kxi = 1:length(kx))
                    D_ = zeros(1, Ny_);
                    for(nypp = 0:Ny_-1)
                        D_(nypp+1) = integral(...
                            @(kyp) D_Integrand_kyp(f, dy, k0, kyp, kx(kxi), tlineup, tlinedown, z0, wslot, nypp), ...
                            lim1, lim2, 'Waypoints', integrationpath);
                    end
                    D(kxi, :) = D_;
                    if(~mod(kxi, 100))
                        send(hDataQueue, nan); % Update progress bar
                    end
                end
                
                fii = length(this.Ds)+1;
                this.Ds{fii} = D;
                this.D_fs(fii) = f;
                this.D_kxs{fii} = real(kx);
            end
            
            function updateWaitbar(~)
                progress = progress + 100;
                waitbar(((fi-1)*N+progress)/(Nf*N), hWaitbar, ...
                    {sprintf('%.1f%% Precomputing D...', ((fi-1)*N+progress)/(Nf*N)*100), ...
                     sprintf('%i/%i frequencies done.', fi-1, Nf), ...
                     sprintf('%i/%i kxs done.', progress, N)});
            end
            delete(hWaitbar);
            dt = toc(tc);
            dispex('D took %.1fs for %i slots, %i frequencies.\n', dt, Ny_, Nf);
        end
        function Zas = GetInputImpedance(this, fs, excitation)
            if(size(excitation, 1) ~= this.Nx || size(excitation, 2) ~= this.Ny)
                error('Invalid excitation matrix supplied. Should be Nx by Ny.');
            end
            % Ensure the appropriate Z matrices have been calculated.
            if(isempty(this.Zmat) || ~isempty(setdiff(fs, this.Z_fs)))
                this.InitializeZMatrix(fs);
            end
            tc = tic;
            
            Nf = length(fs);
            Nx_ = this.Nx;
            Ny_ = this.Ny;
            zfeed_ = this.zfeed;
            
            % Add zeroes to the excitation matrix for the terminations
            excitation_ = zeros(Nx_+2, Ny_);
            excitation_(2:end-1, :) = excitation;
            
            Zas = zeros((Nx_+2)*Ny_, Nf);
            for(fi = 1:Nf)
                f = fs(fi);
                % Retrieve the Z matrix for this frequency
                fii = find(this.Z_fs == f, 1);
                Zmat_ = this.Zmat{fii};
                % Build a matrix with zfeed on the diagonal, but only for non-termination elements
                Zl = repmat([0, ones(1, Nx_)*zfeed_, 0], 1, Ny_);
                Zlmat = diag(Zl);
                % Add the loads to the mutual impedances
                Zp = Zmat_ + Zlmat;
                % Determine current through the elements
                i = Zp\excitation_(:);
                % Determine voltage on the elements themselves.
                v = Zmat_ * i;
                
                Zas(:, fi) = v./i;
            end
            
            % Reshape the output matrix and drop the edge impedances.
            Zas = reshape(Zas, Nx_+2, Ny_, Nf);
            Zas = Zas(2:end-1, :, :);

            dt = toc(tc);
            dispex('Active Z took %.1fms for %ix%i elements, %i frequencies.\n', dt*1e3, Nx_, Ny_, Nf);
        end
        function BuildCST(this, project, parentcomponent)
            if(nargin < 3 || isempty(parentcomponent))
                parentcomponent = '';
            else
                if(~strcmp(parentcomponent(end), '/'))
                    parentcomponent = [parentcomponent, '/'];
                end
            end
            componentname = [parentcomponent, 'Array'];
            
            this.unitcell.BuildCST(project, [componentname, 'UnitCell']);
            
            project.StoreParameter('lambda_min', 'c0/fmin/1e9');
            project.StoreParameter('Nx', num2str(this.Nx, '%.15g'));
            project.StoreParameter('Ny', num2str(this.Ny, '%.15g'));
            project.StoreParameter('edge_distance', num2str(this.dedge*1e3, '%.15g'));
            project.StoreParameter('edge_length', '5/3 * sqr(slot_width * lambda_min)');
            project.StoreParameter('padding_x', 'max(0, lambda_min/4-edge_length)');
            project.StoreParameter('padding_y', 'max(0, lambda_min/4-(dy/2-slot_width/2))');
            
            transform = project.Transform();
            brick = project.Brick();
            solid = project.Solid();
            
            %% Copy the unit cell in the x-direction
            % Only if there's more than 1 feed in X.
            project.NextCommandConditional('Nx > 1');
                transform.Reset();
                transform.Name([componentname, 'UnitCell']);
                transform.Vector('dx', 0, 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Nx-1');
                transform.Transform('Shape', 'Translate');
            
            %% Create the terminations
            brick.Reset();
            brick.Name('Termination1');
            brick.Component([componentname, 'UnitCell/Slot']);
            brick.Xrange('-edge_distance-edge_length',  '-dx/2');
            brick.Yrange('-dy/2', 'dy/2');
            brick.Zrange('0', '0');
            brick.Material('PEC');
            brick.Create();
            
            brick.Reset();
            brick.Name('Termination1_Slot');
            brick.Component([componentname, 'UnitCell/Slot']);
            brick.Xrange('-edge_distance',  '-dx/2');
            brick.Yrange('-slot_width/2', 'slot_width/2');
            brick.Zrange('0', '0');
            brick.Material('PEC');
            brick.Create();
            
            brick.Reset();
            brick.Name('Termination2');
            brick.Component([componentname, 'UnitCell/Slot']);
            brick.Xrange('(Nx-1)*dx+dx/2',  '(Nx-1)*dx+edge_distance+edge_length');
            brick.Yrange('-dy/2', 'dy/2');
            brick.Zrange('0', '0');
            brick.Material('PEC');
            brick.Create();
            
            brick.Reset();
            brick.Name('Termination2_Slot');
            brick.Component([componentname, 'UnitCell/Slot']);
            brick.Xrange('(Nx-1)*dx+dx/2',  '(Nx-1)*dx+edge_distance');
            brick.Yrange('-slot_width/2', 'slot_width/2');
            brick.Zrange('0', '0');
            brick.Material('PEC');
            brick.Create();
            
            solid.Subtract([componentname, 'UnitCell/Slot:Termination1'], [componentname, 'UnitCell/Slot:Termination1_Slot']);
            solid.Subtract([componentname, 'UnitCell/Slot:Termination2'], [componentname, 'UnitCell/Slot:Termination2_Slot']);
            
            %% Copy the slot in the y-direction
            % Only if there's more than 1 slot in Y.
            project.NextCommandConditional('Ny > 1');
                transform.Reset();
                transform.Name([componentname, 'UnitCell']);
                transform.Vector(0, 'dy', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Ny-1');
                transform.Transform('Shape', 'Translate');
            
            
            %% Copy the port in the x-direction
            % Only if there's more than 1 feed in X.
            project.NextCommandConditional('Nx > 1');
                transform.Reset();
                transform.Name('port1 (SlotFeed)');
                transform.Vector('dx', 0, 0);
                transform.MultipleObjects(1);
                transform.Repetitions('Nx-1');
                transform.Transform('Port', 'Translate');

            %% Copy the port in the x-direction
            % Only if there's more than 1 slot in Y.
            project.NextCommandConditional('Ny > 1');
                % Copy the port of each nx.
                project.NextCommandLoop('nxi', '1', 'Nx')
                    transform.Reset();
                    transform.Name('port" & nxi & " (SlotFeed)');
                    transform.Vector(0, 'dy', 0);
                    transform.MultipleObjects(1);
                    transform.Repetitions('Ny-1');
                    transform.Transform('Port', 'Translate');
            
            %% Add the padding to ensure lambda/4 spacing from the slot to the open boundary
            project.NextCommandConditional('padding_x > 0');
                brick.Reset();
                brick.Name('Padding_x1');
                brick.Component(componentname);
                brick.Xrange('-edge_distance-edge_length-padding_x',  '-edge_distance-edge_length');
                brick.Yrange('-dy/2', '(Ny-1)*dy+dy/2');
                brick.Zrange('0', '0');
                brick.Material('PEC');
                brick.Create();
            project.NextCommandConditional('padding_x > 0');
                brick.Reset();
                brick.Name('Padding_x2');
                brick.Component(componentname);
                brick.Xrange('(Nx-1)*dx+edge_distance+edge_length',  '(Nx-1)*dx+edge_distance+edge_length+padding_x');
                brick.Yrange('-dy/2', '(Ny-1)*dy+dy/2');
                brick.Zrange('0', '0');
                brick.Material('PEC');
                brick.Create();
            project.NextCommandConditional('padding_y > 0');
                brick.Reset();
                brick.Name('Padding_y1');
                brick.Component(componentname);
                brick.Xrange('-edge_distance-edge_length-padding_x', '(Nx-1)*dx+edge_distance+edge_length+padding_x');
                brick.Yrange('-dy/2-padding_y',  '-dy/2');
                brick.Zrange('0', '0');
                brick.Material('PEC');
                brick.Create();
            project.NextCommandConditional('padding_y > 0');
                brick.Reset();
                brick.Name('Padding_y2');
                brick.Component(componentname);
                brick.Xrange('-edge_distance-edge_length-padding_x', '(Nx-1)*dx+edge_distance+edge_length+padding_x');
                brick.Yrange('(Ny-1)*dy+dy/2',  '(Ny-1)*dy+dy/2+padding_y');
                brick.Zrange('0', '0');
                brick.Material('PEC');
                brick.Create();
        end
    end
end

function v = Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp)
    % Retrieve the D matrix for this frequency
    fii = find(this.D_fs == f, 1);
    Ds = this.Ds{fii};
    D_kxs = this.D_kxs{fii};
    
    % Interpolate the D matrix at the desired kx
    D = interp1(D_kxs, Ds, real(kx));
    int2 = zeros(1, length(kx));
    for(kxi = 1:length(kx))
        if(Ny_ == 1) % Indexing D(kxi, :) fails when Ny_ is 1.
            int2(kxi) = integral(@(ky) Z_Integrand_ky(ky, Ny_, ny, nyp, dy, D(kxi)), -pi/dy, pi./dy);
        else
            int2(kxi) = integral(@(ky) Z_Integrand_ky(ky, Ny_, ny, nyp, dy, D(kxi, :)), -pi/dy, pi./dy);
        end
    end
    v = int2 .* Fp(kx) .* F(-kx) .* exp(-1j .* kx .* (x-xp));
end

function v = Z_Integrand_ky(ky, Ny_, ny, nyp, dy, D)
    nypp = (0:Ny_-1).';
    sumval = sum(D.' .* exp(1j .* (nypp * ky)  .* dy), 1);
    v = exp(-1j .* ky .* (ny-nyp) .* dy) ./ sumval;
end

function v = D_Integrand_kyp(f, dy, k0, kyp, kx, tlineup, tlinedown, z0, wslot, nypp)
    Vtm = 1;
    Vte = 1;

    kr = sqrt(kx.^2 + kyp.^2);
    isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kx, kyp, Vtm, Vte, itmup, iteup);
    Gxxup = Ghm.xx;

    isTE = 1;    zdownte = tlinedown.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    zdowntm = tlinedown.GetInputImpedance(isTE, f, k0, kr);

    itedown = 1 ./ zdownte;
    itmdown = 1 ./ zdowntm;

    [Ghm] = SpectralGF.hm(z0, k0, kx, kyp, Vtm, Vte, itmdown, itedown);
    Gxxdown = Ghm.xx;

    Gxx = Gxxup + Gxxdown;

    v = Gxx.*besselj(0, kyp.*wslot./2).*exp(-1j.*kyp.*nypp.*dy);
end

function [nxs, nys, nxps, nyps] = Z_Indexing(Nx_, Ny_)
    % Interaction of ny=0:Ny_-1 and nyp=0 % All slots
        % nx=-1:Nx_ and nxp=-1, edge to all
        %                                               nx,      ny,      nxp, nyp
        [nxs_mat, nys_mat, nxps_mat, nyps_mat] = ndgrid(-1:Nx_,  0:Ny_-1, -1,  0);
        nxs = nxs_mat(:);
        nys = nys_mat(:);
        nxps = nxps_mat(:);
        nyps = nyps_mat(:);
        % nx=0:Nx_-1 and nxp=0, feed to feeds
        %                                               nx,      ny,      nxp, nyp
        [nxs_mat, nys_mat, nxps_mat, nyps_mat] = ndgrid(0:Nx_-1, 0:Ny_-1, 0,   0);
        nxs = [nxs; nxs_mat(:)];
        nys = [nys; nys_mat(:)];
        nxps = [nxps; nxps_mat(:)];
        nyps = [nyps; nyps_mat(:)];
end

function Zmat = Z_Matrix(Zvec, Nx_, Ny_)
    % Generate the blocks
    blocks = cell(1,Ny_);
    for(ny = 0:Ny_-1)
        % First the Zvec has Nx_+2 elements which give the coupling between -1 to nx (edge to all)
        % This repeats Ny_ times
        ind = (1:Nx_+2)+ny*(Nx_+2);
        block = toeplitz(real(Zvec(ind)))+1j .* toeplitz(imag(Zvec(ind)));
        % Then we have the the coupling between 0 to nx (feed to feed)
        ind = (1:Nx_) + Ny_*(Nx_+2) + ny*Nx_;
        block(2:end-1, 2:end-1) = toeplitz(real(Zvec(ind)))+1j .* toeplitz(imag(Zvec(ind)));
        blocks{ny+1} = block;
    end
    
    Zmat = cell2mat(blocks(toeplitz(1:Ny_)));
end



%     pltx = kx;
%     plty = v;
%     [hFig, hAx] = figureex;
%         title(hAx, 'v');
%         hAx.ColorOrder = lines(min(size(plty)));
%         plot(real(pltx), real(plty));
%         plot(real(pltx), imag(plty), '--');
%         plot(real(pltx), abs(plty), ':');
%         alignplot(hFig, 5, 3, [], hFig.Number, 1);