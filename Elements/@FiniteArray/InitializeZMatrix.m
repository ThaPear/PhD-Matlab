function InitializeZMatrix(this, fs)
    newfs = setdiff(fs, this.Z_fs);
    if(length(newfs) < 1)
        return;
    end

    % Ensure the appropriate D integrals have been calculated.
    if(isempty(this.Ds) || ~isempty(setdiff(fs, this.D_fs)))
        this.InitializeDs(fs);
    end
    tc = tic;

    dispex('Calculating Z matrix for %i frequencies, %ix%i elements.\n', length(newfs), this.Nx, this.Ny);

    Nf = length(newfs);
    Nx_ = this.Nx;
    Ny_ = this.Ny;
    dedge_ = this.dedge;

    dx = this.unitcell.dx;
    dy = this.unitcell.dy;
    wslot = this.unitcell.wslot;
    dslot = this.unitcell.dslot;
    walled = this.unitcell.walled;
    tlinedown_ = this.tlinedown;

    c0 = Constants.c0;

    %% Initialize progress bar
    hDataQueue = parallel.pool.DataQueue;
    hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
    afterEach(hDataQueue, @updateWaitbar);

%     Zmat_ = zeros((Nx_+2)*Ny_, (Nx_+2)*Ny_, Nf);
%     global fi;
    for(fi = 1:Nf)
        f = newfs(fi);
%         dispex('%i / %i.\n', fi, length(fs));

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
        delta = 0.01.*k0;
        % Straight path, go to +-inf
        lim1 = -100.*k0-1j.*delta;
        lim2 = -lim1;
        % Deform integration path around branch cuts.
        integrationpath = [(-1-1j).*delta, (1+1j).*delta];
        % Deformed path, go to +-1j*inf
        lim1_deformed = -5j.*k0-1.*delta;
        lim2_deformed = -5j.*k0+1.*delta+1.5*k0;
        % Deform integration path around branch cuts.
        integrationpath_deformed = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];

        %% Set up indexing.
        [nxs, nys, nxps, nyps] = Z_Indexing(Nx_, Ny_);

        N = length(nxs);
        progress = -1; send(hDataQueue, nan);

        Zvec = zeros(1, N);

        xs = [-dedge_, (0:Nx_-1)*dx, (Nx_-1)*dx+dedge_];
        basisfunctions = cell(1,Nx_+2);
        basisfunctions{1} = Fleft;
        basisfunctions(2:Nx_+1) = {Ffeed};
        basisfunctions{end} = Fright;

        parfor(in = 1:N) % PARFOR
            nx = nxs(in);
            ny = nys(in);
            nxp = nxps(in);
            nyp = nyps(in);
%             nx = 0;
%             nxp = 0;
%             ny = 0;
%             nyp = 0;

            x = xs(nx+2); %#ok<PFBNS> "The entire array or structure 'xs' is a broadcast variable."
            xp = xs(nxp+2);
            F = basisfunctions{nx+2}; %#ok<PFBNS> "The entire array or structure 'basisfunctions' is a broadcast variable."
            Fp = basisfunctions{nxp+2};


%             dispex('nx %i, ny %i, nxp %i, nyp %i\n   x %f, xp %f\n', nx, ny, nxp, nyp, x, xp);

            if(abs(x - xp) < 2*g)
                % For close elements use the straight integration path.
                z = -dy./(2*pi).* ...
                    fastintegral(...                                                    % deformedpath = 0, straight path
                        @(kx) Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp, 0, integrationpath, tlinedown_, walled, wslot), ...
                        lim1, lim2, 'Waypoints', integrationpath);
            else
                % For far-away elements deform the path downwards.
                z = -dy./(2*pi).* ...
                    fastintegral(...                                                    % deformedpath = 1, deformed path
                        @(kx) Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp, 1, integrationpath_deformed, tlinedown_, walled, wslot), ...
                        lim1_deformed, lim2_deformed, 'Waypoints', integrationpath_deformed);
%                 z2 = -dy./(2*pi).* ...
%                     fastintegral(...                                                    % deformedpath = 0, straight path
%                         @(kx) Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp, 0, integrationpath), ...
%                         lim1, lim2, 'Waypoints', integrationpath);
%                 dispex('%.2f%% error for nx=%i ny=%i nxp=%i nyp=%i.\n', abs(z2-z)/abs(z)*100, nx, ny, nxp, nyp);
%                 
%                 [~, hAx] = figureex(200);
%                 cla(hAx);
%                 [~, hAx] = figureex(201);
%                 cla(hAx);
%                 [~, hAx] = figureex(202);
%                 cla(hAx);
%                 [~, hAx] = figureex(203);
%                 cla(hAx);
            end
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
function v = Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp, deformedpath, integrationpath, tlinedown, walled, wslot)
    % Retrieve the correct D index for this frequency.
    fii = find(this.D_fs == f, 1);
    
    if(~deformedpath)
        % If the path is straight, there's no need to unfold.
        interpolationkx = real(kx);
    else
        % If the path is deformed, unfold it to the real axis before interpolation.
        interpolationkx = real(this.UnfoldKVector(kx, integrationpath));
    end
    
    % Interpolate the D vectors at the desired kx values.
    D = zeros(length(kx), Ny_);
    for(nypp = 0:Ny_-1)
        D(:, nypp+1) = this.D_interpolants{deformedpath+1, fii, nypp+1}(interpolationkx);
%         D(:, nypp+1) = interp1(this.D_kxs{typei, fii, nypp+1}, this.Ds{typei, fii, nypp+1}, real(kx));
    end
    
    
    int2 = zeros(1, length(kx));
    for(kxi = 1:length(kx))
        % If we have vertical walls, calculate the Ddown using Floquet modes.
        if(walled)
            Ddown = 0;
            for(my = -20:20)
                k0 = 2*pi*f/3e8;
                z0 = Constants.z0;
                
                Vtm = 1;
                Vte = 1;
                
                kym = -2*pi*my/dy;
                kr = sqrt(kx(kxi).^2 + kym.^2);
                isTE = 1;    ztedown = tlinedown.GetInputImpedance(isTE, f, k0, kr);
                isTE = 0;    ztmdown = tlinedown.GetInputImpedance(isTE, f, k0, kr);

                itedown = 1 ./ ztedown;
                itmdown = 1 ./ ztmdown;

                [Ghm] = SpectralGF.hm(z0, k0, kx(kxi), kym, Vtm, Vte, itmdown, itedown);
                Ghm_xx_down = Ghm.xx;

                Ddown = Ddown + 2*pi/dy .* Ghm_xx_down .* besselj(0, -2*pi*my/dy*wslot/2);
            end
        else
            Ddown = 0;
        end
        
        if(length(kx) == 1) % Indexing D(kxi, :) fails when Ny_ is 1.
            Dparam = D(kxi);
        else
            Dparam = D(kxi, :);
        end
        int2(kxi) = fastintegral(@(ky) Z_Integrand_ky(ky, Ny_, ny, nyp, dy, Dparam, Ddown), -pi/dy, pi./dy);
    end
    v = int2 .* Fp(kx) .* F(-kx) .* exp(-1j .* kx .* abs(x-xp));
    
%     global doplot;
%     if(doplot)
%         [~, hAx] = figureex(200);
%             plot(hAx, real(kx), imag(kx), '.');
%         k0 = 2*pi*f/3e8;
%         [~, hAx] = figureex(201);
%             pkx = this.UnfoldKVector(kx, integrationpath)/k0;
%             plot(hAx, real(pkx), real(int2), 'k.');
%             title(hAx, sprintf('ny %i, nyp %i', ny, nyp));
%             plot(hAx, real(pkx), imag(int2), 'r.');
%     end
end
function v = Z_Integrand_ky(ky, Ny_, ny, nyp, dy, D, Ddown)
%     sumD = zeros(size(ky));
%     for(nypp = -ny:(Ny_-1-ny))
%         sumD = sumD + D(abs(nypp)+1) .* exp(1j .* (nypp * ky)  .* dy);
%     end

    nypp = (0:Ny_-1).';
    sumD = sum(D.' .* exp(1j .* (nypp * ky)  .* dy), 1);
    
    % NOTE: Ddown is zero when there are no vertical walls since the Down stratification will then
    % be included in the D
    v = exp(-1j .* ky .* abs(ny-nyp) .* dy) ./ (sumD + Ddown);
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
