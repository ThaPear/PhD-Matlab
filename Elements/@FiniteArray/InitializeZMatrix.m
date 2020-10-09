% Matrix shape, Znx,nxp,ny,nyp
% Z-1,-1,0,0     Z 0,-1,0,0 ...     |   Z-1,-1,0,1     Z 0,-1,0,1
% Z-1, 0,0,0     Z 0, 0,0,0 ...     |   Z-1, 0,0,1     Z 0, 0,0,1
% Z-1, 1,0,0     Z 0, 1,0,0 ...     |
%   :               :               |
% ----------------------------------------------------------------
% Z-1,-1,1,0     Z 0,-1,1,0 ...     |   Z-1,-1,1,1     Z 0,-1,1,1
% Z-1, 0,1,0     Z 0, 0,1,0 ...     |   Z-1, 0,1,1     Z 0, 0,1,1

% Differently formatted
% Within block nx and nxp change (nx,nxp below)
% -1,-1    0,-1    1,-1
% -1, 0    0, 0    1, 0
% -1, 1    0, 1    1, 1
% Between blocks ny and nyp change (ny,nyp below)
% -1,-1 |  0,-1 |  1,-1
% ---------------------
% -1, 0 |  0, 0 |  1, 0
% ---------------------
% -1, 1 |  0, 1 |  1, 1

% To get element Znx,nxp,ny,nyp you use:
% i1 = (nx +1) + (Nx+2)*ny;
% i2 = (nxp+1) + (Nx+2)*nyp;
% array.Zmat(i1, i2)

function InitializeZMatrix(this, fs)
    newfs = setdiff(fs, this.Z_fs);
    if(length(newfs) < 1)
        return;
    end
    % Ensure the appropriate D integrals have been calculated.
    this.InitializeDs(fs);
    % Ensure the appropriate Ky integrals have been calculated.
    this.InitializeKyInts(fs);
    
    dispex('Z matrix: Calculating for %i frequencies, %ix%i elements.\n', length(newfs), this.Nx, this.Ny);
    tc = tic;

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

        parfor(in = 1:N) % parfor
            nx = nxs(in);
            ny = nys(in);
            nxp = nxps(in);
            nyp = nyps(in);
            
            x = xs(nx+2); %#ok<PFBNS> "The entire array or structure 'xs' is a broadcast variable."
            xp = xs(nxp+2);
            F = basisfunctions{nx+2}; %#ok<PFBNS> "The entire array or structure 'basisfunctions' is a broadcast variable."
            Fp = basisfunctions{nxp+2};


%             dispex('nx %i, ny %i, nxp %i, nyp %i\n   x %f, xp %f\n', nx, ny, nxp, nyp, x, xp);

            if(abs(x - xp) < 2*g)
                % For close elements use the straight integration path.
                z = dy./(2*pi).^2.* ...
                    fastintegral(...                                             % deformedpath = 0, straight path
                        @(kx) Z_Integrand_kx(this, f, kx, ny, nyp, x, xp, F, Fp, 0, integrationpath), ...
                        lim1, lim2, 'Waypoints', integrationpath);
            else
                % For far-away elements deform the path downwards.
                z = dy./(2*pi).^2.* ...
                    fastintegral(...                                             % deformedpath = 1, deformed path
                        @(kx) Z_Integrand_kx(this, f, kx, ny, nyp, x, xp, F, Fp, 1, integrationpath_deformed), ...
                        lim1_deformed, lim2_deformed, 'Waypoints', integrationpath_deformed);
                %{
                z2 = dy./(2*pi).^2.* ...
                    fastintegral(...                                                    % deformedpath = 0, straight path
                        @(kx) Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp, 0, integrationpath), ...
                        lim1, lim2, 'Waypoints', integrationpath);
                dispex('%.2f%% error for nx=%i ny=%i nxp=%i nyp=%i.\n', abs(z2-z)/abs(z)*100, nx, ny, nxp, nyp);
                %}
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
    dispex('Z matrix: Completed in %s.\n', fancyduration(dt));
end
function v = Z_Integrand_kx(this, f, kx, ny, nyp, x, xp, F, Fp, deformedpath, integrationpath)
    % Retrieve the correct D index for this frequency.
    fii = find(this.D_fs == f, 1);
    
    if(~deformedpath)
        % If the path is straight, there's no need to unfold.
        interpolationkx = real(kx);
    else
        % If the path is deformed, unfold it to the real axis before interpolation.
        interpolationkx = real(FiniteArray.UnfoldKVector(kx, integrationpath));
    end
    
    % Interpolate the ky integral at the desired kx values.
    dny = abs(ny-nyp);
    int2 = this.KyInt_interpolants{deformedpath+1, fii, dny+1}(interpolationkx);
    
    v = int2 .* Fp(kx) .* F(-kx) .* exp(-1j .* kx .* abs(x-xp));
    
    %{
        [~, hAx] = figureex(200);
            plot(hAx, real(kx), imag(kx), '.');
        k0 = 2*pi*f/3e8;
        [~, hAx] = figureex(201);
            pkx = FiniteArray.UnfoldKVector(kx, integrationpath)/k0;
            plot(hAx, real(pkx), real(int2), 'k.');
            title(hAx, sprintf('ny %i, nyp %i', ny, nyp));
            plot(hAx, real(pkx), imag(int2), 'r.');
    %}

    %{
%         global kxs;
%         global vs;
%         if(nx == 2 && ny == 1 && deformedpath)
%             kxs = [kxs, kx];
%             vs = [vs, int2];
%         end

        global kxs;
        global vs;
        if(nx == 2)
            kxs{dny+1, deformedpath+1} = [kxs{dny+1, deformedpath+1}, kx];
            vs{dny+1, deformedpath+1} = [vs{dny+1, deformedpath+1}, int2];
        end
    %}
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