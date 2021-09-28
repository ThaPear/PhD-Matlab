function InitializeZMatrix(this, fs)
    newfs = setdiff(fs, this.Z_fs);
    if(length(newfs) < 1)
        return;
    end
    % Ensure D has been calculated.
    this.InitializeDs(newfs);
    
    dispex('Z matrix: Calculating for %i frequencies, %ix%i elements.\n', length(newfs), this.Nx, this.Ny);
    tcWaitbar = tic;

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
    hWaitbar = waitbar(0, {'0% Initializing...', '', '', ''});
    afterEach(hDataQueue, @updateWaitbar);

%     Ymat_ = zeros((Nx_+2)*Ny_, (Nx_+2)*Ny_, Nf);
%     global fi;
    for(fi = 1:Nf)
        f = newfs(fi);
        D_fi = find(this.D_fs == f, 1);
%         dispex('%i / %i.\n', fi, length(fs));

        k0 = 2*pi*f/c0;

        lambda = c0/f;
        g = 5/3 * sqrt(wslot * lambda);

        Ffeed = @(kx) sinc(kx .* dslot ./ (2*pi));
        Fright = @(kx) (besselj(0, kx .* g./2) - 1j .* struve(0, kx.*g./2) ...
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

        %% Set up indexing.
        [nxs, nys, nxps, nyps] = Z_Indexing(Nx_, Ny_);

        N = length(nxs); 
        progress = -1; send(hDataQueue, nan);

        basispositions = [-dedge_, (0:Nx_-1)*dx, (Nx_-1)*dx+dedge_];
        basisfunctions = cell(1,Nx_+2);
        basisfunctions{1} = Fleft;
        basisfunctions(2:Nx_+1) = {Ffeed};
        basisfunctions{end} = Fright;
        
        deformedpath = 0;
        
        %% Prepare the variables for slicing in the parfor.
        xs = zeros(1,N);
        xps = zeros(1,N);
        Fs = cell(1,N);
        Fps = cell(1,N);
        Dinv_interpolants = cell(1,N);
        for(ni = 1:N)
            nx  = nxs(ni);
            ny  = nys(ni);
            nxp = nxps(ni);
            nyp = nyps(ni);

            xs(ni)  = basispositions(nx +2);
            xps(ni) = basispositions(nxp+2);
            Fs{ni}  = basisfunctions{nx +2};
            Fps{ni} = basisfunctions{nxp+2};

            % Broadcasting 'this' to the parfor below caused frequency points to go from 8 minutes
            % to over 2 hours.
            Dinv_interpolants{ni} = this.Dinv_interpolants{deformedpath+1, D_fi, ny+1, nyp+1};
        end
        %% Run the integrals to fill the Z-matrix.
        Zvec = zeros(1, N);
        parfor(ni = 1:N) % parfor
            x  = xs(ni);
            xp = xps(ni);
            F  = Fs{ni};
            Fp = Fps{ni};
            Dinv_interpolant = Dinv_interpolants{ni};
            
            Zvec(ni) = -1/(2*pi) .* integral(@(kx) Z_integrand_kx(kx, F, Fp, x, xp, Dinv_interpolant, deformedpath, integrationpath), ...
                                         lim1, lim2, 'Waypoints', integrationpath);

            send(hDataQueue, nan); % Update progress bar
        end
        % Reshape the Y vector into the correct matrix.
        this.Zmat{end+1} = Z_Matrix(Zvec, Nx_, Ny_);
        this.Z_fs(end+1) = f;
    end
    function updateWaitbar(~)
        progress = progress + 1;
        
        % Estimate ETA.
        dt = toc(tcWaitbar);
        fractiondone = ((fi-1)*N+progress)/(Nf*N);
        eta = dt / fractiondone * (1-fractiondone);
        if(isnan(eta) || isinf(eta))
            eta = 0;
        end
        
        waitbar(fractiondone, hWaitbar, ...
            {sprintf('%.1f%% Calculating Z Matrix...', fractiondone*100), ...
             sprintf('%i/%i frequencies done.', fi-1, Nf), ...
             sprintf('%i/%i impedances done.', progress, N), ...
             sprintf('%s (~%s remaining)', fancyduration(dt), fancyduration(eta))});
    end
    delete(hWaitbar);

    dt = toc(tcWaitbar);
    dispex('Z matrix: Completed in %s.\n', fancyduration(dt));
end
function integrand = Z_integrand_kx(kx, F, Fp, x, xp, Dinv_interpolant, deformedpath, integrationpath)%, f, Ny, ny, nyp, dy, wslot)

    if(deformedpath)
        interpolationkx = real(FiniteArray.UnfoldKVector(kx, integrationpath));
    else
        interpolationkx = real(kx);
    end
    Dinv_nynyp = Dinv_interpolant(interpolationkx);

    %{
    % Compare interpolated version to actual Dinv of free-space
    c0 = 3e8;
    k0 = 2*pi*f / c0;
    z0 = 120*pi;
    Dinvnynyp = zeros(1, length(kx));
    for(ki = 1:length(kx))
        K = -1i*sqrt(-(k0.^2-kx(ki).^2));
        
        % Free-space up & down (hence *2)
        Dself = 2*-0.5/k0/z0*(k0.^2-kx(ki).^2).*besselj(0,wslot/4*K).*besselh(0,2,wslot/4*K);
        Dmutual = @(dny) 2*-0.5/k0/z0*(k0.^2-kx(ki).^2).*besselh(0,2,K*dy*abs(dny));
        
        Dvec = [Dself, Dmutual(1:Ny-1)];
        Dmat = complextoeplitz(Dvec);

        Dinv = pinv(Dmat);
        Dinvnynyp(ki) = Dinv(ny+1, nyp+1);
    end

    [hFig, hAx] = figureex;
        repeatcolormap(hAx, 2);
        plot(hAx, real(kx)./k0, real(Dinvnynyp));
        plot(hAx, real(kx)./k0, imag(Dinvnynyp), '--');
        plot(hAx, real(kx)./k0, real(Dinv_nynyp));
        plot(hAx, real(kx)./k0, imag(Dinv_nynyp), '--');
    [hFig2, hAx2] = figureex;
        figurenextto(hFig2, hFig);
        repeatcolormap(hAx2, 2);
        plot(real(kx)./k0, (real(Dinv_nynyp) - real(Dinvnynyp))./real(Dinvnynyp)*100);
        plot(real(kx)./k0, (imag(Dinv_nynyp) - imag(Dinvnynyp))./imag(Dinvnynyp)*100, '--');
    %}
    
    
%     integrand = Dinv_nynyp .* Fp(kx) .* F(-kx) .* exp(-1j .* kx .* abs(xp-x));
    integrand = Dinv_nynyp .* Fp(kx) .* exp(1j .* kx .* xp) .* F(-kx) .* exp(-1j .* kx .* x);
end
function [nxs, nys, nxps, nyps] = Z_Indexing(Nx, Ny)
    % TODO: Preallocate?
    nxs = [];
    nys = [];
    nxps = [];
    nyps = [];
    for(depth = 0:ceil(Ny/2))
        % Interaction of ny=depth and nyp=depth:Ny-depth
        % Meaning the interaction of the central elements to each other.
            % nx=-1 and nxp=-1:Nx, edge to all
            %                                               nx, ny,        nxp,    nyp
            [nxs_mat, nys_mat, nxps_mat, nyps_mat] = ndgrid(-1, depth, -1:Nx,  depth:Ny-1-depth);
            nxs = [nxs; nxs_mat(:)]; %#ok<AGROW>
            nys = [nys; nys_mat(:)]; %#ok<AGROW>
            nxps = [nxps; nxps_mat(:)]; %#ok<AGROW>
            nyps = [nyps; nyps_mat(:)]; %#ok<AGROW>
            % nx=0 and nxp=0:Nx-1, feed to feeds
            %                                               nx, ny,    nxp,     nyp
            [nxs_mat, nys_mat, nxps_mat, nyps_mat] = ndgrid(0,  depth, 0:Nx-1,  depth:Ny-1-depth);
            nxs = [nxs; nxs_mat(:)]; %#ok<AGROW>
            nys = [nys; nys_mat(:)]; %#ok<AGROW>
            nxps = [nxps; nxps_mat(:)]; %#ok<AGROW>
            nyps = [nyps; nyps_mat(:)]; %#ok<AGROW>
    end
%     nymat = Z_Matrix(nys, Nx, Ny);
%     nxmat = Z_Matrix(nxs, Nx, Ny);
%     nypmat = Z_Matrix(nyps, Nx, Ny);
%     nxpmat = Z_Matrix(nxps, Nx, Ny);
end
function Zmat = Z_Matrix(Zvec, Nx, Ny)
    Zmat = zeros((Nx+2)*Ny, (Nx+2)*Ny);
    
    ind0 = 0;
    for(depth = 0:ceil(Ny/2)-1)
        % Generate the blocks
        blocks = cell(1,Ny-2*depth);
        for(ny = 0:Ny-1-2*depth)
            % First the Zvec has Nx_+2 elements which give the coupling between -1 to nx (edge to all)
            % This repeats Ny times
            ind = ind0 + (1:Nx+2)+ny*(Nx+2);
            block = complextoeplitz(Zvec(ind));
            % Then we have the the coupling between 0 to nx (feed to feed)
            ind = ind0 + (1:Nx) + (Ny-depth*2)*(Nx+2) + ny*Nx;
            block(2:end-1, 2:end-1) = complextoeplitz(Zvec(ind));
            blocks{ny+1} = block;
        end
        ind0 = max(ind);

        Zmat(1+((Nx+2)*depth:((Nx+2)*(Ny-depth)-1)), 1+((Nx+2)*depth:((Nx+2)*(Ny-depth)-1))) = cell2mat(blocks(toeplitz(1:Ny-2*depth)));
    end
end

% function [nxs, nys, nxps, nyps] = Z_Indexing_Old(Nx, Ny)
%     % Interaction of ny=0:Ny_-1 and nyp=0 % All slots
%     % nx=-1:Nx_ and nxp=-1, edge to all
%     %                                               nx,      ny,      nxp, nyp
%     [nxs_mat, nys_mat, nxps_mat, nyps_mat] = ndgrid(-1:Nx,  0:Ny-1, -1,  0);
%     nxs = nxs_mat(:);
%     nys = nys_mat(:);
%     nxps = nxps_mat(:);
%     nyps = nyps_mat(:);
%     % nx=0:Nx_-1 and nxp=0, feed to feeds
%     %                                               nx,      ny,      nxp, nyp
%     [nxs_mat, nys_mat, nxps_mat, nyps_mat] = ndgrid(0:Nx-1, 0:Ny-1, 0,   0);
%     nxs = [nxs; nxs_mat(:)];
%     nys = [nys; nys_mat(:)];
%     nxps = [nxps; nxps_mat(:)];
%     nyps = [nyps; nyps_mat(:)];
% end
% function Zmat = Z_Matrix_Old(Zvec, Nx_, Ny_)
%     % Generate the blocks
%     blocks = cell(1,Ny_);
%     for(ny = 0:Ny_-1)
%         % First the Zvec has Nx_+2 elements which give the coupling between -1 to nx (edge to all)
%         % This repeats Ny_ times
%         ind = (1:Nx_+2)+ny*(Nx_+2);
%         block = toeplitz(real(Zvec(ind)))+1j .* toeplitz(imag(Zvec(ind)));
%         % Then we have the the coupling between 0 to nx (feed to feed)
%         ind = (1:Nx_) + Ny_*(Nx_+2) + ny*Nx_;
%         block(2:end-1, 2:end-1) = toeplitz(real(Zvec(ind)))+1j .* toeplitz(imag(Zvec(ind)));
%         blocks{ny+1} = block;
%     end
%     
%     Zmat = cell2mat(blocks(toeplitz(1:Ny_)));
% end
% function [nxs, nys, nxps, nyps] = Z_Indexing_Full(Nx, Ny)
%     %                                               nx,     ny,     nxp,     nyp
%     [nxs_mat, nys_mat, nxps_mat, nyps_mat] = ndgrid(-1:Nx,  0:Ny-1, -1:Nx,  0:Ny-1);
%     nxs = nxs_mat(:);
%     nys = nys_mat(:);
%     nxps = nxps_mat(:);
%     nyps = nyps_mat(:);
% end
% function Zmat = Z_Matrix_Full(Zvec, Nx, Ny)
%     Zmat = reshape(Zvec, ((Nx+2)*Ny), []);
% end