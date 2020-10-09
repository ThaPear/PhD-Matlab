function [x, v] = Voltage(this, fs, excitation, ny)
    this.InitializeDs(fs);
    this.InitializeZMatrix(fs);

    dedge_ = this.dedge;
    Nx_ = this.Nx;
    Ny_ = this.Ny;
    dx = this.unitcell.dx;
    dy = this.unitcell.dy;
    wslot = this.unitcell.wslot;
    dslot = this.unitcell.dslot;
    zfeed_ = this.zfeed;

    % Do not deform for now.
    deformedpath = 0;

    % Add zeroes to the excitation matrix for the terminations
    excitation_ = zeros(Nx_+2, Ny_);
    excitation_(2:end-1, :) = excitation;

    basispositions = [-dedge_, (0:Nx_-1)*dx, (Nx_-1)*dx+dedge_];
%     x = linspace(min(basispositions)-g, max(basispositions)+g, 20);
%     x = linspace(-(max(basispositions)+g), max(basispositions)+g, 20);
    N = 30;
    Nf = length(fs);

    %% Initialize progress bar
    hDataQueue = parallel.pool.DataQueue;
    hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
    afterEach(hDataQueue, @updateWaitbar);
    progress = 0;

    v = zeros(Nf, N);
    for(fi = 1:Nf)
        f = fs(fi);
        lambda = Constants.c0/f;
        k0 = 2*pi/lambda;

        g = 5/3 * sqrt(wslot * lambda);
        x = linspace(min(basispositions), max(basispositions), N);
%         x = linspace(-0.02, 0.02, numXpoints);

        %% Solve the equiv circuit
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

        i = reshape(i, Nx_+2, Ny_);

        %% Basis functions
        Ffeed = @(kx) sinc(kx .* dslot ./ (2*pi));
        Fright = @(kx) ( besselj(0, kx .* g./2) - 1j .* struve(0, kx.*g./2) ...
                        - 2/pi .* sinc(kx.*g./(4*pi)).*exp(-1j.*kx.*g./4)) .* exp(1j.*kx.*g./2);
        % The basis function of the second termination is the reverse of that of the first
        % termination.
        Fleft = @(kx) Fright(-kx); 

        Fs = cell(1,Nx_+2);
        Fs{1} = Fleft;
        Fs(2:Nx_+1) = {Ffeed};
        Fs{end} = Fright;

        %% Integration limits & path
        delta = 0.01.*k0;
        % Straight path, go to +-inf
        lim1 = -100.*k0-1j.*delta;
        lim2 = -lim1;
        % Deform integration path around branch cuts.
        integrationpath = [(-1-1j).*delta, (1+1j).*delta];

        parfor(xi = 1:length(x)) % PARFOR
            v(xi) = 1/2 .* fastintegral(...
                @(kx) v_Integrand_kx(this, f, dy, kx, x(xi), Nx_, Ny_, ny, Fs, basispositions, i, deformedpath, integrationpath), ...
                lim1, lim2, 'Waypoints', integrationpath);

%             v(xi) = fastintegral(@(kx) V_Integrand_ky(kx, 1, Nx_, Ny_, ny, dy, 1, Fs, basispositions, i, x(xi), dslot), lim1, lim2, 'Waypoints', integrationpath);
%             v(xi) = fastintegral(@(kx) sinc(kx .* dslot ./ (2*pi)) .* exp(-1j .* kx .* x(xi)), lim1, lim2, 'Waypoints', integrationpath);

            send(hDataQueue, nan); % Update progress bar
        end
    end

    function updateWaitbar(~)
        progress = progress + 1;
        waitbar(((fi-1)*N+progress)/(Nf*N), hWaitbar, ...
            {sprintf('%.1f%% Calculating Voltage Distribution...', ((fi-1)*N+progress)/(Nf*N)*100), ...
             sprintf('%i/%i frequencies done.', fi-1, Nf), ...
             sprintf('%i/%i positions done.', progress, N)});
    end
    delete(hWaitbar);
end
function v = v_Integrand_kx(this, f, dy, kx, x, Nx_, Ny_, ny, Fs, xs, i, deformedpath, integrationpath)
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
    
    int = zeros(size(kx));
    for(kxi = 1:length(kx))
        if(Ny_ == 1) % Indexing D(kxi, :) fails when Ny_ is 1.
            Dparam = D(kxi);
        else
            Dparam = D(kxi, :);
        end
        int(kxi) = fastintegral(@(ky) V_Integrand_ky(kx(kxi), ky, Nx_, Ny_, ny, dy, Dparam, Fs, xs, i), -pi/dy, pi./dy);
    end
    v = dy/(2*pi) .* int .* exp(-1j .* kx .* x);
end
