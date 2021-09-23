function [x, v] = Current(this, fs, excitation, ny)
    this.InitializeDs(fs);
    this.InitializeZMatrix(fs);

    dedge_ = this.dedge;
    Nx_ = this.Nx;
    Ny_ = this.Ny;
    dx = this.unitcell.dx;
    dy = this.unitcell.dy;
    wslot_ = this.unitcell.wslot;

    Nf = length(fs);
    N = 500;
    basispositions = [-dedge_, (0:Nx_-1)*dx, (Nx_-1)*dx+dedge_];
    x = linspace(min(basispositions), max(basispositions), N);
    
    lambda = Constants.c0/fs;
    g = 5/3 * sqrt(wslot_ * lambda);
    x = linspace(min(basispositions)-2*g, max(basispositions)+2*g, N);

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

        %% Integration limits & path
        delta = 0.01.*k0;
        % Straight path, go to +-inf
        lim1 = -100.*k0-1j.*delta;
        lim2 = -lim1;
        % Deform integration path around branch cuts.
        integrationpath = [(-1-1j).*delta, (1+1j).*delta];

        for(xi = 1:length(x)) % PARFOR
            v(xi) = 1/2 .* fastintegral(...
                @(kx) v_Integrand_kx(this, f, kx, x(xi), excitation, ny), ...
                lim1, lim2, 'Waypoints', integrationpath);

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
function v = v_Integrand_kx(this, f, kx, x, excitation, ny)
    Vny = this.CurrentSpectrum(f, excitation, kx, ny);
    
    v =1/(2*pi) .* Vny .* exp(-1j .* kx .* x);
end
