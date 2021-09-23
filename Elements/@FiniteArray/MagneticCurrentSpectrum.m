function M = MagneticCurrentSpectrum(this, fs, excitation, kx, ky)
    dispex('Calculating Magnetic Current Spectrum for %i slots.\n', this.Ny);
    if(length(kx) ~= length(ky))
        error('Expected same-length kx and ky vectors. Expecting kx, ky pairs as input.\n');
    end
    
    dedge_ = this.dedge;
    Nx_ = this.Nx;
    Ny_ = this.Ny;
    dx = this.unitcell.dx;
    dy = this.unitcell.dy;
    wslot = this.unitcell.wslot;
    dslot = this.unitcell.dslot;
    zfeed_ = this.zfeed;
    
    % TODO: Loop frequency
    fi = 1;
    Nf = length(fs);
    
    %% Initialize progress bar
    hDataQueue = parallel.pool.DataQueue;
    hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
    afterEach(hDataQueue, @updateWaitbar);
    progress = 0;
    
    %% Calculate sum
    sumVJ = 0;
    for(ny = 0:Ny_-1) % parfor does NOT work here.
        Vny = this.VoltageSpectrum(fs, excitation, kx, ny);
        sumVJ = sumVJ + Vny .* besselj(0, ky .* wslot ./ 2) .* exp(1j .* ky .* ny .* dy);
        
        send(hDataQueue, nan); % Update progress bar
    end
    M = sumVJ;
    
    %% Close progress bar
    delete(hWaitbar);
    return;
    function updateWaitbar(~)
        progress = progress + 1;
        waitbar(((fi-1)*Ny_+progress)/(Nf*Ny_), hWaitbar, ...
            {sprintf('%.1f%% Calculating Magnetic Current Spectrum...', ((fi-1)*Ny_+progress)/(Nf*Ny_)*100), ...
             sprintf('%i/%i frequencies done.', fi-1, Nf), ...
             sprintf('%i/%i slots done.', progress, Ny_)});
    end
end