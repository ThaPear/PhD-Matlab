function M = MagneticCurrentSpectrum(this, fs, excitation, kx, ky)
    if(length(kx) ~= length(ky))
        error('Expected same-length kx and ky vectors. Expecting kx, ky pairs as input.\n');
    end
    dispex('Magnetic Current Spectrum: Calculating for %i slots.\n', this.Ny);
    tcWaitbar = tic;
    
    dedge_ = this.dedge;
    Nx_ = this.Nx;
    Ny_ = this.Ny;
    dx = this.unitcell.dx;
    dy = this.unitcell.dy;
    wslot = this.unitcell.wslot;
    dslot = this.unitcell.dslot;
    zfeed_ = this.zfeed;
    
    Nf = length(fs);
    % TODO: Loop frequency
    assert(Nf == 1);
    fi = 1;
    
    %% Initialize progress bar
    hDataQueue = parallel.pool.DataQueue;
    hWaitbar = waitbar(0, {'0% Initializing...', '', '', ''});
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
    
    function updateWaitbar(~)
        progress = progress + 1;
        
        % Estimate ETA.
        dt = toc(tcWaitbar);
        fractiondone = ((fi-1)*Ny_+progress)/(Nf*Ny_);
        eta = dt / fractiondone * (1-fractiondone);
        if(isnan(eta) || isinf(eta))
            eta = 0;
        end
        
        waitbar(fractiondone, hWaitbar, ...
            {sprintf('%.1f%% Calculating Magnetic Current Spectrum...', fractiondone*100), ...
             sprintf('%i/%i frequencies done.', fi-1, Nf), ...
             sprintf('%i/%i slots done.', progress, Ny_), ...
             sprintf('%s (~%s remaining)', fancyduration(dt), fancyduration(eta))});
    end
    delete(hWaitbar);

    dt = toc(tcWaitbar);
    dispex('Magnetic Current Spectrum: Completed in %s.\n', fancyduration(dt));
end