function Zas = GetInputImpedance(this, fs, excitation)
    if(size(excitation, 1) ~= this.Nx || size(excitation, 2) ~= this.Ny)
        error('Invalid excitation matrix supplied. Should be Nx by Ny.');
    end
    % Ensure the appropriate Z matrices have been calculated.
    this.InitializeZMatrix(fs);
    
    dispex('Active Z: Calculating for %i frequencies, %ix%i elements.\n', length(fs), this.Nx, this.Ny);
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
    dispex('Active Z: Completed in %s.\n', fancyduration(dt));
end