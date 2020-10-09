function Zred = ReducedZMatrix(this, fs)
    tc = tic;
    Zred = cell(1, length(fs));
    for(fi = 1:length(fs))
        Zmati = this.Zmat{fi};

        % Load the edge terminations with shorts, then define the loaded section of the
        % matrix as the elements with a short, and the nonloaded section as the elements
        % without a short.
        % The reduced Z-matrix can then be calculated with the input impedance equation
        % Zin = Z11 - Z12*Z21/(Z22+ZL)

        % Build logical index vector denoting whether or not it's a nonloaded element
        % E.g. for 4x4: % 0 1 1 0  0 1 1 0  0 1 1 0  0 1 1 0
        iNL = logical(repmat([0, ones(1, this.Nx), 0], 1, this.Ny));

        % Nonloaded - The part of the Z-matrix that does not include the terminations
        ZNL = Zmati(iNL, iNL);

        % Loaded - The inverse of the nonloaded part of the matrix
        ZL = Zmati(~iNL, ~iNL);

        % Nonloaded-loaded
        ZNLL = Zmati(iNL, ~iNL);

        % Loaded-nonloaded
        ZLNL = Zmati(~iNL, iNL);

        % Matrix version of the input impedance calculation
        % Yields the reduced Z-matrix without the terminations
        Zred{fi} = ZNL - ZNLL * pinv(ZL) * ZLNL;
    end

    dt = toc(tc);
    dispex('Reducing Z took %.1fms for %ix%i elements, %i frequencies.\n', dt*1e3, Nx_, Ny_, Nf);
end
