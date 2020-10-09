function [M] = VoltageSpectrum(this, fs, excitation, kx, ky)
    this.InitializeDs(fs);
    this.InitializeZMatrix(fs);

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

    % Do not deform for now.
    deformedpath = 0;

    % Add zeroes to the excitation matrix for the terminations
    excitation_ = zeros(Nx_+2, Ny_);
    excitation_(2:end-1, :) = excitation;

    basispositions = [-dedge_, (0:Nx_-1)*dx, (Nx_-1)*dx+dedge_];
    x = linspace(min(basispositions), max(basispositions), 20);


    M = zeros(length(fs), length(kx));
    for(fi = 1:length(fs)) % PARFOR?
        f = fs(fi);
        lambda = Constants.c0/f;

        g = 5/3 * sqrt(wslot * lambda);

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
        % Reshape to Nx by Ny
        i = reshape(i, Nx_+2, Ny_);

        % Determine the D values for the given kxs.
        fii = find(this.D_fs == f, 1);
        % If the path is straight, there's no need to unfold.
        interpolationkx = real(kx);
        % Interpolate the D vectors at the desired kx values.
        D = zeros(length(kx), Ny_);
        for(nypp = 0:Ny_-1)
            D(:, nypp+1) = this.D_interpolants{deformedpath+1, fii, nypp+1}(interpolationkx);
    %         D(:, nypp+1) = interp1(this.D_kxs{typei, fii, nypp+1}, this.Ds{typei, fii, nypp+1}, real(kx));
        end

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

        sumVJ = 0;
        for(ny = 0:Ny_-1)
            V_kx_ny = zeros(size(kx));
            for(kxi = 1:length(kx)) % PARFOR?
                if(Ny_ == 1) % Indexing D(kxi, :) fails when Ny_ is 1.
                    Dparam = D(kxi);
                else
                    Dparam = D(kxi, :);
                end
                V_kx_ny(kxi) = dy ./ (2*pi) .* fastintegral(@(ky) V_Integrand_ky(kx(kxi), ky, Nx_, Ny_, ny, dy, Dparam, Fs, basispositions, i), -pi/dy, pi./dy);
            end

            sumVJ = sumVJ + V_kx_ny .* besselj(0, ky .* wslot ./ 2) .* exp(1j .* ky .* ny .* dy);
        end

        M(fi, :) = sumVJ;
    end
end
function v = V_Integrand_ky(kx, ky, Nx_, Ny_, ny, dy, D, Fs, xs, i)
    nypp = (0:Ny_-1).';
    sumD = sum(D.' .* exp(1j .* (nypp * ky)  .* dy), 1);
    
%     sumD2 = 0;
%     for(nypp = 0:Ny_-1)
%         sumD2 = sumD2 + D(nypp+1) .* exp(1j .* (nypp * ky)  .* dy);
%     end
    
%     sumIF = 0;
%     for(nyp = (0:Ny_-1))
%         for(nxp = (-1:Nx_))
% %             sumIF = sumIF + i(nxp+2, nyp+1) .* Fs{nxp+2}(kx) .* exp(1j .* kx .* xs(nxp+2)) ...
% %                                                              .* exp(-1j .* ky .* abs(ny-nyp) .* dy);
%             sumIF = sumIF + i(nxp+2, nyp+1) .* exp(1j .* kx .* xs(nxp+2)) ...
%                                                              .* exp(-1j .* ky .* abs(ny-nyp) .* dy);
%         end
%     end
    sumIF = 0;
    nyp = (0:Ny_-1).';
    for(nxp = (-1:Nx_))
        sumIF = sumIF + sum(i(nxp+2, :).' .* exp(1j .* kx .* xs(nxp+2)) .* Fs{nxp+2}(kx) .* exp(-1j .* (abs(ny-nyp) * ky) .* dy), 1);
    end
    
    v = -1./sumD .* sumIF;% .* exp(-1j .* ky .* dy .* ny); % Exponent moved into absolute above
end
