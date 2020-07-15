classdef FiniteArray < handle
    properties
        unitcell  % Unit cell in the finite array
        
        Nx % Number of unit cells in x
        Ny % Number of unit cells in y
        dedge % Distance to the slot termination
        zfeed % Impedance of the feeding network
        
        Zmat % Mutual impedance matrix
        
        Ds % Precomputed Ds, contains 1/sum(D(nypp,kx))
        D_kxs % Kx values for precomputed Ds
        D_fs  % Frequencies for precomputed Ds
    end
    methods
        function this = FiniteArray(unitcell, Nx, Ny, dedge, zfeed)
            this.unitcell = unitcell;
            this.Nx = Nx;
            this.Ny = Ny;
            this.dedge = dedge;
            this.zfeed = zfeed;
            
        end
        function InitializeZMatrix(this, fs)
            if(isempty(this.Ds) || ~isempty(setdiff(fs, this.D_fs)))
                this.PrecomputeDs(fs);
            end
            tc = tic;
            
            Nf = length(fs);
            Nx_ = this.Nx;
            Ny_ = this.Ny;
            dedge_ = this.dedge;
            
            dx = this.unitcell.dx;
            dy = this.unitcell.dy;
            wslot = this.unitcell.wslot;
            dslot = this.unitcell.dslot;
            tlineup = this.unitcell.tlineup;
            tlinedown = this.unitcell.tlinedown;
            walled = this.unitcell.walled;
            
            z0 = Constants.z0;
            c0 = Constants.c0;
            
%             this = parallel.pool.Constant(this);
                
            %% Initialize progress bar
            hDataQueue = parallel.pool.DataQueue;
            hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
            afterEach(hDataQueue, @updateWaitbar);
            
            Zmat_ = zeros((Nx_+2)*Ny_, (Nx_+2)*Ny_, Nf);
            global fi;
            for(fi = 1:Nf)
                f = fs(fi);
%                 dispex('%i / %i.\n', fi, length(fs));
                
                k0 = 2*pi*f/c0;
                
                lambda = c0/f;
                g = 5/3 * sqrt(wslot * lambda);
                
                Ffeed = @(kx) sinc(kx .* dslot ./ (2*pi));
                Fedge = @(kx) ( besselj(0, kx .* g./2) - 1j .* struve(0, kx.*g./2) ...
                                - 2/pi .* sinc(kx.*g./(4*pi)).*exp(-1j.*kx.*g./4)) .* exp(1j.*kx.*g./2);
                % The basis function of the second termination is the reverse of that of the first
                % termination.
                Fedge2 = @(kx) Fedge(-kx); 
                
                %% Integration limits & path
                % Go to +-inf
                delta = 0.01.*k0;
                lim1 = -100.*k0-1j.*delta;
                lim2 = -lim1;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta];
                
                %% Set up indexing.
                [nxs, nys, nxps, nyps] = ndgrid(-1:Nx_, 0:Ny_-1, -1:Nx_, 0:Ny_-1);
                nxs = nxs(:);
                nys = nys(:);
                nxps = nxps(:);
                nyps = nyps(:);
                
                N = length(nxs);
                progress = -1; send(hDataQueue, nan);
                
                Z = zeros(1, N);
                
                xs = [-dedge_, (0:Nx_-1)*dx, (Nx_-1)*dx+dedge_];
                Fs = cell(1,Nx_+2);
                Fs{1} = Fedge;
                Fs(2:Nx_+1) = {Ffeed};
                Fs{end} = Fedge2;
                
                parfor(in = 1:N)
                    nx = nxs(in);
                    ny = nys(in);
                    nxp = nxps(in);
                    nyp = nyps(in);
%                     nx = 0;
%                     nxp = 0;
%                     ny = 0;
%                     nyp = 0;
                    
                    x = xs(nx+2); %#ok<PFBNS> "The entire array or structure 'xs' is a broadcast variable."
                    xp = xs(nxp+2);
                    F = Fs{nx+2}; %#ok<PFBNS> "The entire array or structure 'Fs' is a broadcast variable."
                    Fp = Fs{nxp+2};
                    
%                     dispex('nx %i, ny %i, nxp %i, nyp %i\n   x %f, xp %f\n', nx, ny, nxp, nyp, x, xp);
                    
                    z = -dy./(2*pi).* ...
                        integral(...
                            @(kx) Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp), ...
                            lim1, lim2, 'Waypoints', integrationpath);
                    Z(in) = z;
                    
                    send(hDataQueue, nan); % Update progress bar
                end
                % Reshape Z such that the rows are the (nx, ny) combinations, and the columns are
                % the (nxp, nyp) for those indices.
                Zmat_(:,:,fi) = reshape(Z, (Nx_+2)*Ny_, []);
            end
            this.Zmat = Zmat_;
            
            function updateWaitbar(~)
                progress = progress + 1;
                waitbar(((fi-1)*N+progress)/(Nf*N), hWaitbar, ...
                    {sprintf('%.1f%% Calculating Z Matrix...', ((fi-1)*N+progress)/(Nf*N)*100), ...
                     sprintf('%i/%i frequencies done.', fi-1, Nf), ...
                     sprintf('%i/%i impedances done.', progress, N)});
            end
            delete(hWaitbar);
            
            toc(tc);
        end
        function PrecomputeDs(this, fs)
            Nf = length(fs);
            Ny_ = this.Ny;
            
            dy = this.unitcell.dy;
            wslot = this.unitcell.wslot;
            tlineup = this.unitcell.tlineup;
            tlinedown = this.unitcell.tlinedown;
            
            z0 = Constants.z0;
            c0 = Constants.c0;
            
            %% Initialize progress bar
            hDataQueue = parallel.pool.DataQueue;
            hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
            afterEach(hDataQueue, @updateWaitbar);
            
            newfs = setdiff(fs, this.D_fs);
            Nf = length(newfs);
            for(fi = 1:Nf)
                f = newfs(fi);
                k0 = 2*pi*f / c0;
                
                % Go to +-inf
                delta = 0.01.*k0;
                lim1 = -100.*k0-1j.*delta;
                lim2 = -lim1;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta];
                
                kx = [linspace(lim1, integrationpath(1), 2950), ...
                      linspace(integrationpath(1), integrationpath(2), 102), ...
                      linspace(integrationpath(2), lim2, 2950)];
                % Remove the duplicate elements
                kx(2951) = [];
                kx(end-2950) = [];
                
                N = length(kx);
                
                progress = -100; send(hDataQueue, nan);
                
                D = zeros(length(kx), Ny_);
                parfor(kxi = 1:length(kx))
                    D_ = zeros(1, Ny_);
                    for(nypp = 0:Ny_-1)
                        D_(nypp+1) = integral(...
                            @(kyp) D_Integrand_kyp(f, dy, k0, kyp, kx(kxi), tlineup, tlinedown, z0, wslot, nypp), ...
                            lim1, lim2, 'Waypoints', integrationpath);
                    end
                    D(kxi, :) = D_;
                    if(~mod(kxi, 100))
                        send(hDataQueue, nan); % Update progress bar
                    end
                end
                
                fii = length(this.Ds)+1;
                this.Ds{fii} = D;
                this.D_fs(fii) = f;
                this.D_kxs{fii} = real(kx);
            end
            
            function updateWaitbar(~)
                progress = progress + 100;
                waitbar(((fi-1)*N+progress)/(Nf*N), hWaitbar, ...
                    {sprintf('%.1f%% Precomputing D...', ((fi-1)*N+progress)/(Nf*N)*100), ...
                     sprintf('%i/%i frequencies done.', fi-1, Nf), ...
                     sprintf('%i/%i kxs done.', progress, N)});
            end
            delete(hWaitbar);
        end
        function Zas = GetInputImpedance(this, fs, excitation)
            if(size(excitation, 1) ~= this.Nx || size(excitation, 2) ~= this.Ny)
                error('Invalid excitation matrix supplied. Should be Nx by Ny.');
            end
            
            Nf = length(fs);
            Nx_ = this.Nx;
            Ny_ = this.Ny;
            zfeed_ = this.zfeed;
            
            % Add zeroes to the excitation matrix for the terminations
            excitation_ = zeros(Nx_+2, Ny_);
            excitation_(2:end-1, :) = excitation;
            
            Zas = zeros((Nx_+2)*Ny_, length(fs));
            for(fi = 1:length(fs))
                Zmat_ = this.Zmat(:, :, fi);
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
%             Zas = Zas(2:end-1, :, :); % Don't drop yet, for debugging
        end
    end
end

function v = Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp)
    % Retrieve the D matrix for this frequency
    fii = find(this.D_fs == f, 1);
    Ds = this.Ds{fii};
    D_kxs = this.D_kxs{fii};
    
    % Interpolate the D matrix at the desired kx
    D = interp1(D_kxs, Ds, real(kx));
    int2 = zeros(1, length(kx));
    for(kxi = 1:length(kx))
        if(Ny_ == 1) % Indexing D(kxi, :) fails when Ny_ is 1.
            int2(kxi) = integral(@(ky) Z_Integrand_ky(ky, Ny_, ny, nyp, dy, D(kxi)), -pi/dy, pi./dy);
        else
            int2(kxi) = integral(@(ky) Z_Integrand_ky(ky, Ny_, ny, nyp, dy, D(kxi, :)), -pi/dy, pi./dy);
        end
    end
    v = int2 .* Fp(kx) .* F(-kx) .* exp(-1j .* kx .* (x-xp));
    return
    %%
    sum(v) %#ok<UNRCH>
    pltx = kx;
    plty = D./(2*pi);
    [hFig, hAx] = figureex;
        title(hAx, 'D./(2*pi)');
        hAx.ColorOrder = lines(min(size(plty)));
        plot(real(pltx), real(plty));
        plot(real(pltx), imag(plty), '--');
        plot(real(pltx), abs(plty), ':');
        alignplot(hFig, 5, 3, [], hFig.Number, 1);
    pltx = kx;
    plty = int2;
    [hFig, hAx] = figureex;
        title(hAx, 'int2');
        hAx.ColorOrder = lines(min(size(plty)));
        plot(real(pltx), real(plty));
        plot(real(pltx), imag(plty), '--');
        plot(real(pltx), abs(plty), ':');
        alignplot(hFig, 5, 3, [], hFig.Number, 1);
    pltx = kx;
    plty = v;
    [hFig, hAx] = figureex;
        title(hAx, 'v');
        hAx.ColorOrder = lines(min(size(plty)));
        plot(real(pltx), real(plty));
        plot(real(pltx), imag(plty), '--');
        plot(real(pltx), abs(plty), ':');
        alignplot(hFig, 5, 3, [], hFig.Number, 1);
end

function v = Z_Integrand_ky(ky, Ny_, ny, nyp, dy, D)
%     sumval = 0;
%     for(nypp = 0:Ny_-1)
%         sumval = sumval + D(nypp+1) .* exp(1j .* ky .* nypp .* dy);
%     end

    nypp = (0:Ny_-1).';
    sumval = sum(D.' .* exp(1j .* (nypp * ky)  .* dy), 1);
    v = exp(-1j .* ky .* abs(ny-nyp) .* dy) ./ sumval;
    
%     sumval = sum(D.', 1);
%     v = ones(size(ky)) ./ sumval;

    return;
    %%
    sum(v) %#ok<UNRCH>
    pltx = ky;
    plty = D.' .* exp(1j .* (nypp * ky)  .* dy);
    [hFig, hAx] = figureex;
        title(hAx, 'D.'' .* exp(1j .* (nypp * ky)  .* dy)');
        hAx.ColorOrder = lines(min(size(plty)));
        plot(real(pltx), real(plty));
        plot(real(pltx), imag(plty), '--');
        plot(real(pltx), abs(plty), ':');
        alignplot(hFig, 5, 3, [], hFig.Number, 1);
    pltx = ky;
    plty = sumval;
    [hFig, hAx] = figureex;
        title(hAx, 'sumval');
        hAx.ColorOrder = lines(min(size(plty)));
        plot(real(pltx), real(plty));
        plot(real(pltx), imag(plty), '--');
        plot(real(pltx), abs(plty), ':');
        alignplot(hFig, 5, 3, [], hFig.Number, 1);
    pltx = ky;
    plty = v;
    [hFig, hAx] = figureex;
        title(hAx, 'v');
        hAx.ColorOrder = lines(min(size(plty)));
        plot(real(pltx), real(plty));
        plot(real(pltx), imag(plty), '--');
        plot(real(pltx), abs(plty), ':');
        alignplot(hFig, 5, 3, [], hFig.Number, 1);
end

function v = D_Integrand_kyp(f, dy, k0, kyp, kx, tlineup, tlinedown, z0, wslot, nypp)
    Vtm = 1;
    Vte = 1;

    kr = sqrt(kx.^2 + kyp.^2);
    isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kx, kyp, Vtm, Vte, itmup, iteup);
    Gxxup = Ghm.xx;

    isTE = 1;    zdownte = tlinedown.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    zdowntm = tlinedown.GetInputImpedance(isTE, f, k0, kr);

    itedown = 1 ./ zdownte;
    itmdown = 1 ./ zdowntm;

    [Ghm] = SpectralGF.hm(z0, k0, kx, kyp, Vtm, Vte, itmdown, itedown);
    Gxxdown = Ghm.xx;

    Gxx = Gxxup + Gxxdown;

    v = Gxx.*besselj(0, kyp.*wslot./2).*exp(+1j.*kyp.*(-nypp).*dy);
    return
    %%
    sum(v)
    pltx = kyp;
    plty = v;
    [hFig, hAx] = figureex;
        title(hAx, 'v');
        hAx.ColorOrder = lines(min(size(plty)));
        plot(real(pltx), real(plty));
        plot(real(pltx), imag(plty), '--');
        plot(real(pltx), abs(plty), ':');
        alignplot(hFig, 5, 3, [], hFig.Number, 1);
end













