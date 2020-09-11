classdef FiniteArray < handle
    properties
        unitcell  % Unit cell in the finite array
        
        tlineup   % Stratification above the array
        tlinedown % Stratification below the array
        
        Nx % Number of unit cells in x
        Ny % Number of unit cells in y
        dedge % Distance to the slot termination
        zfeed % Impedance of the feeding network
        
        Zmat % Mutual impedance matrix
        Z_fs % Frequencies for impedance matrix
        
        Ds    % Precomputed Ds
        D_kxs % Kx values for precomputed Ds
        D_fs  % Frequencies for precomputed Ds
        D_interpolants % Interpolant to index at any given kx
    end
    methods
        function this = FiniteArray(unitcell, tlineup, tlinedown, Nx, Ny, dedge, zfeed)
            this.unitcell = unitcell;
            this.tlineup = tlineup;
            this.tlinedown = tlinedown;
            this.Nx = Nx;
            this.Ny = Ny;
            this.dedge = dedge;
            this.zfeed = zfeed;
        end
        function InitializeZMatrix(this, fs)
            newfs = setdiff(fs, this.Z_fs);
            if(length(newfs) < 1)
                return;
            end
            
            % Ensure the appropriate D integrals have been calculated.
            if(isempty(this.Ds) || ~isempty(setdiff(fs, this.D_fs)))
                this.InitializeDs(fs);
            end
            tc = tic;
            
            dispex('Calculating Z matrix for %i frequencies, %ix%i elements.\n', length(newfs), this.Nx, this.Ny);
            
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
            
%             Zmat_ = zeros((Nx_+2)*Ny_, (Nx_+2)*Ny_, Nf);
%             global fi;
            for(fi = 1:Nf)
                f = newfs(fi);
%                 dispex('%i / %i.\n', fi, length(fs));
                
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
                
                parfor(in = 1:N) % PARFOR
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
                    F = basisfunctions{nx+2}; %#ok<PFBNS> "The entire array or structure 'basisfunctions' is a broadcast variable."
                    Fp = basisfunctions{nxp+2};
                    
                    
%                     dispex('nx %i, ny %i, nxp %i, nyp %i\n   x %f, xp %f\n', nx, ny, nxp, nyp, x, xp);
                    
                    if(abs(x - xp) < 2*g)
                        % For close elements use the straight integration path.
                        z = -dy./(2*pi).* ...
                            fastintegral(...                                                    % deformedpath = 0, straight path
                                @(kx) Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp, 0, integrationpath, tlinedown_, walled, wslot), ...
                                lim1, lim2, 'Waypoints', integrationpath);
                    else
                        % For far-away elements deform the path downwards.
                        z = -dy./(2*pi).* ...
                            fastintegral(...                                                    % deformedpath = 1, deformed path
                                @(kx) Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp, 1, integrationpath_deformed, tlinedown_, walled, wslot), ...
                                lim1_deformed, lim2_deformed, 'Waypoints', integrationpath_deformed);
%                         z2 = -dy./(2*pi).* ...
%                             fastintegral(...                                                    % deformedpath = 0, straight path
%                                 @(kx) Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp, 0, integrationpath), ...
%                                 lim1, lim2, 'Waypoints', integrationpath);
%                         dispex('%.2f%% error for nx=%i ny=%i nxp=%i nyp=%i.\n', abs(z2-z)/abs(z)*100, nx, ny, nxp, nyp);
%                         
%                         [~, hAx] = figureex(200);
%                         cla(hAx);
%                         [~, hAx] = figureex(201);
%                         cla(hAx);
%                         [~, hAx] = figureex(202);
%                         cla(hAx);
%                         [~, hAx] = figureex(203);
%                         cla(hAx);
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
            dispex('Z matrix took %.1fs for %ix%i elements, %i frequencies.\n', dt, Nx_, Ny_, Nf);
        end
        function InitializeDs(this, fs)
            newfs = setdiff(fs, this.D_fs);
            if(length(newfs) < 1)
                return;
            end
            
            dispex('Calculating D for %i frequencies.\n', length(newfs));
            
            tc = tic;
            Ny_ = this.Ny;

            dy = this.unitcell.dy;
            wslot = this.unitcell.wslot;
            tlineup_ = this.tlineup;
            tlinedown_ = this.tlinedown;
            
            if(this.unitcell.walled)
                tlinedown_ = [];
            end

            z0 = Constants.z0;
            c0 = Constants.c0;

            % Create a list of all f and nypp combinations
            [fimat, nyppmat] = meshgrid(1:length(newfs), 0:Ny_-1);
            fivec = fimat(:);
            nyppvec = nyppmat(:);
            N = length(fivec);
            
            deformedpath = 0;

            %% Initialize progress bar
            hDataQueue = parallel.pool.DataQueue;
            hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
            afterEach(hDataQueue, @updateWaitbar);
            progress = -1; send(hDataQueue, nan);
            
            % deformedpath: 0 means straight integration path, 1 means deformed integration path
            for(deformedpath = 0:1)
                Dlist = cell(1,N);
                kxlist = cell(1,N);
                parfor(ni = 1:N) % PARFOR
                    fi = fivec(ni);
                    f = newfs(fi); %#ok<PFBNS> % Broadcast variable.
                    nypp = nyppvec(ni);
                    
                    lambda = Constants.c0 / f;
                    g = 5/3 * sqrt(wslot * lambda);

                    %% Determine integration path
                    k0 = 2*pi*f / c0;
                    % Go to +-inf
                    delta_ky = 0.01.*k0;
                    lim1_ky = -100.*k0-1j.*delta_ky;
                    lim2_ky = -lim1_ky;
                    % Deform integration path around branch cuts.
                    integrationpath_ky = [(-1-1j).*delta_ky, (1+1j).*delta_ky];
                    
                    %% Initial sample points
                    if(~deformedpath)
                        %% Straight path
                        % Go to +-inf
                        delta = 0.01.*k0;
                        lim1 = -100.*k0-1j.*delta;
                        lim2 = -lim1;
                        % Deform integration path around branch cuts.
                        integrationpath = [(-1-1j).*delta, (1+1j).*delta];
                        
                        % Integration path (dotted)
                        %               |
                        %               |       N1
                        %               |...............
                        % --------------:---------------
                        % .............:|\
                        %      N1       | N2
                        %               |

                        N1 = 100; % Initial number of points on the horizontal sections
                        N2 = 10;  % Initial number of points on the diagonal section through the origin
                        newkx = [linspace(lim1, integrationpath(1), N1), ...
                                 linspace(integrationpath(1), integrationpath(2), N2+2), ...
                                 linspace(integrationpath(2), lim2, N1)];
                          
                        % Remove the duplicate elements
                        newkx(N1) = [];
                        newkx(end-N1) = [];
                    else
                        %% Deformed path
                        % Go to -inf * 1j
                        delta = 0.01*k0;
                        lim1 = -5j.*k0-1.*delta;
                        lim2 = -5j.*k0+1.*delta+1.5*k0;
                        % Deform integration path around branch cuts.
                        integrationpath = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];
                        
                        % Integration path (dotted)
                        %                |      N3
                        %            N2  |.............  N2
                        %              \ :             :
                        % --------------:+--------------:-
                        %              : |              :
                        %           N1 : |              : N1
                        %              : |              :
                        %              V |              V
                        
                        N1 = 100; % Initial number of points on the vertical sections
                        N2 = 10;  % Initial number of points on the diagonal sections
                        N3 = 30; % Initial number of points on the horizontal section
                        newkx = [linspace(lim1,               integrationpath(1), N1),   ...
                                 linspace(integrationpath(1), integrationpath(2), N2+2), ...
                                 linspace(integrationpath(2), integrationpath(3), N3),   ...
                                 linspace(integrationpath(3), integrationpath(4), N2+2), ...
                                 linspace(integrationpath(4), lim2,               N1)];
                        % Remove the duplicate elements
                        newkx(N1+N2+2) = [];     % Same as newkx(N1+N2+3)
                        newkx(N1) = [];          % Same as newkx(N1+1)
                        newkx(end-N1-N2-2) = []; % Same as newkx(end-N1-N2-1)
                        newkx(end-N1) = [];      % Same as newkx(end-N1+1)
                    end

                    kx = [];
                    D = [];

%                     if(deformedpath)
%                     [hFig, hAx] = figureex;
%                     [hFig2, hAx2] = figureex;
%                     end
                    it = 0;
                    while(~isempty(newkx) && it < 10 && length(kx) < 15e3)
                        it = it + 1;
                        % Make room for the new D(newkx)
                        D = [D zeros(1, length(newkx))]; %#ok<AGROW> Cannot preallocate due to unknown number of kx.
                        % Remember the indices to be calculated.
                        kxis = length(kx)+1:(length(kx) + length(newkx));
                        % Append new kxs.
                        kx = [kx newkx]; %#ok<AGROW> Cannot preallocate due to unknown number of kx.
                        % Calculate D for the new kxs.
                        for(kxi = kxis)
                            D(kxi) = fastintegral(...
                                @(kyp) D_Integrand_kyp_old(f, dy, k0, kyp, kx(kxi), tlineup_, tlinedown_, z0, wslot, nypp), ...
                                lim1_ky, lim2_ky, 'Waypoints', integrationpath_ky);
                        end

                        % Sort kx based on the real part.
                        [kx, I] = sort(kx, 'ComparisonMethod', 'real');
                        D = D(I);
                        if(deformedpath)
                            % The sort sorts by the real part first, then the imaginary part.
                            % Since the imaginary part of the right tail is negative, it's sorted in
                            % reverse, so flip it now.
                            rightTail = logical((real(kx) > 0) & (imag(kx) <= imag(integrationpath(end))));
                            D(rightTail) = fliplr(D(rightTail));
                            kx(rightTail) = fliplr(kx(rightTail));
                        end

                        errorbound = 0.005;
                        % Calculate first and second derivative of resulting D vector.
                        % The first is an actual derivative, the second shows the change in slope
                        % between sample points, which is an indicator of how well-sampled it is.
                        if(~deformedpath)
                            dD = (D(2:end) - D(1:end-1)) ./ abs(kx(2:end) - kx(1:end-1));
                        else
                            % For the deformed path, take the exponent into account.
                            % The exponent will make it much smoother since large imaginary values
                            % of kx will pull the result to zero.
                            % The exponent -1j*kx*g is the worst case since for (x-xp) < 2g the
                            % straight deformation is used.
                            dD = (D(2:end) - D(1:end-1)) ./ abs(kx(2:end) - kx(1:end-1)) .* exp(-1j .* kx(2:end) .* 2 .* g);
                        end
                        ddD = (dD(2:end) - dD(1:end-1));
                        % Determine the points where D must be sampled more.
                        % The +1 arises from the indexing used for the derivatives.
                        ind = find(abs(ddD) > errorbound) + 1;

                        % Since we will sample after each point, add the previous points as well.
                        ind = unique([ind-1, ind]);
                        % Determine the new kx values to sample at.
                        newkx = (kx(ind) + kx(ind+1))/2;
%                         dispex('Calculating %i new values on iteration %i.\n', length(newkx), it);
    
                        %% Plots to show the calculated D
                        % Shows D (real, imag) and shows the integration path with calculated points.
                        % Updates on each iteration to show the points which are sampled well enough (black) and the points
                        % which are not (red)
                        if(0)
                            if(~deformedpath) %#ok<UNRCH>
                                ukx = kx;
                            else
                                ukx = UnfoldKVector(kx, integrationpath);
                            end
                            % Plots for debugging
                                cla(hAx);
    %                             plot(hAx, real(ukx)./k0, real(D .* exp(-1j .* kx .* 2 .* g)), 'k');
    %                             plot(hAx, real(ukx)./k0, imag(D .* exp(-1j .* kx .* 2 .* g)), 'r');
    %                             plot(hAx, real(ukx(ind))./k0, imag(D(ind) .* exp(-1j .* kx(ind) .* 2 .* g)), 'rx');
                                plot(hAx, real(ukx)./k0, real(D), 'k');
                                plot(hAx, real(ukx)./k0, imag(D), 'r');
                                plot(hAx, real(ukx(ind))./k0, imag(D(ind)), 'rx');
                                title(hAx, sprintf('nypp = %i, deformed = %i', nypp, deformedpath));
                                cla(hAx2);
                                plot(hAx2, real(kx)./k0, imag(kx)./k0, 'ko');
                                plot(hAx2, real(kx(ind))./k0, imag(kx(ind))./k0, 'ro');
                                plot(hAx2, real(ukx)./k0, imag(ukx)./k0, 'kx');
                                plot(hAx2, real(ukx(ind))./k0, imag(ukx(ind)./k0), 'rx');
                                title(hAx2, sprintf('nypp = %i, deformed = %i', nypp, deformedpath));
                        end
                        
                        % If there's no new indices, we're done.
                        if(isempty(ind))
                            break;
                        end
                    end
%                     dispex('Finished ni=%i, nypp=%i, fi=%i after calculating %i values in %i iterations (%.2fs).\n', ni, nypp, fi, length(kx), it, toc(tc));

                    Dlist{ni} = D;
                    kxlist{ni} = kx;
                    send(hDataQueue, nan);
                end
                
                if(~deformedpath)
                    % Store the frequency points in the this.D_fs vector.
                    % Only once, so do it on the first iteration, where deformedpath is false
                    this.D_fs = [this.D_fs, newfs];
                    % Update progress bar.
                    waitbar(1, hWaitbar, ...
                        {sprintf('%.1f%% Precomputing D...', 50), ...
                         sprintf('Interpolating...'), ''});
                else
                    % Update progress bar.
                    waitbar(1, hWaitbar, ...
                        {sprintf('%.1f%% Precomputing D...', 100), ...
                         sprintf('Finalizing...'), ''});
                end
                % Store the calculated Ds and their kx values in the this.Ds and this.D_kxs cells.
                for(ni = 1:N)
                    fii = find(this.D_fs == newfs(fivec(ni)));
                    nypp = nyppvec(ni);
                    
                    % The D vector is pre-interpolated on a fixed set of points to speed up later
                    % interpolation.
                    tempkxs = kxlist{ni};
                    tempDs = Dlist{ni};
                    
%                     tempkxsold = tempkxs;

                    %% Determine integration path
                    f = this.D_fs(fii);
                    k0 = 2*pi*f / c0;
                    delta = 0.01*k0;
                    
                    if(~deformedpath)
                        %% Straight integration path
                        % Go to +-inf
                        delta = 0.01.*k0;
                        lim1 = -100.*k0-1j.*delta;
                        lim2 = -lim1;
                        
                        % Determine the points on which to interpolate.
                        newkxs = linspace(real(lim1), real(lim2), 5000);
                    else
%                         if(nypp == 0)
%                             continue;
%                         end
                        %% Deformed integration path
                        lim1 = -5j.*k0-1.*delta;
                        lim2 = -5j.*k0+1.*delta+1.5*k0;
                        % Deform integration path around branch cuts.
                        integrationpath = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];
                        
%                         tempkxsold = tempkxs;
                        tempkxs = UnfoldKVector(tempkxs, integrationpath);
                        
                        % Unfold the limits of integration such that any point on the path is
                        % available in the resulting interpolated vector.
                        ulims = UnfoldKVector([lim1, lim2], integrationpath);
                        % Determine the points on which to interpolate.
                        newkxs = linspace(real(ulims(1)), real(ulims(2)), 5000);
                    end

%                     % Debug plots
%                     if(deformedpath)
%                         dkx = min(real(tempkxs(2:end) - tempkxs(1:end-1)));
%                         dispex('%i elements, %.1f ideally.\n', length(tempkxs), 2*real(lim2)/dkx);
%                         [hFig, hAx] = figureex(101); 
%                                 plot(hAx, real(tempkxs)./k0, abs(tempDs - interp1(newkxs, interp1(real(tempkxs), tempDs, newkxs), real(tempkxs)))./abs(tempDs));
%                         [hFig, hAx] = figureex(102); 
%                                 plot(hAx, real(tempkxs)./k0, real(tempDs .* exp(-1j .* tempkxsold .* 2 .* g)));
%                                 plot(hAx, real(tempkxs)./k0, imag(tempDs .* exp(-1j .* tempkxsold .* 2 .* g)), '--');
%                                 plot(hAx, real(newkxs)./k0, real(interp1(real(tempkxs), tempDs .* exp(-1j .* tempkxsold .* 2 .* g), newkxs)), ':');
%                                 plot(hAx, real(newkxs)./k0, imag(interp1(real(tempkxs), tempDs .* exp(-1j .* tempkxsold .* 2 .* g), newkxs)), '-.');
%                         [hFig, hAx] = figureex(103);
%                                 plot(hAx, real(tempkxsold), imag(tempkxsold), '.');
%                         figureex(104);
%                         unfolded = UnfoldKVector(tempkxsold, integrationpath);
%                             plot(hAx, real(unfolded), imag(unfolded), '.')
%                     end

                    tempDs = interp1(real(tempkxs), tempDs, newkxs);
                    tempkxs = newkxs;
                    

                    this.Ds{deformedpath+1, fii, nypp+1} = tempDs;
                    this.D_kxs{deformedpath+1, fii, nypp+1} = real(tempkxs);
                    this.D_interpolants{deformedpath+1, fii, nypp+1} = griddedInterpolant(this.D_kxs{deformedpath+1, fii, nypp+1}, this.Ds{deformedpath+1, fii, nypp+1});
                end
                if(~deformedpath)
                    progress = 0;
                    dispex('Finished half in %.2fs.\n', toc(tc));
                end
            end
            function updateWaitbar(~)
                progress = progress + 1;
                waitbar(progress/N, hWaitbar, ...
                    {sprintf('%.1f%% Precomputing D...', progress/N/2*100+deformedpath*50), ...
                     sprintf('%i/%i iterations done.', progress, N), ...
                     sprintf('%i/%i deformations done.', deformedpath, 2)});
            end
            delete(hWaitbar);
            dt = toc(tc);
            dispex('D took %.1fs for %i slots, %i frequencies.\n', dt, Ny_, length(newfs));
        end
        function Zas = GetInputImpedance(this, fs, excitation)
            if(size(excitation, 1) ~= this.Nx || size(excitation, 2) ~= this.Ny)
                error('Invalid excitation matrix supplied. Should be Nx by Ny.');
            end
            % Ensure the appropriate Z matrices have been calculated.
            this.InitializeZMatrix(fs);
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
            dispex('Active Z took %.1fms for %ix%i elements, %i frequencies.\n', dt*1e3, Nx_, Ny_, Nf);
        end
        function BuildCST(this, project, parentcomponent)
            if(nargin < 3 || isempty(parentcomponent))
                parentcomponent = '';
            else
                if(~strcmp(parentcomponent(end), '/'))
                    parentcomponent = [parentcomponent, '/'];
                end
            end
            componentname = [parentcomponent, 'Array'];
            
            this.unitcell.BuildCST(project, [componentname, '/UnitCell']);
            
            wcs = project.WCS();
%             obj.tlineup.BuildCST(project); 
            wcs.Enable(); 
            wcs.Store('Pre-Slot'); 
            wcs.RotateWCS('u', 180); 
             
            % Build down-stratification. 
            this.tlinedown.BuildCST(project, [componentname, '/Stratification']); 
            
            wcs.Restore('Pre-Slot'); 
            
            % Build up-stratification. 
            this.tlineup.BuildCST(project, [componentname, '/Stratification']); 
            
            wcs.Restore('Pre-Slot'); 
            wcs.Delete('Pre-Slot'); 
            wcs.Disable(); 
            
            project.StoreParameter('lambda_min', 'c0/fmin/1e9');
            project.StoreParameter('Nx', num2str(this.Nx, '%.15g'));
            project.StoreParameter('Ny', num2str(this.Ny, '%.15g'));
            project.StoreParameter('edge_distance', num2str(this.dedge*1e3, '%.15g'));
            project.StoreParameter('edge_length', '5/3 * sqr(slot_width * lambda_min)');
%             project.StoreParameter('padding_x', 'max(0, lambda_min/4-edge_length)');
%             project.StoreParameter('padding_y', 'max(0, lambda_min/4-(dy/2-slot_width/2))');
            % Round up the X padding to a multiple of dx, to ensure a good fit of the stratification
            % above the unit cells.
%                 goal_x = max(-Int(-(lambda_min/4+edge_distance)/dx)*dx, -Int(-(edge_length+edge_distance-dx/2)/dx)*dx)
            %     padding_x = goal_x - edge_length - edge_distance
            project.StoreParameter('padding_x', 'max(-Int(-(lambda_min/4)/dx)*dx, -Int(-(edge_length+edge_distance-dx/2)/dx)*dx) - edge_length - edge_distance + dx/2');
            % Round up the Y padding to a multiple of dy, to ensure a good fit of the stratification
            % above the unit cells.
%                 goal_y = max(-Int(-(lambda_min/4)/dx)*dx, -Int(-(dy/2-slot_width/2)/dx)*dx)
            %     padding_y = goal_y
            project.StoreParameter('padding_y', 'max(-Int(-(lambda_min/4)/dx)*dx, -Int(-(dy/2-slot_width/2)/dx)*dx)');
            
            transform = project.Transform();
            brick = project.Brick();
            solid = project.Solid();
            
            %% Copy the unit cell in the x-direction along the array
            % Only if there's more than 1 feed in X.
            project.NextCommandConditional('Nx > 1');
                transform.Reset();
                transform.Name([componentname, '/UnitCell']);
                transform.Vector('dx', 0, 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Nx-1');
                transform.Transform('Shape', 'Translate');
            project.NextCommandConditional('Nx > 1');
                transform.Reset();
                transform.Name([componentname, '/Stratification']);
                transform.Vector('dx', 0, 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Nx-1');
                transform.Transform('Shape', 'Translate');
                
            %% Copy the stratification above and below along the X-terminations and X-padding.
            transform.Reset();
            transform.Name([componentname, '/Stratification']);
            transform.Vector('-dx', 0, 0);
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Repetitions('-Int(0.01-(edge_distance+edge_length+padding_x-dx/2)/dx)'); % -Int(-x) = Ceiling(x)
            transform.Transform('Shape', 'Translate');
            
            transform.Reset();
            transform.Name([componentname, '/Stratification']);
            transform.Vector('dx', 0, 0);
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Repetitions('-Int(0.01-(edge_distance+edge_length+padding_x-dx/2)/dx)'); % -Int(-x) = Ceiling(x)
            transform.Transform('Shape', 'Translate');
            
            %% Create the terminations
            brick.Reset();
            brick.Name('Termination1');
            brick.Component([componentname, '/UnitCell/Slot']);
            brick.Xrange('-edge_distance-edge_length',  '-dx/2');
            brick.Yrange('-dy/2', 'dy/2');
            brick.Zrange('0', '0');
            brick.Material('PEC');
            brick.Create();
            
            brick.Reset();
            brick.Name('Termination1_Slot');
            brick.Component([componentname, '/UnitCell/Slot']);
            brick.Xrange('-edge_distance',  '-dx/2');
            brick.Yrange('-slot_width/2', 'slot_width/2');
            brick.Zrange('0', '0');
            brick.Material('PEC');
            brick.Create();
            
            brick.Reset();
            brick.Name('Termination2');
            brick.Component([componentname, '/UnitCell/Slot']);
            brick.Xrange('(Nx-1)*dx+dx/2',  '(Nx-1)*dx+edge_distance+edge_length');
            brick.Yrange('-dy/2', 'dy/2');
            brick.Zrange('0', '0');
            brick.Material('PEC');
            brick.Create();
            
            brick.Reset();
            brick.Name('Termination2_Slot');
            brick.Component([componentname, '/UnitCell/Slot']);
            brick.Xrange('(Nx-1)*dx+dx/2',  '(Nx-1)*dx+edge_distance');
            brick.Yrange('-slot_width/2', 'slot_width/2');
            brick.Zrange('0', '0');
            brick.Material('PEC');
            brick.Create();
            
            solid.Subtract([componentname, '/UnitCell/Slot:Termination1'], [componentname, '/UnitCell/Slot:Termination1_Slot']);
            solid.Subtract([componentname, '/UnitCell/Slot:Termination2'], [componentname, '/UnitCell/Slot:Termination2_Slot']);
            
            %% Copy the slot in the y-direction
            % Only if there's more than 1 slot in Y.
            project.NextCommandConditional('Ny > 1');
                transform.Reset();
                transform.Name([componentname, '/UnitCell']);
                transform.Vector(0, 'dy', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Ny-1');
                transform.Transform('Shape', 'Translate');
            project.NextCommandConditional('Ny > 1');
                transform.Reset();
                transform.Name([componentname, '/Stratification']);
                transform.Vector(0, 'dy', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Ny-1');
                transform.Transform('Shape', 'Translate');
                
            %% Copy the stratification above and below along the Y-terminations and Y-padding.
            transform.Reset();
            transform.Name([componentname, '/Stratification']);
            transform.Vector(0, '-dy', 0);
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Repetitions('-Int(-padding_y/dy)'); % -Int(-x) = Ceiling(x)
            transform.Transform('Shape', 'Translate');
            
            transform.Reset();
            transform.Name([componentname, '/Stratification']);
            transform.Vector(0, 'dy', 0);
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Repetitions('-Int(-padding_y/dy)'); % -Int(-x) = Ceiling(x)
            transform.Transform('Shape', 'Translate');
            
            %% Copy the port in the x-direction
            % Only if there's more than 1 feed in X.
            project.NextCommandConditional('Nx > 1');
                transform.Reset();
                transform.Name('port1 (SlotFeed)');
                transform.Vector('dx', 0, 0);
                transform.MultipleObjects(1);
                transform.Repetitions('Nx-1');
                transform.Transform('Port', 'Translate');

            %% Copy the port in the x-direction
            % Only if there's more than 1 slot in Y.
            project.NextCommandConditional('Ny > 1');
                % Copy the port of each nx.
                project.NextCommandLoop('nxi', '1', 'Nx')
                    transform.Reset();
                    transform.Name('port" & nxi & " (SlotFeed)');
                    transform.Vector(0, 'dy', 0);
                    transform.MultipleObjects(1);
                    transform.Repetitions('Ny-1');
                    transform.Transform('Port', 'Translate');
            
            %% Add the padding to ensure lambda/4 spacing from the slot to the open boundary
            project.NextCommandConditional('padding_x > 0');
                brick.Reset();
                brick.Name('Padding_x1');
                brick.Component(componentname);
                brick.Xrange('-edge_distance-edge_length-padding_x',  '-edge_distance-edge_length');
                brick.Yrange('-dy/2', '(Ny-1)*dy+dy/2');
                brick.Zrange('0', '0');
                brick.Material('PEC');
                brick.Create();
            project.NextCommandConditional('padding_x > 0');
                brick.Reset();
                brick.Name('Padding_x2');
                brick.Component(componentname);
                brick.Xrange('(Nx-1)*dx+edge_distance+edge_length',  '(Nx-1)*dx+edge_distance+edge_length+padding_x');
                brick.Yrange('-dy/2', '(Ny-1)*dy+dy/2');
                brick.Zrange('0', '0');
                brick.Material('PEC');
                brick.Create();
            project.NextCommandConditional('padding_y > 0');
                brick.Reset();
                brick.Name('Padding_y1');
                brick.Component(componentname);
                brick.Xrange('-edge_distance-edge_length-padding_x', '(Nx-1)*dx+edge_distance+edge_length+padding_x');
                brick.Yrange('-dy/2-padding_y',  '-dy/2');
                brick.Zrange('0', '0');
                brick.Material('PEC');
                brick.Create();
            project.NextCommandConditional('padding_y > 0');
                brick.Reset();
                brick.Name('Padding_y2');
                brick.Component(componentname);
                brick.Xrange('-edge_distance-edge_length-padding_x', '(Nx-1)*dx+edge_distance+edge_length+padding_x');
                brick.Yrange('(Ny-1)*dy+dy/2',  '(Ny-1)*dy+dy/2+padding_y');
                brick.Zrange('0', '0');
                brick.Material('PEC');
                brick.Create();
                
            % Create the vertical walls if applicable.
            if(this.unitcell.walled)
                brick.Reset()
                brick.Name('Wall')
                brick.Component(componentname);
                brick.Xrange('-edge_distance-edge_length-padding_x', '(Nx-1)*dx+edge_distance+edge_length+padding_x');
                brick.Yrange('-dy/2','-dy/2');
                brick.Zrange('-hback', '0');
                brick.Material('PEC');
                brick.Create();
                
            project.NextCommandConditional('Ny > 1');
                transform.Reset();
                transform.Name([componentname, ':Wall']);
                transform.Vector(0, 'dy', 0);
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Repetitions('Ny');
                transform.Transform('Shape', 'Translate');
            end
                
%             %% Mask the entire stratification above and below to ensure it's exactly above the padding & edges
%             brick.Reset();
%             brick.Name('Mask');
%             brick.Component(componentname);
%             brick.Xrange('-edge_distance-edge_length-padding_x', '(Nx-1)*dx+edge_distance+edge_length+padding_x');
%             brick.Yrange('-dy/2-padding_y',  '(Ny-1)*dy+dy/2+padding_y');
%             brick.Zrange(-100, 100);
%             brick.Create();
%             brick.Reset();
%             brick.Name('InverseMask');
%             brick.Component(componentname);
%             brick.Xrange('-edge_distance-edge_length-padding_x-10', '(Nx-1)*dx+edge_distance+edge_length+padding_x+10');
%             brick.Yrange('-dy/2-padding_y-10',  '(Ny-1)*dy+dy/2+padding_y+10');
%             brick.Zrange(-100, 100);
%             brick.Create();
%             
%             solid.Subtract([componentname, ':InverseMask'], [componentname, ':Mask']);
        end
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
%             x = linspace(min(basispositions)-g, max(basispositions)+g, 20);
%             x = linspace(-(max(basispositions)+g), max(basispositions)+g, 20);
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
%                 x = linspace(-0.02, 0.02, numXpoints);
                
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
                    
%                     v(xi) = fastintegral(@(kx) V_Integrand_ky(kx, 1, Nx_, Ny_, ny, dy, 1, Fs, basispositions, i, x(xi), dslot), lim1, lim2, 'Waypoints', integrationpath);
%                     v(xi) = fastintegral(@(kx) sinc(kx .* dslot ./ (2*pi)) .* exp(-1j .* kx .* x(xi)), lim1, lim2, 'Waypoints', integrationpath);
                    
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
        function ff = Farfield(this, f, excitation, ths, phs, r)
            this.InitializeDs(f);
            this.InitializeZMatrix(f);
            
            dispex('Calculating Farfield for %i angles.\n', max(length(ths), length(phs)));
            tc = tic;
            
            
            N = length(ths);
            
            % If only a single theta or phi is given, repeat it to match the length of the other.
            if(length(ths) == 1)
                ths = repmat(ths, 1, length(phs));
            end
            if(length(phs) == 1)
                phs = repmat(phs, 1, length(ths));
            end
            
            if(any(ths < -pi/2) || any(ths > pi/2))
                error('theta outside [-pi/2, pi/2] not supported.\n');
            end
            
            z0 = Constants.z0;
            
            % Only the up-transmission line is relevant since we're calculating the pattern above
            % the array.
            tlineup_ = this.tlineup;
%             tlinedown_ = this.tlinedown;
            
%             ff = zeros(1, N);
%             for(ni = 1:N) % PARFOR
%                 th = ths(ni);
%                 ph = phs(ni);
%                 
%                 [k0, kx, ky, kz] = k(f, 1, th, ph);
%                 
%                 V = this.VoltageSpectrum(f, excitation, kx, ky);
%                 
%                 Ghm_xx = GreensFunction(f, k0, kx, ky, tlineup_, [], z0);
%                 
%                 ff(ni) = 1j .* Ghm_xx .* kz .* V .* exp(-1j .* k0 .* r) ./ (2.*pi.*r);
%             end
            
            [k0s, kxs, kys, kzs] = k(f, 1, ths, phs);

            V = this.VoltageSpectrum(f, excitation, kxs, kys);

            Ghm = GreensFunction3(f, k0s, kxs, kys, tlineup_, z0);
            ff_x = 1j .* Ghm.xx .* kzs .* V .* exp(-1j .* k0s .* r) ./ (2.*pi.*r);
            ff_y = 1j .* Ghm.yx .* kzs .* V .* exp(-1j .* k0s .* r) ./ (2.*pi.*r);
            ff_z = 1j .* Ghm.zx .* kzs .* V .* exp(-1j .* k0s .* r) ./ (2.*pi.*r);
            
            ff_theta = ff_x .* cos(ths) .* cos(phs) + ff_y .* cos(ths) .* sin(phs) - ff_z .* sin(ths);
            ff_phi   = -ff_x .* sin(phs) + ff_y .* cos(phs);
            
            ff = sqrt(abs(ff_theta) .^ 2 + abs(ff_phi) .^ 2);
            dt = toc(tc);
            dispex('Farfield took %.1fs for %i angles.\n', dt, max(length(ths), length(phs)));
        end
        function ff = NormalizedFarfield(this, f, excitation, ths, phs)
            r = 1;
            ff = Farfield(this, f, excitation, ths, phs, r);
            
            ff = ff ./ max(abs(ff));
        end
    end
end

function v = v_Integrand_kx(this, f, dy, kx, x, Nx_, Ny_, ny, Fs, xs, i, deformedpath, integrationpath)
    % Retrieve the correct D index for this frequency.
    fii = find(this.D_fs == f, 1);
    
    if(~deformedpath)
        % If the path is straight, there's no need to unfold.
        interpolationkx = real(kx);
    else
        % If the path is deformed, unfold it to the real axis before interpolation.
        interpolationkx = real(UnfoldKVector(kx, integrationpath));
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
function v = Z_Integrand_kx(this, f, dy, kx, Ny_, ny, nyp, x, xp, F, Fp, deformedpath, integrationpath, tlinedown, walled, wslot)
    % Retrieve the correct D index for this frequency.
    fii = find(this.D_fs == f, 1);
    
    if(~deformedpath)
        % If the path is straight, there's no need to unfold.
        interpolationkx = real(kx);
    else
        % If the path is deformed, unfold it to the real axis before interpolation.
        interpolationkx = real(UnfoldKVector(kx, integrationpath));
    end
    
    % Interpolate the D vectors at the desired kx values.
    D = zeros(length(kx), Ny_);
    for(nypp = 0:Ny_-1)
        D(:, nypp+1) = this.D_interpolants{deformedpath+1, fii, nypp+1}(interpolationkx);
%         D(:, nypp+1) = interp1(this.D_kxs{typei, fii, nypp+1}, this.Ds{typei, fii, nypp+1}, real(kx));
    end
    
    
    int2 = zeros(1, length(kx));
    for(kxi = 1:length(kx))
        % If we have vertical walls, calculate the Ddown using Floquet modes.
        if(walled)
            Ddown = 0;
            for(my = -20:20)
                k0 = 2*pi*f/3e8;
                z0 = Constants.z0;
                
                Vtm = 1;
                Vte = 1;
                
                kym = -2*pi*my/dy;
                kr = sqrt(kx(kxi).^2 + kym.^2);
                isTE = 1;    ztedown = tlinedown.GetInputImpedance(isTE, f, k0, kr);
                isTE = 0;    ztmdown = tlinedown.GetInputImpedance(isTE, f, k0, kr);

                itedown = 1 ./ ztedown;
                itmdown = 1 ./ ztmdown;

                [Ghm] = SpectralGF.hm(z0, k0, kx(kxi), kym, Vtm, Vte, itmdown, itedown);
                Ghm_xx_down = Ghm.xx;

                Ddown = Ddown + 2*pi/dy .* Ghm_xx_down .* besselj(0, -2*pi*my/dy*wslot/2);
            end
        else
            Ddown = 0;
        end
        
        if(length(kx) == 1) % Indexing D(kxi, :) fails when Ny_ is 1.
            Dparam = D(kxi);
        else
            Dparam = D(kxi, :);
        end
        int2(kxi) = fastintegral(@(ky) Z_Integrand_ky(ky, Ny_, ny, nyp, dy, Dparam, Ddown), -pi/dy, pi./dy);
    end
    v = int2 .* Fp(kx) .* F(-kx) .* exp(-1j .* kx .* abs(x-xp));
    
%     global doplot;
%     if(doplot)
%         [~, hAx] = figureex(200);
%             plot(hAx, real(kx), imag(kx), '.');
%         k0 = 2*pi*f/3e8;
%         [~, hAx] = figureex(201);
%             pkx = UnfoldKVector(kx, integrationpath)/k0;
%             plot(hAx, real(pkx), real(int2), 'k.');
%             title(hAx, sprintf('ny %i, nyp %i', ny, nyp));
%             plot(hAx, real(pkx), imag(int2), 'r.');
%     end
end
function v = Z_Integrand_ky(ky, Ny_, ny, nyp, dy, D, Ddown)
%     sumD = zeros(size(ky));
%     for(nypp = -ny:(Ny_-1-ny))
%         sumD = sumD + D(abs(nypp)+1) .* exp(1j .* (nypp * ky)  .* dy);
%     end

    nypp = (0:Ny_-1).';
    sumD = sum(D.' .* exp(1j .* (nypp * ky)  .* dy), 1);
    
    % NOTE: Ddown is zero when there are no vertical walls since the Down stratification will then
    % be included in the D
    v = exp(-1j .* ky .* abs(ny-nyp) .* dy) ./ (sumD + Ddown);
end
function v = D_Integrand_kyp(f, dy, k0, kyp, kx, tlineup, tlinedown, z0, wslot, nypp)
    Vtm = 1;
    Vte = 1;
    
    tlineupfs = FreeSpace();
    tlinedownfs = FreeSpace();
    error
    kr = sqrt(kx.^2 + kyp.^2);
    
    % Use the real Green's function inside this interval, otherwise use FreeSpace
    % The ratio of abs(kyp)/k0 varies for different error levels
    %     Values determined for an er=3.0 ADL
    %     Error | 1e-1  | 1e-2  | 1e-3  | 1e-4  | 1e-5  | 1e-6  | 1e-7  | 1e-8  | 1e-9  | 1e-10
    %     ratio |  1.12 |  1.90 |  4.57 |  7.80 | 10.80 | 13.25 | 14.80 | 22.81 | 26.93 | 30.26 
    ratio = 10;
    iRealGF = (abs(kyp) < ratio*k0);
    zteup = zeros(size(kr));    ztmup = zteup;    ztedown = zteup;    ztmdown = zteup;
    
    isTE = 1;    zteup(iRealGF) = tlineup.GetInputImpedance(isTE, f, k0, kr(iRealGF));
    isTE = 0;    ztmup(iRealGF) = tlineup.GetInputImpedance(isTE, f, k0, kr(iRealGF));
    isTE = 1;    zteup(~iRealGF) = tlineupfs.GetInputImpedance(isTE, f, k0, kr(~iRealGF));
    isTE = 0;    ztmup(~iRealGF) = tlineupfs.GetInputImpedance(isTE, f, k0, kr(~iRealGF));

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kx, kyp, Vtm, Vte, itmup, iteup);
    Gxxup = Ghm.xx;

    isTE = 1;    ztedown(iRealGF) = tlinedown.GetInputImpedance(isTE, f, k0, kr(iRealGF));
    isTE = 0;    ztmdown(iRealGF) = tlinedown.GetInputImpedance(isTE, f, k0, kr(iRealGF));
    isTE = 1;    ztedown(~iRealGF) = tlinedownfs.GetInputImpedance(isTE, f, k0, kr(~iRealGF));
    isTE = 0;    ztmdown(~iRealGF) = tlinedownfs.GetInputImpedance(isTE, f, k0, kr(~iRealGF));

    itedown = 1 ./ ztedown;
    itmdown = 1 ./ ztmdown;

    [Ghm] = SpectralGF.hm(z0, k0, kx, kyp, Vtm, Vte, itmdown, itedown);
    Gxxdown = Ghm.xx;

    Gxx = Gxxup + Gxxdown;

    v = Gxx.*besselj(0, kyp.*wslot./2).*exp(-1j.*kyp.*nypp.*dy);
end
function v = D_Integrand_kyp_old(f, dy, k0, kyp, kx, tlineup, tlinedown, z0, wslot, nypp)
    Gxx = GreensFunction(f, k0, kx, kyp, tlineup, tlinedown, z0);
    v = Gxx.*besselj(0, kyp.*wslot./2).*exp(-1j.*kyp.*nypp.*dy);
end
function Ghm_xx = GreensFunction(f, k0, kx, ky, tlineup, tlinedown, z0)
    Vtm = 1;
    Vte = 1;
    
    kr = sqrt(kx.^2 + ky.^2);
    isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    [Ghm] = SpectralGF.hm(z0, k0, kx, ky, Vtm, Vte, itmup, iteup);
    Ghm_xx_up = Ghm.xx;

    if(~isempty(tlinedown))
        isTE = 1;    ztedown = tlinedown.GetInputImpedance(isTE, f, k0, kr);
        isTE = 0;    ztmdown = tlinedown.GetInputImpedance(isTE, f, k0, kr);

        itedown = 1 ./ ztedown;
        itmdown = 1 ./ ztmdown;

        [Ghm] = SpectralGF.hm(z0, k0, kx, ky, Vtm, Vte, itmdown, itedown);
        Ghm_xx_down = Ghm.xx;
    else
        Ghm_xx_down = 0;
    end

    Ghm_xx = Ghm_xx_up + Ghm_xx_down;
end

function Ghm = GreensFunction3(f, k0, kx, ky, tlineup, z0)
    Vtm = 1;
    Vte = 1;
    
    kr = sqrt(kx.^2 + ky.^2);
    isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    Ghm = SpectralGF.hm(z0, k0, kx, ky, Vtm, Vte, itmup, iteup);
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
function kvec = UnfoldKVector(kvec, integrationpath)
    % Unfold the integration path. We assume that integrationpath(1) is the leftmost waypoint, and
    % integrationpath(end) is the rightmost waypoint.
    % The path is unfolded by first converting the imaginary tails into real parts. This is shown in
    % the illustration below, where * are waypoints, and the paths are shown using dots
    % Old Integration path (dotted)             New Integration path (dotted)
    %                |                                         |
    %                |*...........*                            |*............*
    %                :             :                           :              :
    % --------------:+--------------*-          --------------:+---------------*............
    %              * |              :           .............* |
    %              : |              :                          |
    %              : |              :                          |
    %              V |              V                          |
    % This conversion should (hopefully) retain relative distances
    
    % Determine the tails to unfold.
    leftTail = logical((real(kvec) < 0) & (imag(kvec) < imag(integrationpath(1))));
    rightTail = logical((real(kvec) > 0) & (imag(kvec) < imag(integrationpath(end))));
    
    % Convert the relative imaginary distance to a relative real distance
    kvec(leftTail) = integrationpath(1) + imag(kvec(leftTail) - integrationpath(1));
    kvec(rightTail) = integrationpath(end) - imag(kvec(rightTail) - integrationpath(end));
end



%     pltx = kx;
%     plty = v;
%     [hFig, hAx] = figureex;
%         title(hAx, 'v');
%         hAx.ColorOrder = lines(min(size(plty)));
%         plot(real(pltx), real(plty));
%         plot(real(pltx), imag(plty), '--');
%         plot(real(pltx), abs(plty), ':');
%         alignplot(hFig, 5, 3, [], hFig.Number, 1);