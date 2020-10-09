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

%             if(deformedpath)
%             [hFig, hAx] = figureex;
%             [hFig2, hAx2] = figureex;
%             end
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
%                 dispex('Calculating %i new values on iteration %i.\n', length(newkx), it);

                %% Plots to show the calculated D
                % Shows D (real, imag) and shows the integration path with calculated points.
                % Updates on each iteration to show the points which are sampled well enough (black) and the points
                % which are not (red)
                if(0)
                    if(~deformedpath) %#ok<UNRCH>
                        ukx = kx;
                    else
                        ukx = this.UnfoldKVector(kx, integrationpath);
                    end
                    % Plots for debugging
                        cla(hAx);
%                     plot(hAx, real(ukx)./k0, real(D .* exp(-1j .* kx .* 2 .* g)), 'k');
%                     plot(hAx, real(ukx)./k0, imag(D .* exp(-1j .* kx .* 2 .* g)), 'r');
%                     plot(hAx, real(ukx(ind))./k0, imag(D(ind) .* exp(-1j .* kx(ind) .* 2 .* g)), 'rx');
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
%             dispex('Finished ni=%i, nypp=%i, fi=%i after calculating %i values in %i iterations (%.2fs).\n', ni, nypp, fi, length(kx), it, toc(tc));

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

%             tempkxsold = tempkxs;

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
%                 if(nypp == 0)
%                     continue;
%                 end
                %% Deformed integration path
                lim1 = -5j.*k0-1.*delta;
                lim2 = -5j.*k0+1.*delta+1.5*k0;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];

%                 tempkxsold = tempkxs;
                tempkxs = this.UnfoldKVector(tempkxs, integrationpath);

                % Unfold the limits of integration such that any point on the path is
                % available in the resulting interpolated vector.
                ulims = this.UnfoldKVector([lim1, lim2], integrationpath);
                % Determine the points on which to interpolate.
                newkxs = linspace(real(ulims(1)), real(ulims(2)), 5000);
            end

%             % Debug plots
%             if(deformedpath)
%                 dkx = min(real(tempkxs(2:end) - tempkxs(1:end-1)));
%                 dispex('%i elements, %.1f ideally.\n', length(tempkxs), 2*real(lim2)/dkx);
%                 [hFig, hAx] = figureex(101); 
%                         plot(hAx, real(tempkxs)./k0, abs(tempDs - interp1(newkxs, interp1(real(tempkxs), tempDs, newkxs), real(tempkxs)))./abs(tempDs));
%                 [hFig, hAx] = figureex(102); 
%                         plot(hAx, real(tempkxs)./k0, real(tempDs .* exp(-1j .* tempkxsold .* 2 .* g)));
%                         plot(hAx, real(tempkxs)./k0, imag(tempDs .* exp(-1j .* tempkxsold .* 2 .* g)), '--');
%                         plot(hAx, real(newkxs)./k0, real(interp1(real(tempkxs), tempDs .* exp(-1j .* tempkxsold .* 2 .* g), newkxs)), ':');
%                         plot(hAx, real(newkxs)./k0, imag(interp1(real(tempkxs), tempDs .* exp(-1j .* tempkxsold .* 2 .* g), newkxs)), '-.');
%                 [hFig, hAx] = figureex(103);
%                         plot(hAx, real(tempkxsold), imag(tempkxsold), '.');
%                 figureex(104);
%                 unfolded = this.UnfoldKVector(tempkxsold, integrationpath);
%                     plot(hAx, real(unfolded), imag(unfolded), '.')
%             end

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