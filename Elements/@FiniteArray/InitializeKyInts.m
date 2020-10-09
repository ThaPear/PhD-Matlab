function InitializeKyInts(this, fs)
    newfs = setdiff(fs, this.KyInt_fs);
    if(length(newfs) < 1)
        return;
    end
    this.InitializeDs(fs);

    dispex('Ky Integral: Calculating for %i frequencies.\n', length(newfs));

    tc = tic;
    Ny_ = this.Ny;
    
    dy = this.unitcell.dy;
    wslot = this.unitcell.wslot;
    tlineup_ = this.tlineup;
    tlinedown_ = this.tlinedown;
    walled = this.unitcell.walled;

    z0 = Constants.z0;
    c0 = Constants.c0;

    % Create a list of all f and nypp combinations
    [fimat, dnymat] = meshgrid(1:length(newfs), 0:Ny_-1);
    fivec = fimat(:);
    dnyvec = dnymat(:);
    N = length(fivec);

    deformedpath = 0;

    %% Initialize progress bar
    hDataQueue = parallel.pool.DataQueue;
    hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
    afterEach(hDataQueue, @updateWaitbar);
    progress = -1; send(hDataQueue, nan);

    % deformedpath: 0 means straight integration path, 1 means deformed integration path
    for(deformedpath = 0:1)
        kyintlist = cell(1,N);
        kxlist = cell(1,N);
        parfor(ni = 1:N) % parfor
            fi = fivec(ni);
            f = newfs(fi); %#ok<PFBNS> % Broadcast variable.
            dny = dnyvec(ni);

            lambda = Constants.c0 / f;
            g = 5/3 * sqrt(wslot * lambda);

            %% Determine integration path
            k0 = 2*pi*f / c0;
            % Go to +-inf
%             delta_ky = 0.01.*k0;
            lim1_ky = -pi/dy;
            lim2_ky = -lim1_ky;
            % No deformation for now.
            integrationpath_ky = [];

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
            kyint = [];

            %{
                [hFig, hAx] = figureex;
                [hFig2, hAx2] = figureex;
            %}
            it = 0;
            while(~isempty(newkx) && it < 15 && length(kx) < 15e3)
                it = it + 1;
                % Make room for the new kyint(newkx)
                kyint = [kyint zeros(1, length(newkx))]; %#ok<AGROW> Cannot preallocate due to unknown number of kx.
                % Remember the indices to be calculated.
                kxis = length(kx)+1:(length(kx) + length(newkx));
                % Append new kxs.
                kx = [kx newkx]; %#ok<AGROW> Cannot preallocate due to unknown number of kx.
                % Determine sample points for D
                if(deformedpath)
                    interpolationkx = real(FiniteArray.UnfoldKVector(newkx, integrationpath));
                else
                    interpolationkx = real(newkx);
                end
                
                % Interpolate D at the appropriate kx values.
                fii = find(this.D_fs == f, 1);
                D = zeros(length(newkx), Ny_);
                for(nypp = 0:(Ny_-1))
                    D(:, nypp+1) = this.D_interpolants{deformedpath+1, fii, nypp+1}(interpolationkx);
                end
                kxi0 = min(kxis)-1;
                
                % Calculate ky integral for the new kxs.
                for(kxi = kxis)
                    if(Ny_ == 1) % Indexing D(kxi, :) fails when Ny_ is 1.
                        Dparam = D(kxi-kxi0);
                    else
                        Dparam = D(kxi-kxi0, :);
                    end
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
                            isTE = 1;    ztedown = tlinedown_.GetInputImpedance(isTE, f, k0, kr);
                            isTE = 0;    ztmdown = tlinedown_.GetInputImpedance(isTE, f, k0, kr);

                            itedown = 1 ./ ztedown;
                            itmdown = 1 ./ ztmdown;

                            [Ghm] = SpectralGF.hm(z0, k0, kx(kxi), kym, Vtm, Vte, itmdown, itedown);
                            Ghm_xx_down = Ghm.xx;

                            Ddown = Ddown + -1/dy .* Ghm_xx_down .* besselj(0, -2*pi*my/dy*wslot/2);
                        end
                    else
                        Ddown = 0;
                    end
                    kyint(kxi) = fastintegral(@(ky) Z_Integrand_ky(ky, Ny_, dny, dy, Dparam, Ddown), ...
                        lim1_ky, lim2_ky, 'Waypoints', integrationpath_ky);
                end

                % Sort kx based on the real part.
                [kx, I] = sort(kx, 'ComparisonMethod', 'real');
                kyint = kyint(I);
                if(deformedpath)
                    % The sort sorts by the real part first, then the imaginary part.
                    % Since the imaginary part of the right tail is negative, it's sorted in
                    % reverse, so flip it now.
                    rightTail = logical((real(kx) > 0) & (imag(kx) <= imag(integrationpath(end))));
                    kyint(rightTail) = fliplr(kyint(rightTail));
                    kx(rightTail) = fliplr(kx(rightTail));
                end

                errorbound = 0.002;
                % Calculate first and second derivative of resulting D vector.
                % The first is an actual derivative, the second shows the change in slope
                % between sample points, which is an indicator of how well-sampled it is.
                if(~deformedpath)
                    dD = (kyint(2:end) - kyint(1:end-1)) ./ abs(kx(2:end) - kx(1:end-1));
                else
                    % For the deformed path, take the exponent into account.
                    % The exponent will make it much smoother since large imaginary values
                    % of kx will pull the result to zero.
                    % The exponent -1j*kx*g is the worst case since for (x-xp) < 2g the
                    % straight deformation is used.
                    dD = (kyint(2:end) - kyint(1:end-1)) ./ abs(kx(2:end) - kx(1:end-1)) .* exp(-1j .* kx(2:end) .* 2 .* g);
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

                %% Plots to show the calculated Ky int
                % Shows Ky (real, imag) and shows the integration path with calculated points.
                % Updates on each iteration to show the points which are sampled well enough (black) and the points
                % which are not (red)
                %{
                    if(~deformedpath) %#ok<UNRCH>
                        ukx = kx;
                    else
                        ukx = FiniteArray.UnfoldKVector(kx, integrationpath);
                    end
                    % Plots for debugging
                    cla(hAx);
%                     plot(hAx, real(ukx)./k0, real(kyint .* exp(-1j .* kx .* 2 .* g)), 'k');
%                     plot(hAx, real(ukx)./k0, imag(kyint .* exp(-1j .* kx .* 2 .* g)), 'r');
%                     plot(hAx, real(ukx(ind))./k0, imag(kyint(ind) .* exp(-1j .* kx(ind) .* 2 .* g)), 'rx');
                    plot(hAx, real(ukx)./k0, real(kyint), 'ko');
                    plot(hAx, real(ukx)./k0, imag(kyint), 'kx');
                    plot(hAx, real(ukx(ind))./k0, imag(kyint(ind)), 'rx');
                    title(hAx, sprintf('dny = %i, deformed = %i', dny, deformedpath));
                    cla(hAx2);
                    plot(hAx2, real(ukx(2:end-1))./k0, abs(ddD));
%                     plot(hAx2, real(kx)./k0, imag(kx)./k0, 'ko');
%                     plot(hAx2, real(kx(ind))./k0, imag(kx(ind))./k0, 'ro');
%                     plot(hAx2, real(ukx)./k0, imag(ukx)./k0, 'kx');
%                     plot(hAx2, real(ukx(ind))./k0, imag(ukx(ind)./k0), 'rx');
                    title(hAx2, sprintf('dny = %i, deformed = %i', dny, deformedpath));
                %}

                % If there's no new indices, we're done.
                if(isempty(ind))
                    break;
                end
            end
%             dispex('Finished ni=%i, dny=%i, fi=%i after calculating %i values in %i iterations (%.2fs).\n', ni, dny, fi, length(kx), it, toc(tc));

            kyintlist{ni} = kyint;
            kxlist{ni} = kx;
            send(hDataQueue, nan);
        end

        if(~deformedpath)
            % Store the frequency points in the this.D_fs vector.
            % Only once, so do it on the first iteration, where deformedpath is false
            this.KyInt_fs = [this.KyInt_fs, newfs];
            % Update progress bar.
            waitbar(1, hWaitbar, ...
                {sprintf('%.1f%% Precomputing Ky Integral...', 50), ...
                 sprintf('Interpolating...'), ''});
        else
            % Update progress bar.
            waitbar(1, hWaitbar, ...
                {sprintf('%.1f%% Precomputing Ky Integral...', 100), ...
                 sprintf('Finalizing...'), ''});
        end
        % Store the calculated Ds and their kx values in the this.Ds and this.D_kxs cells.
        for(ni = 1:N)
            fii = find(this.KyInt_fs == newfs(fivec(ni)));
            dny = dnyvec(ni);
            
            f = this.KyInt_fs(fii);
            k0 = 2*pi*f / c0;

            % The D vector is pre-interpolated on a fixed set of points to speed up later
            % interpolation.
            tempkxs = kxlist{ni};
            tempKyInts = kyintlist{ni};

%             tempkxsold = tempkxs;
            %{
            %% Determine integration path
            delta = 0.01*k0;

            if(~deformedpath)
                %% Straight integration path
                % Go to +-inf
                lim1 = -100.*k0-1j.*delta;
                lim2 = -lim1;

                % Determine the points on which to interpolate.
                newkxs = linspace(real(lim1), real(lim2), 5000);
            else
                %% Deformed integration path
                lim1 = -5j.*k0-1.*delta;
                lim2 = -5j.*k0+1.*delta+1.5*k0;
                % Deform integration path around branch cuts.
                integrationpath = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];

%                 tempkxsold = tempkxs;
                tempkxs = FiniteArray.UnfoldKVector(tempkxs, integrationpath);

                % Unfold the limits of integration such that any point on the path is
                % available in the resulting interpolated vector.
                ulims = FiniteArray.UnfoldKVector([lim1, lim2], integrationpath);
                % Determine the points on which to interpolate.
                newkxs = linspace(real(ulims(1)), real(ulims(2)), 5000);
            end

%             % Debug plots
%             if(deformedpath)
                dkx = min(real(tempkxs(2:end) - tempkxs(1:end-1)));
                dispex('%i elements, %.1f ideally.\n', length(tempkxs), 2*real(lim2)/dkx);
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
%                 unfolded = FiniteArray.UnfoldKVector(tempkxsold, integrationpath);
%                     plot(hAx, real(unfolded), imag(unfolded), '.')
%             end

            tempKyInts = interp1(real(tempkxs), tempKyInts, newkxs);
            tempkxs = newkxs;
            %}
            %{*
                % Do not interpolate first.
                tempkxs = kxlist{ni};
                if(deformedpath)
                    % Define deformation path again.
                    delta = 0.01.*k0;
                    integrationpath = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];
                    
                    tempkxs = FiniteArray.UnfoldKVector(tempkxs, integrationpath);
                end
                tempKyInts = kyintlist{ni};
            %}

            this.KyInts{deformedpath+1, fii, dny+1} = tempKyInts;
            this.KyInt_kxs{deformedpath+1, fii, dny+1} = real(tempkxs);
            this.KyInt_interpolants{deformedpath+1, fii, dny+1} = griddedInterpolant(this.KyInt_kxs{deformedpath+1, fii, dny+1}, this.KyInts{deformedpath+1, fii, dny+1});
        end
        if(~deformedpath)
            progress = 0;
            dispex('Ky Integral: Completed straight integration path in %s.\n', fancyduration(toc(tc)));
            halfdt = toc(tc);
        else
            dispex('Ky Integral: Completed deformed integration path in %s.\n', fancyduration(toc(tc) - halfdt));
        end
    end
    function updateWaitbar(~)
        progress = progress + 1;
        waitbar(progress/N, hWaitbar, ...
            {sprintf('%.1f%% Precomputing Ky Integrals...', progress/N/2*100+deformedpath*50), ...
             sprintf('%i/%i iterations done.', progress, N), ...
             sprintf('%i/%i deformations done.', deformedpath, 2)});
    end
    delete(hWaitbar);
    dt = toc(tc);
    dispex('Ky Integral: Completed in %.1fs for %i slots, %i frequencies.\n', dt, Ny_, length(newfs));
    
end

function v = Z_Integrand_ky(ky, Ny_, dny, dy, D, Ddown)
    % NOTE: The loop to calculate sumD appears to be faster than the vector version.
    sumD = zeros(size(ky));
    for(nypp = 0:(Ny_-1))
        sumD = sumD + D(abs(nypp)+1) .* exp(1j .* (nypp * ky)  .* dy);
    end

%     nypp = (0:Ny_-1).';
%     sumD = sum(D.' .* exp(1j .* (nypp * ky)  .* dy), 1);
    
    % NOTE: Ddown is zero when there are no vertical walls since the Down stratification will then
    % be included in the D
    v = exp(-1j .* ky .* abs(dny) .* dy) ./ (sumD + Ddown);
end





































