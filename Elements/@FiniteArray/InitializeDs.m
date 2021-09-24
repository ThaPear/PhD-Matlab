function InitializeDs(this, fs)
    newfs = setdiff(fs, this.D_fs);
    if(length(newfs) < 1)
        return;
    end

    dispex('Ds: Calculating for %i frequencies, %i slots, 2 integration paths.\n', length(newfs), this.Ny);
    tc = tic;
    
    Ny_ = this.Ny;

    dy = this.unitcell.dy;
    wslot = this.unitcell.wslot;
    walled = this.unitcell.walled;
    tlineup_ = this.tlineup;
    tlinedown_ = this.tlinedown;


    z0 = Constants.z0;
    c0 = Constants.c0;

    Nf = length(newfs);

    deformedpath = 0;

    %% Initialize progress bar
    hDataQueue = parallel.pool.DataQueue;
    hWaitbar = waitbar(0, {'0% Initializing...', '', ''});
    afterEach(hDataQueue, @updateWaitbar);
    progress = -1; send(hDataQueue, nan);

    % deformedpath: 0 means straight integration path, 1 means deformed integration path
    for(deformedpath = 0)
        Dinvlist = cell(1,Nf);
        kxlist = cell(1,Nf);
        for(fi = 1:Nf)
%             dispex('Calculating %i.\n', ni);
            f = newfs(fi); %#ok<PFBNS> % Broadcast variable.

            %% Determine integration path
            k0 = 2*pi*f / c0;
            % Go to +-inf
            delta_ky = 0.01.*k0;
            lim1_ky = -5.*k0-1j.*delta_ky;
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

                N1 = 50; % Initial number of points on the horizontal sections
                N2 = 10;  % Initial number of points on the diagonal section through the origin
                newkx = [linspace(lim1, integrationpath(1), N1), ...
                         linspace(integrationpath(1), integrationpath(2), N2+2), ...
                         linspace(integrationpath(2), lim2, N1)];

                % Remove the duplicate elements
                newkx(N1) = [];
                newkx(end-N1) = [];
                realnewkx = real(newkx);
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
                realnewkx = real(FiniteArray.UnfoldKVector(newkx, integrationpath));
            end
            
            [kxadapt, realkxadapt, Dinvadapt] = adaptivecalc(@(kx) FiniteArray.CalculateDinv(f, dy, k0, kx, tlineup_, tlinedown_, z0, wslot, Ny_, walled), newkx, realnewkx);
            
            %{
            % Plot all Dinv.
            hFig = figureex;
            for(ny = 1:Ny_)
                for(nyp = 1:Ny_)
                    hAx = subplot(Ny_, Ny_, ny+(nyp-1)*Ny_);
                    hold(hAx, 'on');
                    grid(hAx, 'on');
                    box(hAx, 'on');
                    plot(real(kxadapt)./k0, real(Dinvadapt(:, ny, nyp)));
                    plot(real(kxadapt)./k0, imag(Dinvadapt(:, ny, nyp)));
                end
            end
            %}
            
            Dinvlist{fi} = Dinvadapt;
            kxlist{fi} = kxadapt;
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
        for(fi = 1:Nf)
            f = newfs(fi);
            k0 = 2*pi*f/c0;
            fii = find(this.D_fs == newfs(fi), 1);

            % The D vector is pre-interpolated on a fixed set of points to speed up later
            % interpolation.
            tempkxs = kxlist{fi};
            tempDs = Dinvlist{fi};
            
            if(~deformedpath)
                realtempkxs = real(tempkxs);
            else
                delta = 0.01*k0;
                integrationpath = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];
                realtempkxs = real(FiniteArray.UnfoldKVector(tempkxs, integrationpath));
            end

            for(ny = 0:Ny_-1)
                for(nyp = 0:Ny_-1)
                    this.Dinv_interpolants{deformedpath+1, fii, ny+1, nyp+1} = griddedInterpolant(realtempkxs, tempDs(:, ny+1, nyp+1));
                end
            end
            this.Dinvs{deformedpath+1, fii} = tempDs;
            this.Dinv_kxs{deformedpath+1, fii} = realtempkxs;
        end
        if(~deformedpath)
            progress = 0;
            dispex('Ds: Completed straight integration path in %s.\n', fancyduration(toc(tc)));
            halfdt = toc(tc);
        else
            dispex('Ds: Completed deformed integration path in %s.\n', fancyduration(toc(tc) - halfdt));
        end
    end
    
    function updateWaitbar(~)
        progress = progress + 1;
        waitbar(progress/Nf, hWaitbar, ...
            {sprintf('%.1f%% Precomputing D...', progress/Nf/2*100+deformedpath*50), ...
             sprintf('%i/%i iterations done.', progress, Nf), ...
             sprintf('%i/%i integration paths done.', deformedpath, 2)});
    end
    delete(hWaitbar);
    
    dt = toc(tc);
    dispex('Ds: Completed in %s.\n', fancyduration(dt));
end