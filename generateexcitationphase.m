function [excitation, groups] = generateexcitationphase(Nx, Ny, fs, dx, dy, th, ph, excitationtype)
    % Possible excitation types: 'full', 'tiles', 'tiles-vertical', 'tiles-shifted', 'street-tiles', 'random'
    if(nargin < 8)
        excitationtype = 'full';
    end
    groups = {};
    dispex('Generating ''%s'' excitation phases.\n', excitationtype);
    if(~isinf(Nx) && ~isinf(Ny))
        Nf = length(fs);
        % Old fully sampled
%         excitation = zeros(Nx, Ny, Nf);
%         for(fi = 1:Nf)
%             f = fs(fi);
%             [~, kx0, ky0, ~] = k(f, 1, th, ph);
%             excitation(:, :, fi) = exp(-1j .* kx0 .* (1:Nx) .* dx).' * exp(-1j .* ky0 .* (1:Ny) .* dy);
%         end
        % Various different excitation types.
        Nxindex = -(Nx-1)/2:(Nx-1)/2;
        Nyindex = -(Ny-1)/2:(Ny-1)/2;
        excitation = zeros(Nx, Ny, Nf);
        for(fi = 1:Nf)
            f = fs(fi);
            [~, kx0, ky0, ~] = k(f, 1, th, ph);
            % Start with fully sampled.
            betax = (kx0 .* Nxindex.' .* dx) * ones(1, Ny);
            betay = ones(1, Nx).' * (ky0 .* Nyindex .* dy);
            switch(excitationtype)
                case 'full'
                    % No changes from fully sampled.
                case 'tiles'
                    % Set of 2 tiles have the same phase, which is the average between their respective phases.
                    % +---+---+---+---+---+---+---+---+
                    % |       |       |       |       |
                    % +---+---+---+---+---+---+---+---+
                    % |       |       |       |       |
                    % +---+---+---+---+---+---+---+---+
                    % |       |       |       |       |
                    % +---+---+---+---+---+---+---+---+

                    % Make every second tile the average of the first and second.
                    betax(2:2:end, :) = 0.5*(betax(1:2:end-1, :) + betax(2:2:end, :));
                    % Set the first to be equal to the second.
                    betax(1:2:end-1, :) = betax(2:2:end, :);

                    groups = cell(1, floor(Nx/2)*Ny);
                    ind = 1;
                    for(nx = 1:2:Nx-1)
                        for(ny = 1:Ny)
                            groups{ind} = [nx ny nx+1 ny];
                            ind = ind + 1;
                        end
                    end
                    if(ind ~= floor(Nx/2)*Ny+1)
                        error('Incorrect cell size.');
                    end
                case 'tiles-vertical'
                    % Set of 2 tiles have the same phase, which is the average between their respective phases.

                    % Make every second tile the average of the first and second.
                    betay(:, 2:2:end) = 0.5*(betay(:, 1:2:end) + betay(:, 2:2:end));
                    % Set the first to be equal to the second.
                    betay(:, 1:2:end) = betay(:, 2:2:end);
                    
                    groups = cell(1, Nx*floor(Ny/2));
                    ind = 1;
                    for(nx = 1:Nx)
                        for(ny = 1:2:Ny-1)
                            groups{ind} = [nx ny nx ny+1];
                            ind = ind + 1;
                        end
                    end
                    if(ind ~= Nx*floor(Ny/2)+1)
                        error('Incorrect cell size.');
                    end
                case 'tiles-shifted'
                    % Sets of 2 tiles have the same phase, which is the average between their respective phases.
                    % Every row in y is shifted by 1 tile versus the previous one.
                    % +---+---+---+---+---+---+---+---+
                    % |       |       |       |       |
                    % +---+---+---+---+---+---+---+---+
                    %     |       |       |       |    
                    % +---+---+---+---+---+---+---+---+
                    % |       |       |       |       |
                    % +---+---+---+---+---+---+---+---+

                    % Even rows: Make every set of tiles the average of the first and second.
                    betax(2:2:end, 1:2:end) = 0.5*(betax(1:2:end, 1:2:end) + betax(2:2:end, 1:2:end));
                    betax(1:2:end, 1:2:end) = betax(2:2:end, 1:2:end);

                    % Odd rows: Make every set of tiles the average of the first and second, starting at the
                    % second tile.
                    betax(2:2:end-1, 2:2:end) = 0.5*(betax(2:2:end-1, 2:2:end) + betax(3:2:end, 2:2:end));
                    betax(3:2:end, 2:2:end) = betax(2:2:end-1, 2:2:end);

                    groups = cell(1, floor(Nx/2)*floor(Ny/2) + floor((Nx-1)/2)*ceil(Ny/2));
                    ind = 1;
                    for(nx = 1:2:Nx-1)
                        for(ny = 1:Ny)
                            if(~mod(ny, 2)) % Even rows
                                groups{ind} = [nx ny nx+1 ny];
                            else % Odd rows
                                if(nx == Nx-1)
                                    continue;
                                end
                                groups{ind} = [nx+1 ny nx+2 ny];
                            end
                            ind = ind + 1;
                        end
                    end
                    if(ind ~= floor(Nx/2)*floor(Ny/2) + floor((Nx-1)/2)*ceil(Ny/2)+1)
                        error('Incorrect cell size.');
                    end
                case 'street-tiles'
                    % +---+---+---+   +---+---+---+   +
                    % |       |   |   |       |   |   |
                    % +---+---+   +---+---+---+   +---+
                    %     |   |   |       |   |   |    
                    % +---+   +---+---+---+   +---+---+
                    % |   |   |       |   |   |       |
                    % +   +---+---+---+   +---+---+---+
                    % |   |       |   |   |       |   |
                    % +---+---+---+   +---+---+---+   +
                    % |       |   |   |       |   |   |
                    % +---+---+   +---+---+---+   +---+

                    % 1st row: Make every second set of 2 tiles a pair
            %         betax(2:4:end, 1:4:end) = 1;%0.5*(betax(1:4:end, 1:4:end) + betax(2:4:end, 1:4:end));
            %         betax(1:4:end, 1:4:end) = 1;%betax(2:4:end, 1:4:end);
                    for(x = 1:4)
                        betax((x+1):4:end, x:4:end) = 0.5*(betax(x:4:end-1, x:4:end) + betax((x+1):4:end, x:4:end));
                        betax(x:4:end-1, x:4:end) = betax((x+1):4:end, x:4:end);
                    end

                    for(x = 1:4)
                        x0 = mod(x+3-1, 4)+1;
                        betax(x0:4:end, (x+1):4:end) = 0.5*(betax(x0:4:end, (x+1):4:end) + betax(x0:4:end, x:4:end-1));
                        betax(x0:4:end, x:4:end-1) = betax(x0:4:end, (x+1):4:end);
                    end

                    for(x = 1:4)
                        betay((x+1):4:end, x:4:end) = 0.5*(betay(x:4:end-1, x:4:end) + betay((x+1):4:end, x:4:end));
                        betay(x:4:end-1, x:4:end) = betay((x+1):4:end, x:4:end);
                    end

                    for(x = 1:4)
                        x0 = mod(x+3-1, 4)+1;
                        betay(x0:4:end, (x+1):4:end) = 0.5*(betay(x0:4:end, (x+1):4:end) + betay(x0:4:end, x:4:end-1));
                        betay(x0:4:end, x:4:end-1) = betay(x0:4:end, (x+1):4:end);
                    end
                    
                    warning('Groups not implemented for street-tiles.');
                case 'random'
                    % Ensure we get the same rand every run.
                    rng(1000);

                    rand01 = rand(Nx,Ny);

                    figR = zeros(Nx, Ny);

%                     factor = 3.65; % 33.3% reduction
%                     factor = 2; % 46.9% reduction
                    factor = 6; % 24.4% reduction
                    factor = 8; % 20.1% reduction
                    factor = 13; % 15.1% reduction
                    factor = 18.5; % % reduction
                    lim1 = 1/factor;
                    lim2 = 2/factor;

                    groups = {};
                    for nx = 1:Nx
                        for ny = 1:Ny
                            % Ensure we can make a horizontal tile.
                            if nx < Nx && rand01(nx,ny)>=0 && rand01(nx,ny)<lim1 && rand01(nx+1,ny)<=1
                                % Make horizontal tile.
                                rand01(nx,ny)=2;
                                rand01(nx+1,ny)=2;


                                val = round(rand(1,1)*5);
                                figR(nx,ny)= val;
                                figR(nx+1,ny)= val;

                                betax(nx, ny) = 0.5*(betax(nx, ny) + betax(nx+1, ny));
                                betax(nx+1, ny) = betax(nx, ny);
                                betay(nx, ny) = 0.5*(betay(nx, ny) + betay(nx+1, ny));
                                betay(nx+1, ny) = betay(nx, ny);
                                groups{end+1} = [nx ny nx+1 ny]; %#ok<AGROW>
                            end
                            % Ensure we can make a vertical tile.
                            if ny < Ny && rand01(nx,ny)>=lim1 && rand01(nx,ny)<lim2 && rand01(nx,ny+1)<=1
                                % Make vertical tile.
                                rand01(nx,ny)=4;
                                rand01(nx,ny+1)=4;

                                val = round(rand(1,1)*5);
            %                     val = rand(1,1)*100;
                                figR(nx,ny)=val;
                                figR(nx,ny+1)=val;

                                betax(nx, ny) = 0.5*(betax(nx, ny) + betax(nx, ny+1));
                                betax(nx, ny+1) = betax(nx, ny);
                                betay(nx, ny) = 0.5*(betay(nx, ny) + betay(nx, ny+1));
                                betay(nx, ny+1) = betay(nx, ny);
                                groups{end+1} = [nx ny nx ny+1]; %#ok<AGROW>
                            end
                        end
                    end
                    rand01(rand01<=1) = 0;

            %         figureex
            %         hold off;
            %         pcolor(rand01)
            %         axis square
            %         colormap(jet)
            %         colormap([0 0 1; 0 1 0; 1 0 0])
            %         caxis([0,4])
            %         
            %         figureex
            %         hold off;
            %         pcolor(figR)
            %         colormap jet

                    numcontrols = sum(rand01(:)==0) + sum(rand01(:)==2)/2 + sum(rand01(:) ==4)/2;
                    if(fi == Nf)
                        dispex('Generated random tiling with a %.1f%% reduction.\n', 100*(1-numcontrols/(Nx*Ny)));
                    end
                otherwise
                    error('Invalid excitation type ''%s''.\n', excitationtype);
            end

            excitation(:, :, fi) = exp(-1i*betax) .* exp(-1i*betay);
        end
    else
        error('Not implemented.\n');
    end
    
    groups = cell2mat(groups.');
end