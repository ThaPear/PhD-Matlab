clear;
% close all;
for(i = 1:99)
    if(ishandle(i))
        close(i);
    end
end


f0 = 31e9;
l0 = 3e8/f0;
k0 = 2 * pi * f0 / 3e8;

Nx = 10;
Ny = 10;

dx = l0 * 0.45;
dy = l0 * 0.45;



% Observation points
th_lin = linspace(0, pi/2, 310);
% th_lin = linspace(0, -pi/2, 300);
% ph_lin = [0 pi/2];
Nph = 300;
ph_lin = linspace(0, 2*pi, Nph);
[ths, phs] = meshgrid(th_lin, ph_lin);

[~, kxs, kys, ~] = k(f0, 1, ths, phs);

kxs_lin = kxs(1,:);
kys_lin = kys(1,:);

% Scanning angles
% th0s = [eps, 20, 40, 60];
th0s = [60];
for(th0 = th0s * pi/180)
% th0 = 60 * pi/180;
ph0 = 90 * pi/180;


folder = 'h:\Git\PhD-Matlab\Validations\finitearray_ESA\';
[parameters, S] = CST.LoadData([folder, '10x10_esa.s100p']);
fsCST = parameters.frequencies;

zfeed = 80;
for(fi = 1:length(fsCST))
    f0 = fsCST(fi);


    [~, kx0, ky0, ~] = k(f0, 1, th0, ph0);


    Nxindex = -(Nx-1)/2:(Nx-1)/2;
    Nyindex = -(Ny-1)/2:(Ny-1)/2;


    excitationtype = 'full';
%     excitationtype = 'tiled';
%     excitationtype = 'tiled-vertical';
%     excitationtype = 'tiled-shifted';
%     excitationtype = 'street-tiles';

    switch(excitationtype)
        case 'full'
            betax = (kx0 .* Nxindex.' .* dx) * ones(1, Ny);
            betay = ones(1, Nx).' * (ky0 .* Nyindex .* dy);

        case 'tiled'
            % Set of 2 tiles have the same phase, which is the average between their respective phases.
            % +---+---+---+---+---+---+---+---+
            % |   :   |   :   |   :   |   :   |
            % +---+---+---+---+---+---+---+---+
            % |   :   |   :   |   :   |   :   |
            % +---+---+---+---+---+---+---+---+
            % |   :   |   :   |   :   |   :   |
            % +---+---+---+---+---+---+---+---+

            % Start with fully sampled.
            betax = (kx0 .* Nxindex.' .* dx) * ones(1, Ny);
            betay = ones(1, Nx).' * (ky0 .* Nyindex .* dy);
            % Make every second tile the average of the first and second.
            betax(2:2:end, :) = 0.5*(betax(1:2:end, :) + betax(2:2:end, :));
            % Set the first to be equal to the second.
            betax(1:2:end, :) = betax(2:2:end, :);

        case 'tiled-vertical'
            % Set of 2 tiles have the same phase, which is the average between their respective phases.

            % Start with fully sampled.
            betax = (kx0 .* Nxindex.' .* dx) * ones(1, Ny);
            betay = ones(1, Nx).' * (ky0 .* Nyindex .* dy);
            % Make every second tile the average of the first and second.
            betay(:, 2:2:end) = 0.5*(betay(:, 1:2:end) + betay(:, 2:2:end));
            % Set the first to be equal to the second.
            betay(:, 1:2:end) = betay(:, 2:2:end);


        case 'tiled-shifted'
            % Sets of 2 tiles have the same phase, which is the average between their respective phases.
            % Every row in y is shifted by 1 tile versus the previous one.
            % +---+---+---+---+---+---+---+---+
            % |   :   |   :   |   :   |   :   |
            % +---+---+---+---+---+---+---+---+
            % :   |   :   |   :   |   :   |   :
            % +---+---+---+---+---+---+---+---+
            % |   :   |   :   |   :   |   :   |
            % +---+---+---+---+---+---+---+---+
            % Start with fully sampled.
            betax = (kx0 .* Nxindex.' .* dx) * ones(1, Ny);
            betay = ones(1, Nx).' * (ky0 .* Nyindex .* dy);

            % Even rows: Make every set of tiles the average of the first and second.
            betax(2:2:end, 1:2:end) = 0.5*(betax(1:2:end, 1:2:end) + betax(2:2:end, 1:2:end));
            betax(1:2:end, 1:2:end) = betax(2:2:end, 1:2:end);

            % Odd rows: Make every set of tiles the average of the first and second, starting at the
            % second tile.
            betax(2:2:end-1, 2:2:end) = 0.5*(betax(2:2:end-1, 2:2:end) + betax(3:2:end, 2:2:end));
            betax(3:2:end, 2:2:end) = betax(2:2:end-1, 2:2:end);

        case 'street-tiles'
            % +---+---+---+   +---+---+---+   +
            % |   :   |   |   |   :   |   |   |
            % +---+---+   +---+---+---+   +---+
            % :   |   |   |   :   |   |   |   :
            % +---+   +---+---+---+   +---+---+
            % |   |   |   :   |   |   |   :   |
            % +   +---+---+---+   +---+---+---+
            % |   |   :   |   |   |   :   |   |
            % +---+---+---+   +---+---+---+   +
            % |   :   |   |   |   :   |   |   |
            % +---+---+   +---+---+---+   +---+
            % Start with fully sampled.
            betax = (kx0 .* Nxindex.' .* dx) * ones(1, Ny);
            betay = ones(1, Nx).' * (ky0 .* Nyindex .* dx);

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


        otherwise
    end

    amplitude = ones(Nx, Ny);
    excitations(fi, :,:) = amplitude .* exp(1i*betax) .* exp(1i*betay);
end
% figure; pcolor(unwrap(angle(betax.'), [], 2)); colorbar; colormap jet;
% figure; pcolor(unwrap(angle(betay.'), [], 1)); colorbar; colormap jet;
% figure; pcolor(unwrap(angle(excitations.'), [], 2)); colorbar; colormap jet;
% figure; pcolor(angle(excitations.')); colorbar; colormap jet;

%%
Z = zeros(size(S));
Zas = zeros(Nx*Ny, length(fsCST));
for(fi = 1:length(fsCST))
    Z(:,:,fi) = inv(eye(100)-S(:,:,fi)) * (eye(100)+S(:,:,fi)) .* 80;
    
    Zmat_ = Z(:,:,fi);
    % Build a matrix with zfeed on the diagonal, but only for non-termination elements
    Zl = repmat([ones(1, Nx)*zfeed], 1, Ny);
    Zlmat = diag(Zl);
    % Add the loads to the mutual impedances
    Zp = Zmat_ + Zlmat;
    % Determine current through the elements
    i = Zp\squeeze(excitations(fi, :)).';
    % Determine voltage on the elements themselves.
    v = Zmat_ * i;

    Zas(:, fi) = v./i;
end


portindexing = zeros(Nx,Ny);
portindexing(:,1) = 1:Nx;
for(nx = 1:Nx)
    portindexing(nx, 2:end) = Nx + (nx-1)*(Ny-1) + (1:Ny-1);
end

%%

hFig = figureex;
for(nx = 1:Nx)
    for(ny = 1:Ny)
        idx = portindexing(nx, ny);
        hAx = subplot(Nx, Ny, ny+(nx-1)*Ny);
        hold(hAx, 'on');
        grid(hAx, 'on');
        box(hAx, 'on');
        plot(hAx, fsCST/1e9, squeeze(real(Zas(idx, :))), 'k');
        plot(hAx, fsCST/1e9, squeeze(imag(Zas(idx, :))), 'k--');
    end
end
legend(hAx, {'reCST', 'imCST', 'reMATLAB', 'imMATLAB'});

Sact = (Zas - zfeed) ./ (Zas + zfeed);

Sfixed = zeros(Nx, Ny, length(fsCST));
for(nx = 1:Nx)
    for(ny = 1:Ny)
        idx = portindexing(nx, ny);
        Sfixed(nx, ny, :) = Sact(idx, :);
    end
end

for(fi = 1:length(fsCST))
    Pin = sum(abs(excitations(fi, :)).^2);
    Prefl(fi) = sum(sum(abs(Sfixed(:,:,fi) .* squeeze(excitations(fi, :, :))).^2));
    efficiency(fi) = (1-Prefl(fi)./Pin);
end

figureex(100); plot(fsCST, efficiency);

[hFig, hAx] = figureex;
for(nx = 1:Nx)
    for(ny = 1:Ny)
        idx = portindexing(nx, ny);
%         hAx = subplot(Nx, Ny, ny+(nx-1)*Ny);
%         hold(hAx, 'on');
%         grid(hAx, 'on');
%         box(hAx, 'on');
        plot(hAx, fsCST/1e9, 20*log10(abs(Sact(idx, :))));
    end
end
legend(hAx, {'reCST', 'imCST', 'reMATLAB', 'imMATLAB'});

% 
for(ifig = 1:100)
    if(ishandle(ifig))
        fig = figure(ifig);
        alignplot(fig, 3, 2, [], fig.Number, 2);
    end
end
end