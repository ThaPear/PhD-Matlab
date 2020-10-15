clear;
close all;

f0 = 10e9;
l0 = 3e8/f0;
k0 = 2 * pi * f0 / 3e8;

Nx = 32;
Ny = 32;

dx = l0 * 0.45 * 32/Nx;
dy = l0 * 0.45 * 32/Ny;



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
% th0s = [1.12345]
th0s = [0 60];
for(th0 = th0s * pi/180)
% th0 = 60 * pi/180;
ph0s = [0 45 90 135];
% ph0s = [1.23456789];
for(ph0 = ph0s * pi/180)
    if(th0 == 0 && ph0 ~= 0)
        continue;
    end

[~, kx0, ky0, ~] = k(f0, 1, th0, ph0);


Nxindex = -(Nx-1)/2:(Nx-1)/2;
Nyindex = -(Ny-1)/2:(Ny-1)/2;


% excitationtype = 'full';
% excitationtype = 'tiled';
% excitationtype = 'tiled-vertical';
% excitationtype = 'tiled-shifted';
% excitationtype = 'street-tiles';
excitationtype = 'random';

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
        betay(2:2:end, :) = 0.5*(betay(1:2:end, :) + betay(2:2:end, :));
        % Set the first to be equal to the second.
        betax(1:2:end, :) = betax(2:2:end, :);
        betay(1:2:end, :) = betay(2:2:end, :);

    case 'tiled-vertical'
        % Set of 2 tiles have the same phase, which is the average between their respective phases.

        % Start with fully sampled.
        betax = (kx0 .* Nxindex.' .* dx) * ones(1, Ny);
        betay = ones(1, Nx).' * (ky0 .* Nyindex .* dy);
        % Make every second tile the average of the first and second.
        betax(:, 2:2:end) = 0.5*(betax(:, 1:2:end) + betax(:, 2:2:end));
        betay(:, 2:2:end) = 0.5*(betay(:, 1:2:end) + betay(:, 2:2:end));
        % Set the first to be equal to the second.
        betax(:, 1:2:end) = betax(:, 2:2:end);
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
        betay(2:2:end, 1:2:end) = 0.5*(betay(1:2:end, 1:2:end) + betay(2:2:end, 1:2:end));
        betay(1:2:end, 1:2:end) = betay(2:2:end, 1:2:end);
        
        % Odd rows: Make every set of tiles the average of the first and second, starting at the
        % second tile.
        betax(2:2:end-1, 2:2:end) = 0.5*(betax(2:2:end-1, 2:2:end) + betax(3:2:end, 2:2:end));
        betax(3:2:end, 2:2:end) = betax(2:2:end-1, 2:2:end);
        betay(2:2:end-1, 2:2:end) = 0.5*(betay(2:2:end-1, 2:2:end) + betay(3:2:end, 2:2:end));
        betay(3:2:end, 2:2:end) = betay(2:2:end-1, 2:2:end);
        
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
            betay((x+1):4:end, x:4:end) = 0.5*(betay(x:4:end-1, x:4:end) + betay((x+1):4:end, x:4:end));
            betay(x:4:end-1, x:4:end) = betay((x+1):4:end, x:4:end);
        end
        
        for(x = 1:4)
            x0 = mod(x+3-1, 4)+1;
            betax(x0:4:end, (x+1):4:end) = 0.5*(betax(x0:4:end, (x+1):4:end) + betax(x0:4:end, x:4:end-1));
            betax(x0:4:end, x:4:end-1) = betax(x0:4:end, (x+1):4:end);
            betay(x0:4:end, (x+1):4:end) = 0.5*(betay(x0:4:end, (x+1):4:end) + betay(x0:4:end, x:4:end-1));
            betay(x0:4:end, x:4:end-1) = betay(x0:4:end, (x+1):4:end);
        end
        
    case 'random'

        % Ensure we get the same rand every run;
        rng(1000);
        
        % Start with fully sampled.
        betax = (kx0 .* Nxindex.' .* dx) * ones(1, Ny);
        betay = ones(1, Nx).' * (ky0 .* Nyindex .* dx);
        
        rand01 = rand(Nx,Ny);
        
        figR = zeros(Nx, Ny);
        
%         factor = 3.65; % 33.3% reduction
        factor = 2; % 46.9% reduction
        lim1 = 1/factor;
        lim2 = 2/factor;

        for i1 = 1:Nx
            for i2 = 1:Ny
                if i1 < Nx && rand01(i1,i2)>=0 && rand01(i1,i2)<lim1 && rand01(i1+1,i2)<=1
                    rand01(i1,i2)=2;
                    rand01(i1+1,i2)=2;
                    
                    
                    val = round(rand(1,1)*5);
                    figR(i1,i2)= val;
                    figR(i1+1,i2)= val;
                    
                    betax(i1, i2) = 0.5*(betax(i1, i2) + betax(i1+1, i2));
                    betax(i1+1, i2) = betax(i1, i2);
                    betay(i1, i2) = 0.5*(betay(i1, i2) + betay(i1+1, i2));
                    betay(i1+1, i2) = betay(i1, i2);
                end
                if i2 < Ny && rand01(i1,i2)>=lim1 && rand01(i1,i2)<lim2 && rand01(i1,i2+1)<=1
                    rand01(i1,i2)=4;
                    rand01(i1,i2+1)=4;
                    
                    val = round(rand(1,1)*5);
%                     val = rand(1,1)*100;
                    figR(i1,i2)=val;
                    figR(i1,i2+1)=val;
                    
                    betax(i1, i2) = 0.5*(betax(i1, i2) + betax(i1, i2+1));
                    betax(i1, i2+1) = betax(i1, i2);
                    betay(i1, i2) = 0.5*(betay(i1, i2) + betay(i1, i2+1));
                    betay(i1, i2+1) = betay(i1, i2);
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
        
        numcontrols = sum(rand01(:)==0) + sum(rand01(:)==2)/2 + sum(rand01(:) ==4)/2
    otherwise
end

%% Plot the pairs.
%{
% Attempt to generate values that are unique for each pair.
% Works only when theta is set to a small value (1) and phi is set to an irregular small value (e.g. 1.23456)
values = betax + 1000*betay;
clrs = zeros(Nx, Ny);
figureex;
    hold off; 
numcontrols = 0;
for(ny = 1:Ny)
    for(nx = 1:Nx)
        if(clrs(nx, ny) ~= 0)
            continue;
        end
        if(nx < Nx && values(nx, ny) == values(nx+1, ny))
            for(clr = [randperm(5), 6])
                if(   (nx < 2 || clrs(nx-1, ny) ~= clr) ...             % Left
                   && (ny < 2 || clrs(nx, ny-1) ~= clr) ...             % Up
                   && (nx >= Nx-1 || clrs(nx+2, ny) ~= clr) ...         % Right
                   && (nx >= Nx || ny < 2 || clrs(nx+1, ny-1) ~= clr)) % Up-right
                    clrs(nx, ny) = clr;
                    clrs(nx+1, ny) = clr;
                    numcontrols = numcontrols + 1;
                    break;
                end
                if(clr == 6)
                    error('Should not reach 5');
                end
            end
        elseif(ny < Ny && values(nx, ny) == values(nx, ny+1))
            for(clr = [randperm(5), 6])
                if(   (nx < 2 || clrs(nx-1, ny) ~= clr) ...             % Left
                   && (ny < 2 || clrs(nx, ny-1) ~= clr) ...             % Up
                   && (nx >= Nx || clrs(nx+1, ny) ~= clr) ...           % Right
                   && (ny >= Ny || nx < 2 ||  clrs(nx-1, ny+1) ~= clr)) % Down-left
                    clrs(nx, ny) = clr;
                    clrs(nx, ny+1) = clr;
                    numcontrols = numcontrols + 1;
                    break;
                end
                if(clr == 6)
                    error('Should not reach 5');
                end
            end
        else
            numcontrols = numcontrols + 1;
        end
    end
%     imagesc(clrs.');
%     axis square
%     drawnow;
%     pause(0.1)
end
dispex('Array has %i/%i controls (%.1f%%).\n', numcontrols, Nx*Ny, numcontrols/(Nx*Ny)*100);
[hFig, hAx] = figureex;
    hold(hAx, 'off');
    image(hAx, 0:Nx-1, 0:Ny-1, clrs.'+1);
    axis(hAx, 'square');
    map = lines(max(clrs(:))+2);
    map(6, :) = map(1,:);
    map(1,:) = 0;
    colormap(hAx, map);
    caxis(hAx, [0 max(clrs(:))]);
    xticks([]);
    yticks([]);
    alignplot(hFig, 5, 3, hFig.Number, [], 2);
    drawnow;
    cropfigure(hFig);
%}

%%
% numunique = length(unique(betax(:) + 1000*betay(:)))

amplitude = ones(Nx, Ny);
excitations = amplitude .* exp(1i*betax) .* exp(1i*betay);

% figure; pcolor(unwrap(angle(betax.'), [], 2)); colorbar; colormap jet;
% figure; pcolor(unwrap(angle(betay.'), [], 1)); colorbar; colormap jet;
% figure; pcolor(unwrap(angle(excitations.'), [], 2)); colorbar; colormap jet;
% figure; pcolor(angle(excitations.')); colorbar; colormap jet;


AF = zeros(size(kxs));
AFeven = zeros(size(kxs));
AFodd = zeros(size(kxs));
for(kxi = 1:size(kxs, 2))
    for(kyi = 1:size(kys, 1))
        kx = kxs(kyi, kxi);
        ky = kys(kyi, kxi);
        
        AF(kyi, kxi) = sum(sum(excitations .* (exp(-1j .* kx .* Nxindex.' .* dx) * exp(-1j .* ky .* Nyindex .* dx))));
        
        
        AFeven(kyi, kxi) = sum(sum(excitations(:,1:2:Ny) .* (exp(-1j .* kx .* Nxindex.' .* dx) * exp(-1j .* ky .* Nyindex(1:2:end) .* dx))));
        AFodd(kyi, kxi) = sum(sum(excitations(:,2:2:Ny) .* (exp(-1j .* kx .* Nxindex.' .* dx) * exp(-1j .* ky .* Nyindex(2:2:end) .* dx))));
    end
end

dth = th_lin(2) - th_lin(1);
dph = ph_lin(2) - ph_lin(1);

%%
%{
% Short dipole
r = 1;
    I0 = 1;
    IA = I0/2;
    L = l0/20;
    elementpattern_th = 1j .* k0 .* 120.*pi .* L .* IA ./ (4*pi*r) .* exp(-1j .* k0 .* r) .* -cos(ths) .* cos(phs);
    elementpattern_ph = 1j .* k0 .* 120.*pi .* L .* IA ./ (4*pi*r) .* exp(-1j .* k0 .* r) .* sin(phs);
    Za = 20*(pi * L/l0)^2;
    Pin = 1/2 * abs(I0)^2 * Za;

    elementpattern_total = sqrt(abs(elementpattern_th).^2+abs(elementpattern_ph).^2);
%}

% Isotropic source
elementpattern_total = 1;

power_element = 2*sum(sum(abs(elementpattern_total).^2./(240*pi) .* sin(ths) .* dth .* dph));
%
% power_isotropic = 2*sum(sum(1/(4*pi).^2/240/pi .* sin(ths) .* dth .* dph));
integrand = abs(elementpattern_total.*AF).^2/(240*pi) .* sin(ths) .* dth .* dph;
integrand(end, :) = [];
power = sum(sum(integrand));
power_full_broadside = 299.300860;
%%
dispex('Power for %s is %f\n', excitationtype, power);

maxdirectivity = 10*log10(4*pi*Nx*Ny*dx*dy/(l0^2));

% Total array factor
directivity = 4*pi*abs(elementpattern_total.*AF).^2/(240*pi)/power;
% Broadside power: 10.299146
[hFig, hAx] = figureex;
    hold(hAx, 'off');
    surf(kxs./k0, kys./k0, 10*log10(directivity));
    xlabel('kx / k0');
    ylabel('ky / k0');
    grid off
    box on
    view(0,90)
    set(gca, 'YDir', 'normal');
    caxis([-5 35]);
    shading interp;
    colorbar;
    colormap jet;
    hAx.Colorbar.FontSize = 10;
    hAx.Colorbar.Label.String = 'Directivity [dB]';
%     title(excitationtype);
    axis square
    hAx.LineWidth = 1;
    alignplot(hFig, 5, 3, hFig.Number, [], 2);
    cropfigure;
%     
% pks = FastPeakFind(10*log10(directivity), min(10*log10(directivity(:))));
% coords = reshape(pks, 2, []);
% values = 10*log10(directivity(sub2ind(size(directivity), coords(1,:), coords(2,:))));
% sorted = sort(values);
% a = 10;
    
% Partial array factors
%{
% Even rows
directivity = 4*pi*abs(elementpattern_total.*AFeven).^2/(240*pi)/power;
% Broadside power: 10.299146
figure;
    surf(kxs./k0, kys./k0, 10*log10(directivity));
    grid off
    box on
    view(0,90)
    set(gca, 'YDir', 'normal');
    caxis([-5 35]);
    shading interp;
    colorbar;
    colormap jet;
    title(excitationtype);
% Odd rows
directivity = 4*pi*abs(elementpattern_total.*AFodd).^2/(240*pi)/power;
% Broadside power: 10.299146
figure;
    surf(kxs./k0, kys./k0, 10*log10(directivity));
    grid off
    box on
    view(0,90)
    set(gca, 'YDir', 'normal');
    caxis([-5 35]);
    shading interp;
    colorbar;
    colormap jet;
    title('odd');
%}

% phasediff = reshape(phasediff, length(kxs_lin), length(kys_lin));
% Broadside power: 10.299146
% figureex; plot(ths(1,:).*180./pi, angle(AF(1, :) .* exp(-1j .* 2/180*pi)));
%             plot(-ths(end/2,:).*180./pi, angle(AF(end/2, :) .* exp(-1j .* 2/180*pi)));
%             yyaxis right
%           plot(ths(1,:).*180./pi, 20*log10(abs(AF(1,:))));
%           plot(-ths(end/2,:).*180./pi, 20*log10(abs(AF(end/2,:))));
%           axis equal
          
% phasediff = angle(AF .* exp(-1j .* 2/180*pi));
% [hFig, hAx] = figureex;
%     hold(hAx, 'off');
%     surf(kxs./k0, kys./k0, phasediff);
%     grid off
%     box on
%     view(0,90)
%     set(gca, 'YDir', 'normal');
% %     caxis([-5 35]);
%     shading interp;
%     colorbar;
%     colormap jet;
%     title('total');
% integrand = abs(elementpattern_total.*AFodd).^2/(240*pi) .* sin(ths) .* dth .* dph;
% integrand(end, :) = [];
% powerodd = sum(sum(integrand));
% directivityodd = 4*pi*abs(elementpattern_total.*AFodd).^2/(240*pi)/powerodd;
% [hFig, hAx] = figureex;
%     hold(hAx, 'off');
%     surf(kxs./k0, kys./k0, 10*log10(directivityodd));
%     xlabel('kx / k0');
%     ylabel('ky / k0');
%     grid off
%     box on
%     view(0,90)
%     set(gca, 'YDir', 'normal');
%     caxis([-5 35]);
%     shading interp;
%     colorbar;
%     colormap jet;
%     hAx.Colorbar.FontSize = 10;
%     hAx.Colorbar.Label.String = 'Directivity [dB]';
% %     title(excitationtype);
%     axis equal
%     hAx.LineWidth = 1;
%     alignplot(hFig, 5, 3, [], hFig.Number, 2);
%     cropfigure;
% integrand = abs(elementpattern_total.*AFeven).^2/(240*pi) .* sin(ths) .* dth .* dph;
% integrand(end, :) = [];
% powereven = sum(sum(integrand));
% directivityeven = 4*pi*abs(elementpattern_total.*AFeven).^2/(240*pi)/powereven;
% [hFig, hAx] = figureex;
%     hold(hAx, 'off');
%     surf(kxs./k0, kys./k0, 10*log10(directivityeven));
%     xlabel('kx / k0');
%     ylabel('ky / k0');
%     grid off
%     box on
%     view(0,90)
%     set(gca, 'YDir', 'normal');
%     caxis([-5 35]);
%     shading interp;
%     colorbar;
%     colormap jet;
%     hAx.Colorbar.FontSize = 10;
%     hAx.Colorbar.Label.String = 'Directivity [dB]';
% %     title(excitationtype);
%     axis equal
%     hAx.LineWidth = 1;
%     alignplot(hFig, 5, 3, [], hFig.Number, 2);
%     cropfigure;
% phasediff = angle(AFodd);
% [hFig, hAx] = figureex;
%     hold(hAx, 'off');
%     surf(kxs./k0, kys./k0, phasediff);
%     xlabel('kx / k0');
%     ylabel('ky / k0');
%     grid off
%     box on
%     view(0,90)
%     set(gca, 'YDir', 'normal');
% %     caxis([-5 35]);
%     shading interp;
%     colorbar;
%     colormap jet;
%     hAx.Colorbar.FontSize = 10;
%     hAx.Colorbar.Label.String = 'Phase';
% %     title(excitationtype);
%     axis equal
%     hAx.LineWidth = 1;
%     alignplot(hFig, 5, 3, [], hFig.Number, 2);
%     cropfigure;
% phasediff = angle(AFeven);
% [hFig, hAx] = figureex;
%     hold(hAx, 'off');
%     surf(kxs./k0, kys./k0, phasediff);
%     xlabel('kx / k0');
%     ylabel('ky / k0');
%     grid off
%     box on
%     view(0,90)
%     set(gca, 'YDir', 'normal');
% %     caxis([-5 35]);
%     shading interp;
%     colorbar;
%     colormap jet;
%     hAx.Colorbar.FontSize = 10;
%     hAx.Colorbar.Label.String = 'Phase';
%     axis equal
%     hAx.LineWidth = 1;
%     alignplot(hFig, 5, 3, [], hFig.Number, 2);
%     cropfigure;

dispex('Maximum is %f\n', max(max(10*log10(4*pi*abs(elementpattern_total.*AF).^2/(240*pi)/power))));

% [hFig, hAx] = figureex(101); plot(th_lin .* 180./pi, 10*log10(abs(AF).^2/(Nx*Ny*power_element)));
% addlegendentry(excitationtype);
% hFig.Position(1:2) = 1/2 * hFig.Position(1:2);
% hFig.Position(3:4) = 2 * hFig.Position(3:4);
% for(ifig = 1:100)
%     if(ishandle(ifig))
%         fig = figure(ifig);
%         alignplot(fig, 3, 2, [], fig.Number, 2);
%     end
% end
end
end