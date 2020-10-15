clear;
% close all;

f0 = 10e9;
l0 = 3e8/f0;
k0 = 2 * pi * f0 / 3e8;

Nx = 32;
Ny = 32;

dx = l0 * 0.45 * 32/Nx;
dy = l0 * 0.45 * 32/Ny;



% Observation points
th_lin = linspace(0, pi/2, 300);
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


th0s = linspace(0, 90, 20);
Dmax = zeros(1, length(th0s));
for(ith0 = 1:length(th0s))
    th0 = th0s(ith0) * pi/180;
    
% th0 = 60 * pi/180;
ph0s = [0 45 90 135];
ph0s = [10];
for(ph0 = ph0s * pi/180)

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
        
    case 'random'

        % Ensure we get the same rand every run;
        rng(1000);
        
        % Start with fully sampled.
        betax = (kx0 .* Nxindex.' .* dx) * ones(1, Ny);
        betay = ones(1, Nx).' * (ky0 .* Nyindex .* dx);
        
        rand01 = rand(Nx,Ny);
        
        figR = zeros(Nx, Ny);
        
        lim1 = 1/3;
        lim2 = 2/3;

        for i1 = 1:Nx-1
            for i2 = 1:Ny-1
                if rand01(i1,i2)>=0 && rand01(i1,i2)<lim1 && rand01(i1+1,i2)<=1
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
                if rand01(i1,i2)>=lim1 && rand01(i1,i2)<lim2 && rand01(i1,i2+1)<=1
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
        
        % length(unique(betax(:) + 1000*betay(:)))
        numcontrols = sum(rand01(:)==0) + sum(rand01(:)==2)/2 + sum(rand01(:) ==4)/2
    otherwise
end

amplitude = ones(Nx, Ny);
excitations = amplitude .* exp(1i*betax) .* exp(1i*betay);

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

factor = (4*pi)*2*120*pi;

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
% dispex('Power for %s is %f\n', excitationtype, power);

maxdirectivity = 10*log10(4*pi*Nx*Ny*dx*dy/(l0^2));


% Total array factor
directivity = 4*pi*abs(elementpattern_total.*AF).^2/(240*pi)/power;
Dmax(ith0) = max(directivity(:));

% dispex('Maximum is %f\n', max(max(10*log10(4*pi*abs(elementpattern_total.*AF).^2/(240*pi)/power))));

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



figureex; plot(th0s, 10*log10(Dmax));