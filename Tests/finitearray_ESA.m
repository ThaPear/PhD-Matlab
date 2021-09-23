% close all;
closewaitbars;
clearvars -except array;
SetupPath;
clear global;

f0 = 31e9;
c0 = Constants.c0;
l0 = c0/f0;
z0 = Constants.z0;

fs = (12:1:32) * 1e9;
% fs = (12:10:32) * 1e9;
% fs = [fs, (28:1:32)*1e9];
% fs = sort([fs, 31e9]);

th = eps+60 * pi/180;
ph =     0 * pi/180;

% z1 = 90;
% z2 = z0;
% zfeed = 80;
% 
% dx = 4.35e-3; %4.35e-3; %4.354838709677e-3; %c0/f0 * 0.45;
% dy = dx;
% erback = 2.2;
% hback = 0.9e-3;
% hgap = -0.1e-3;
% 
% dslot = 0.3e-3;
% wfeed = 0.3e-3;
% wslot = 0.8e-3;
% 
dedge = 0.25*l0;
% walled = 1;
% 
% C = 0.2e-12;


% p = dx / 2;
% gamma = 0.2;
% N = 2;
% f0match = 19e9;
% f0design = 29e9;
% slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);
% 
% % Apply the hgap parameter.
% slab.elements{1}.elements{1}.L = slab.elements{1}.elements{1}.L + hgap;


c0 = Constants.c0;
z0 = Constants.z0;
f0 = 31e9;          % 13.75 to 31 GHz.

z1 = 80;
z2 = z0;
zfeed = 80;

dx = 4.35e-3; %4.35e-3; %4.354838709677e-3; %c0/f0 * 0.45;
dy = dx;
erback = 2.2;
hback = 1.9e-3;%1.9496e-3; % = (c0/f0)/sqrt(erback*0.7)/4
wslot = 1.4e-3;
dslot = 2e-3;
walled = 1;

C = inf;%0.2e-12;

p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 19e9;
f0design = 29e9;
slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);

tlineup = TerminatedTLine(slab, FreeSpace());
tlinedown = ShortedLine(erback, hback, erback);

% Number of unit cells.
Nx = 32;
Ny = 32;
Nf = length(fs);

%% Generate excitation
% Possible excitation types: 'full', 'tiles', 'tiles-vertical', 'tiles-shifted', 'street-tiles', 'random'
[excitation, groups] = generateexcitationphase(Nx, Ny, fs, dx, dy, th, ph, 'full');

%% HFSS Excitation
%{
%         matfile = load('e:\ HFSS\Feed with phase v2\60 degrees\ Results\s1to1024-lossy.mat');
        matfile = load('e:\ HFSS\Feed with phase v2\60 degrees stripline\ Results\s1to1024-lossy.mat');
%         matfile = load('e:\ HFSS\Feed with phase v2\30 degrees\ Results\s1to1024-lossy.mat');
%         matfile = load('e:\ HFSS\Feed with phase v2\30 degrees stripline\ Results\s1to1024-lossy.mat');
        s1to1024 = matfile.s1to1024;
        fsS1x = matfile.fs;
        S1x = squeeze(s1to1024(1,:, :));
        S1xmat = reshape(S1x(2:end,:), 32, 32, length(fsS1x));
        
        excitation = zeros(Nx, Ny, Nf);
        % fs is from 12 to 32, S1xmat contains 13 to 31.
        excitation(:,:,2:end-1) = S1xmat(:,:,:);
        % Transpose the excitations for the stripline feeds.
%         excitation = permute(excitation, [2 1 3]);
%}


if(Nx > 32)
    excitation([(1:(Nx-32)/2) (end-(Nx-32)/2:end)], :, :) = 0;
end

% Add amplitude.
% Uniform amplitude
%     excitation = ones(Nx, Ny, length(fs)) .* excitation; 
% Random amplitude    
%     powerrange = [0.88 1.12]; % Determine effect of variation in amplitude in feeding network
%     dbrange = 2;
%     powerrange = [(2 - 10^(dbrange/20)) 10^(dbrange/20)];
%     excitation = (powerrange(1) + rand(Nx, Ny) * (powerrange(2) - powerrange(1))) .* excitation;
% Variation per column (along x)
%     dbdiff = 1;
%     factors = [10^((dbdiff/2)/20) (2 - 10^((dbdiff/2)/20))];
%     Nfac = length(factors);
%     for(i = 1:Nfac)
%         excitation(i:Nfac:end, :, :) = excitation(i:Nfac:end, :, :) .* factors(i);
%     end

% tlineup = TerminatedTLine(slab, FreeSpace());
% tlinedown = ShortedLine(erback*0.7, hback);

slot = Slot(dx, dy, wslot, dslot, walled);
%%
infarray = InfiniteArray(slot, tlineup, tlinedown);
Zinf = infarray.GetInputImpedance(fs, th, ph);
Sinf = (Zinf - zfeed) ./ (Zinf + zfeed);
% figureex; plot(fs/1e9, real(Zinf), fs/1e9, imag(Zinf));
% figureex; plot(fs/1e9, 20*log10(abs(Sinf)));
%%
% arrayx = FiniteArrayX(slot, tlineup, tlinedown, Nx, excitation(:,1).', dedge, zfeed);
% Zas = arrayx.GetInputImpedance(fs, th, ph);
% figureex; plot(fs, real(Zas), fs, imag(Zas));
% 
% arrayy = FiniteArrayY(slot, tlineup, tlinedown, Ny, excitation(1,:), zfeed);
% Zas = arrayy.GetInputImpedance(fs, th, ph);
% figureex; plot(fs, real(Zas), fs, imag(Zas));

newarray = 0;
if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny || array.unitcell.walled ~= walled)
    dispex('Creating array variable.\n');
    array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
    newarray = 1;
end

%% Determine path to place files.
filename = sprintf('%ix%i', Nx, Ny);
if(walled)
    filename = [filename, '_walled'];
end
matpath = sprintf('%s/Validations/%s/%s', resultdir_matlab, mfilename, filename);
% matpath = sprintf('E:/data/Sander/%s', filename);
% Ensure MAT file exists.
if(~exist(sprintf('%s.mat', matpath), 'file'))
    tc = tic;
    dispex('Creating .MAT file... ');
    save(sprintf('%s.mat', matpath), 'array', '-v7.3');
    fprintf('Done (%s).\n', fancyduration(toc(tc)));
end

%% Load variable if available.
if(newarray)
    tc = tic;
    dispex('Checking MAT contents... ');
    matContents = who('-file', sprintf('%s.mat', matpath));
    fprintf('Done (%s).\n', fancyduration(toc(tc)));
    if(ismember('array', matContents))
        tc = tic;
        dispex('Loading array variable... ');
        arrayM = load(sprintf('%s.mat', matpath), 'array');
        array = arrayM.array;
        fprintf('Done (%s).\n', fancyduration(toc(tc)));
    end
end

%% Run calculations and save results.
tc = tic;
%{
% Frequencies one-by-one.
lasttoc = 0;
for(fi = 1:length(fs))
    dispex('fi = %i/%i.\n', fi, length(fs))
    f = fs(fi);
    
    array.InitializeDs(f);
    array.InitializeZMatrix(f);
    
    if(toc(tc) > lasttoc+30)
        dispex('Saving array variable... ');
        save(matpath, 'array', '-v7.3');
        fprintf('Done.\n');
        lasttoc = toc(tc);
    end
end
%}
% Frequencies all at once.
array.InitializeDs(fs);
array.InitializeZMatrix(fs);
if(toc(tc) > 30)
    tc2 = tic;
    dispex('Saving array variable... ');
    save(matpath, 'array', '-v7.3');
    fprintf('Done (%s).\n', fancyduration(toc(tc2)));
end

%% Determine active input impedance.
Zact = array.GetInputImpedance(fs, excitation);

Zcap = 1 ./ (1j .* 2.*pi.*fs .* C);
for(fi = 1:length(fs))
    ZactC(:,:,fi) = Zact(:,:,fi) + Zcap(fi); % With series capacitance.
end

Sact = (ZactC - zfeed) ./ (ZactC+ zfeed);

dt_ML = toc(tc);
dispex('MATLAB took %.1fs for %ix%i elements.\n', dt_ML, Nx, Ny);
% save(matpath, 'dt_ML', '-append');

%% All-in-one Z
%{
dispex('Generating input impedance figure.\n');
[hFig, hAx] = figureex;
if(length(hAx.Children) < 2)
    patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none');
    patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
for(nx = 1:Nx)
    for(ny = 1:Ny)
        p = plot(hAx, fs/1e9, squeeze(real(Zas(nx, ny, :))), '-', 'Color', [0 0 0 0.05], 'LineWidth', 0.1);
        p = plot(hAx, fs/1e9, squeeze(imag(Zas(nx, ny, :))), ':', 'Color', [0 0 0 0.05], 'LineWidth', 0.1);
    end
end
plot(hAx, fs/1e9, imag(Zinf), 'r--', 'LineWidth', 3);
plot(hAx, [nan nan], 'k--');
plot(hAx, [nan nan], 'k');
plot(hAx, fs/1e9, real(Zinf), 'r', 'LineWidth', 3);
plot(hAx, [nan nan], 'k');
legend(hAx, hAx.Children(1:4), {'Finite', 'Infinite', 'Real', 'Imaginary'});
xlabel(hAx, 'Frequency [GHz]');
ylabel(hAx, 'Active Input Impedance [dB]');
xlim(hAx, [-inf inf]);
ylim(hAx, [-100 250]);
%}
%% All-in-one S
dispex('Generating reflection coefficient figure.\n');
[hFig, hAx] = figureex;
hFig.Name = 'Active Reflection Coefficient';
if(length(hAx.Children) < 2)
    patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none');
    patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
for(nx = 1:Nx)
    for(ny = 1:Ny)
        p = plot(hAx, fs/1e9, squeeze(20*log10(abs(Sact(nx, ny, :)))), 'k', 'Color', [0 0 0 0.1], 'LineWidth', 0.1);
    end
end
plot(hAx, fs/1e9, 20*log10(abs(Sinf)), 'r', 'LineWidth', 3);
plot(hAx, [nan nan], 'k');
legend(hAx, hAx.Children(1:2), {'Finite', 'Infinite'});
xlabel(hAx, 'Frequency [GHz]');
ylabel(hAx, 'Active Reflection Coefficient [dB]');
xlim(hAx, [-inf inf]);
ylim(hAx, [-30 0]);
%% Grouped S
if(~isempty(groups))
    dispex('Generating grouped coefficient figure.\n');
    % Determine S of the grouped elements.
    groupedelements = zeros(Nx, Ny);
    Zgroups = ZactC;
    Sgroups = Sact;
    for(nx = 1:Nx)
        for(ny = 1:Ny)
            
            
            [ind, ~] = find(groups(:, [1 3]) == nx & groups(:, [2 4]) == ny);
            if(~isempty(ind))
                group = groups(ind, :);
                ng = length(group)/2;
                if(Zgroups(group(1), group(2), 1) == ZactC(group(1), group(2), 1))
                    Ztot = zeros(1, 1, length(fs));
                    for(ig = 1:2:length(group))
                        Ztot = Ztot + ZactC(group(ig), group(ig+1), :);
                    end
                    for(ig = 1:2:length(group))
                        Zgroups(group(ig), group(ig+1), :) = Ztot;
                        Sgroups(group(ig), group(ig+1), :) = (Ztot - ng*zfeed) ./ (Ztot + ng*zfeed);
                        groupedelements(group(ig), group(ig+1)) = ind;
                    end
                end
            end
        end
    end
    
    [hFig, hAx] = figureex;
    [hFig2, hAx2] = figureex;
    if(length(hAx.Children) < 2)
        patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    if(length(hAx2.Children) < 2)
        patch(hAx2, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
        patch(hAx2, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
            'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    for(nx = 1:Nx)
        for(ny = 1:Ny)
            if(groupedelements(nx, ny) ~= 0)
                p = plot(hAx, fs/1e9, squeeze(20*log10(abs(Sgroups(nx, ny, :)))), 'k', 'Color', [0 0 0 0.1], 'LineWidth', 0.1);
                p = plot(hAx2, fs/1e9, squeeze(20*log10(abs(Sact(nx, ny, :)))), 'k', 'Color', [0 0 0 0.1], 'LineWidth', 0.1);
            end
        end
    end
    plot(hAx, fs/1e9, 20*log10(abs(Sinf)), 'r', 'LineWidth', 3);
    plot(hAx, [nan nan], 'k');
    legend(hAx, hAx.Children(1:2), {'Finite', 'Infinite'});
    xlabel(hAx, 'Frequency [GHz]');
    ylabel(hAx, 'Active Reflection Coefficient [dB]');
    xlim(hAx, [-inf inf]);
    ylim(hAx, [-30 0]);
    
    plot(hAx2, fs/1e9, 20*log10(abs(Sinf)), 'r', 'LineWidth', 3);
    plot(hAx2, [nan nan], 'k');
    legend(hAx2, hAx.Children(1:2), {'Finite', 'Infinite'});
    xlabel(hAx2, 'Frequency [GHz]');
    ylabel(hAx2, 'Active Reflection Coefficient [dB]');
    xlim(hAx2, [-inf inf]);
    ylim(hAx2, [-30 0]);
end
%% Input impedance per row
%{
[hFig, hAx] = figureex;
for(nx = 1:Nx)
    hAx = subplot(6,6,nx);
    hold(hAx, 'on');
    grid(hAx, 'on');
    title(hAx, sprintf('nx = %i', nx));
    xlim(hAx, [-inf, inf]);
    ylim(hAx, [-100 150]);
    hAx.Colormap = jet(Ny);
    repeatcolormap(hAx, 2);
    for(ny = 1:Ny)
        plot(hAx, fs/1e9, squeeze(real(Zas(nx, ny, :))));
        plot(hAx, fs/1e9, squeeze(imag(Zas(nx, ny, :))), '--');
    end
end
%}
%% Reflection coefficient
%{
numplotsx = ceil(sqrt(ceil(Nx/2)));
numplotsy = ceil(ceil(Nx/2)/numplotsx);
[hFig, hAx] = figureex;
for(nx = 1:ceil(Nx/2))
    hAx = subplot(numplotsx, numplotsy,nx);
    hold(hAx, 'on');
    grid(hAx, 'on');
    title(hAx, sprintf('nx = %i', nx));
    xlim(hAx, [-inf inf]);
    ylim(hAx, [-30 0]);
    hAx.Colormap = jet(Ny);
    for(ny = 1:Ny)
        plot(hAx, fs/1e9, squeeze(20*log10(abs(Sact(nx, ny, :)))));
    end
end
%}
%% Efficiency
dispex('Generating efficiency figure.\n');
Pin = squeeze(sum(sum(abs(excitation).^2)));
Prefl = squeeze(sum(sum(abs(Sact).^2.*abs(excitation).^2, 1), 2));
Preflinf = abs(Sinf).^2;
Eff = (Pin-Prefl)./(Pin);
Effinf = (1-Preflinf) ./1;

[hFig, hAx] = figureex;
hFig.Name = 'Efficiency';
if(length(hAx.Children) < 2)
    patch(hAx, [28 31 31 28], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none');
    patch(hAx, [13.75 14.5 14.5 13.75], [-1e3 -1e3 1e3 1e3], [0 0 0], ...
        'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
plot(hAx, fs/1e9, Effinf*100, 'r', 'LineWidth', 1);
plot(hAx, fs/1e9, Eff*100, 'k', 'LineWidth', 1);
hAx.LineWidth = 1;
legend(hAx, hAx.Children(1:2), {'Finite', 'Infinite'});
xlabel(hAx, 'Frequency [GHz]');
ylabel(hAx, 'Efficiency [%]');
xlim(hAx, [-inf inf]);
ylim(hAx, [0 100]);
%% Efficiency legend
%{
legend(hAx, 'off');
plot(hAx, [nan nan], 'r');
addlegendentry(hAx, 'Infinite');
plot(hAx, [nan nan], 'k');
addlegendentry(hAx, 'Finite');
plot(hAx, [nan nan], 'k');
addlegendentry(hAx, 'Broadside');
plot(hAx, [nan nan], 'k--');
addlegendentry(hAx, 'H-plane');
plot(hAx, [nan nan], 'k-.');
addlegendentry(hAx, 'E-plane');

%}
%% Reflection coefficient colormap
for(f = [14 31]*1e9)
    fi = find(fs == f);
    if(isempty(fi))
        error('S not calculated for f = %g', f);
    end
    [hFig, hAx] = figureex;
        delete(hAx);
        imagesc(0:Nx-1, 0:Ny-1, squeeze(20*log10(abs(Sact(:, :, fi)))).');
        hAx = hFig.CurrentAxes;
        hAx.YDir = 'normal';
        caxis(hAx, [-10 0])
        if(th < 1 * pi/180) % Broadside
            caxis(hAx, [-15 0])
        end
        xlabel(hAx, 'Nx');
        ylabel(hAx, 'Ny');
        hBar = colorbar(hAx);
        hBar.Label.String = 'S_{ii} [dB]';
        hBar.Label.FontSize = hAx.XLabel.FontSize;
        hFig.Name = sprintf('%.1f GHz', f/1e9);
        axis(hAx, 'square');
        alignplot(hFig, 8, 4, [], hFig.Number, 1);
        colormap(hAx, 'jet');
end


%%
for(iFig = 1:100)
    if(ishandle(iFig))
        hFig = figure(iFig);
        alignplot(hFig, 8, 4, [], hFig.Number, 1);
    end
end
%%
% [hFig4, hAx4] = figureex;
%     plot(hAx4, fs/1e9, squeeze(real(Zas(5,5, :))), 'k');
%     plot(hAx4, fs/1e9, squeeze(imag(Zas(5,5,:))), 'r');
    
% [hFig5, hAx5] = figureex;
% for(nx = 1:Nx)
%     for(ny = 1:Ny)
%         plot(hAx5, fs/1e9, squeeze(20*log10(abs(Sact(nx,ny,:)))), 'r');
%     end
% end


% savefig(hFig, sprintf('H:\\Git\\PhD-Matlab\\Tests\\%s\\figures\\%ix%i', mfilename, Nx, Ny));

% clear hFig hAx;