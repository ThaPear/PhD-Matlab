close all;
clear;
SetupPath;



params = struct('wslot', 0.8, ...
                'dslot', 0.3, ...
                'wfeed', 0.3, ...
                'th', 0, ...
                'ph', 0, ...
                'coax', 1, ...
                'scor', 1, ...
                'Ap', 1.2);
            
            
folder = sprintf('%s/ Reports/WP_2100&2400/gb_coax/', resultdir_matlab);

ths = [eps 60 60] * pi/180;
phs = [0   0  45] * pi/180;
% ths = [eps 60] * pi/180;
% phs = [0   0 ] * pi/180;
th0 = 25;
f0 = 30;
fs_field = 9:33;

for(iangle = 1:length(ths))
    th = ths(iangle);
    ph = phs(iangle);
    th0 = round(th * 180/pi);
    
    %% Load the correct Touchstone file.
    params.th = th;
    params.ph = ph;
    
    filename = struct2string(params);
    filename = regexprep(filename, '[+*]', '_');

    filepath = [folder, filename, '.cst'];

    touchstonefile = [filepath(1:end-4), '.s2p'];
    [parameters, S] = CST.LoadData(touchstonefile);
    fs = parameters.frequencies;
    
    %% Load the fields
    path = [cd, '\PhD-Matlab\ Reports\WP_2100&2400\ra_circular\'];
    load([path, '0.127_ff_', num2str(th0), '_', num2str(ph*180/pi), '.mat']);
    
    EthX = Eth(1,:);
    EphX = Eph(1,:);
    EthY = Eth(2,:);
    EphY = Eph(2,:);

    EthX0 = interp1(fs_field, EthX, f0);
    EphX0 = interp1(fs_field, EphX, f0);
    EthY0 = interp1(fs_field, EthY, f0);
    EphY0 = interp1(fs_field, EphY, f0);
    
    %% Calculate weights
    % Slant Linear
    alphas = pi/180 * (0:15:90);
    alphas = pi/180 * -45;
    for(alpha = alphas)
    ay = 1;    ax = -ay .* (EphY0 .* cos(alpha) - EthY0 .* sin(alpha)) ./ (EphX0 .* cos(alpha) - EthX0 .* sin(alpha));
    if(abs(ax) > 1)
        ay = ay ./ ax;
        ax = ax ./ ax;
    end
    
    %% Combine S-params to get active values
    S11 = squeeze(S(1,1,:));
    S12 = squeeze(S(1,2,:));
    S21 = squeeze(S(2,1,:));
    S12 = S21;
    S21 = S12;
    S22 = squeeze(S(2,2,:));
    
    S11combined = (S11 .* ax + S12 .* ay) ./ ax;
    S22combined = (S22 .* ay + S21 .* ax) ./ ay;
    
    VSWR11 = S2VSWR(S11);
    VSWR22 = S2VSWR(S22);
    VSWR11combined = S2VSWR(S11combined);
    VSWR22combined = S2VSWR(S22combined);
    
    %% Plot results
    [hFig, hAx] = a_figVSWR(iangle);
%         plot(hAx, fs./1e9, VSWR11);
%         plot(hAx, fs./1e9, VSWR22);
        plot(hAx, fs./1e9, VSWR11combined);
        plot(hAx, fs./1e9, VSWR22combined);
        xlim(hAx, [12 32]);
        ylim(hAx, [1 8]);
        legend(hAx, hAx.Children(2:-1:1), {'X-slot', 'Y-slot'});
        movelegend(hAx, 'nw');
        
    Pin = abs(ax).^2 + abs(ay).^2;
    Prefl = abs(S11combined .* ax).^2 + abs(S22combined .* ay).^2;
    efficiency = (1-Prefl./Pin);
    [hFig, hAx] = a_fig(9);
%         plot([min(fs), max(fs)], [Pin, Pin]);
        plot(hAx, fs./1e9, 100*efficiency);
        alignplot(hFig, 8, 4, hFig.Number, [], 2);
        xlim(hAx, [12 32]);
        ylim(hAx, [0 100]);
        xlabel(hAx, 'Frequency [GHz]');
        ylabel(hAx, 'Efficiency [%]');
        
    end
        
end

legend(hAx, hAx.Children(3:-1:1), {'Broadside', 'E-plane 60\circ', 'D-plane 60\circ'});
drawnow; movelegend(hAx, 's');