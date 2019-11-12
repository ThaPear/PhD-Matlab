function [h, s] = DrawADS(ads, f0, h, s, hFig)
    if(nargin < 2)
        disp('f0 not given in DrawADS');
        return;
    end
    metalheight = 3e-5;
    metalhsv = [160 0 120] ./ 255;
    metalclr = hsv2rgb(metalhsv);
    lineclr = [50 50 50] ./ 255;
    dielectrichsv = [135 92 255] ./ 255;
    dielectricclr = hsv2rgb(dielectrichsv);
    
    if(nargin < 5)
        % Open figure window.
        hFig = figureex;
            hold on;
            ax = hFig.CurrentAxes;
            ax.Clipping = 'off';
            view(60, 20);
            hFig.Name = ['Slab at ', num2str(f0/1e9), ' GHz'];
        % Draw top and bottom dielectric.
        p = GetP(ads);
        adsH = ads.GetHeight() / p;
%         box = DrawBox(-1.01, -1.01, 0, 2.02, 2.02, adsH, dielectricclr);
%         set(box, 'FaceLighting', 'none');
%         alpha(box, 0);

        xs = [[-1 1 1 -1 -1] [-1 1 1 -1 -1] [-1 1 1 -1 -1] [-1 1 1 -1 -1] repmat(1, 1, 5) repmat(-1, 1, 5)];
        ys = [[1 1 -1 -1 1] [1 1 -1 -1 1] repmat(1, 1, 5) repmat(-1, 1, 5) [-1 1 1 -1 -1] [-1 1 1 -1 -1]];
        zs = [repmat(0, 1, 5) repmat(adsH, 1, 5) [0 0 adsH adsH 0] [0 0 adsH adsH 0] [0 0 adsH adsH 0] [0 0 adsH adsH 0]];
        for(i = 1:length(xs)-1)
            plot3(xs(i:i+1), ys(i:i+1), zs(i:i+1), 'k-');
        end
%         plot3([-1 1 1 -1 -1], [1 1 -1 -1 1], repmat(0, 1, 5), 'k-');
%         plot3([-1 1 1 -1 -1], [1 1 -1 -1 1], repmat(adsH, 1, 5), 'k-');
%         plot3([-1 1 1 -1 -1], repmat(1, 1, 5), [0 0 adsH adsH 0], 'k-');
%         plot3([-1 1 1 -1 -1], repmat(-1, 1, 5), [0 0 adsH adsH 0], 'k-');
%         plot3(repmat(1, 1, 5), [-1 1 1 -1 -1], [0 0 adsH adsH 0], 'k-');
%         plot3(repmat(-1, 1, 5), [-1 1 1 -1 -1], [0 0 adsH adsH 0], 'k-');
        
        lgt = light('Position', [-4, 4, 8], 'Style', 'local');
        
        h = 0;
        s = 0;
    end
    % If it's not an ads, it's a TLine with ADSs inside. (e.g. TaperedADS)
    if(~isa(ads, 'ADS'))
        % Draw all elements of this line.
        for(i = 1:length(ads.elements))
            [h, s] = DrawADS(ads.elements{i}, f0, h, s, hFig);
        end
        return;
    end
    disp('Drawing ADS');
    
    p = ads.p;
    % Normalize to p.
    ds = ads.ds ./ p;
    ss = ads.ss ./ p;
    ws = ads.ws ./ p;
    if(isempty(ads.erhosts))
        erhost = 1;
    else
        erhost = ads.erhosts(1); % TODO: Fix for varying erhost per layer.
    end
    
    hmetal = metalheight ./ p; % 1 micron

    N = length(ws);
    
    % Draw the dielectric host material.
%     adsH = ads.GetHeight() ./ p;
%     dielectricclr = hsv2rgb([dielectrichsv(1), dielectrichsv(2), dielectrichsv(3) ./ max(1, erhost * 0.9)]);
%     dielectricBox = DrawBox(-1.01, -1.01, h, 2.02, 2.02, adsH, dielectricclr, 'none');
%     alpha(dielectricBox, 0.5);
%     set(dielectricBox, 'FaceLighting', 'none');
%     if(erhost == 1)
%         delete(dielectricBox([1 4])); % Remove top & bottom dielectrics.
%     end
    
%     % Draw the parameter list.
%     plot3([1 1.2 1.2 1], [1 1.2 1.2 1], [h h h+adsH h+adsH], 'k');
    plot3([-1 1 1 -1 -1], [1 1 -1 -1 1], repmat(h, 1, 5), 'k-.');
%     plot3([-1 1 1], [1 1 -1], repmat(h, 1, 3), 'k-.');
% %     plot3([-1 1 1 -1 -1], [1 1 -1 -1 1], repmat(h+adsH, 1, 5), 'k:', 'LineWidth', 1.5);
%     paramstr = ads.PrintParameters(f0);
%     text(1.3, 1.3, h+adsH/2, strsplit(paramstr, ' - '), 'FontSize', 8);
    
    h = h + ds(1);
    for(n = 1:N)
        disp(['Layer ', num2str(n), '/' num2str(N)]);
        s = mod(s+1/2, 1)-1/2; % Ensure -1/2 < s < 1/2
        if(ws(n) < 0)
            breakpoint;
        end
        % |-----------|   |-----------|   |-----------|
        % | top left  |   |    top    |   | top right |
        % |-----------|   |-----------|   |-----------|
        %
        % |-----------|   |-----------|   |-----------|
        % | mid left  |   |    mid    |   | mid right |
        % |-----------|   |-----------|   |-----------|
        %
        % |-----------|   |-----------|   |-----------|
        % | bot left  |   |    bot    |   | bot right |
        % |-----------|   |-----------|   |-----------|
        
        DrawBox(-1/2+ws(n)/2+s, -1/2+ws(n)/2+s, h, 1-ws(n), 1-ws(n), hmetal, metalclr, lineclr); % mid
        DrawBox(-1/2+ws(n)/2+s,  1/2+ws(n)/2+s, h, 1-ws(n), 1-ws(n), hmetal, metalclr, lineclr); % top
        DrawBox( 1/2+ws(n)/2+s,  1/2+ws(n)/2+s, h, 1-ws(n), 1-ws(n), hmetal, metalclr, lineclr); % top right
        DrawBox( 1/2+ws(n)/2+s, -1/2+ws(n)/2+s, h, 1-ws(n), 1-ws(n), hmetal, metalclr, lineclr); % right
        DrawBox( 1/2+ws(n)/2+s, -3/2+ws(n)/2+s, h, 1-ws(n), 1-ws(n), hmetal, metalclr, lineclr); % bottom right
        DrawBox(-1/2+ws(n)/2+s, -3/2+ws(n)/2+s, h, 1-ws(n), 1-ws(n), hmetal, metalclr, lineclr); % bottom
        DrawBox(-3/2+ws(n)/2+s, -3/2+ws(n)/2+s, h, 1-ws(n), 1-ws(n), hmetal, metalclr, lineclr); % bottom left
        DrawBox(-3/2+ws(n)/2+s, -1/2+ws(n)/2+s, h, 1-ws(n), 1-ws(n), hmetal, metalclr, lineclr); % left
        DrawBox(-3/2+ws(n)/2+s,  1/2+ws(n)/2+s, h, 1-ws(n), 1-ws(n), hmetal, metalclr, lineclr); % top left
        
        if(n < N)
            s = s + ss(n);
            h = h + ds(n+1);
        end
    end
    h = h + ds(end);
    s = s + ss(end);
    
    xlim([-1/2 1/2]*2);
    ylim([-1/2 1/2]*2);
    zlim([0 h]);
%     axis equal;
    axis equal;
    axis off;
end

function p = GetP(slab)
    if(isprop(slab, 'p'))
        p = slab.p;
    else
        p = GetP(slab.elements{1});
    end
end

function box = DrawBox(x, y, z, dx, dy, dz, clr, lineclr)
    %          2 ------ 3
    %         /|       /|
    %        / 6 -----/ 7
    %       1 ------ 4 /
    %       |/       |/
    %       5 ------ 8
    if(nargin < 8)
        lineclr = [0 0 0];
    end
    box = [];
    if(   x > 1 ...
       || x + dx < -1 ...
       || y > 1 ...
       || y + dy < -1)
        return;
    end
    
    if(dx < 0 || dy < 0)
        breakpoint;
    end
    
    % Adjust parameters so that it fits within [-1 1]
    if(x < -1)
        delta = -0.99-x;
        x = x + delta;
        dx = dx - delta;
    end
    if(y < -1)
        delta = -0.99-y;
        y = y + delta;
        dy = dy - delta;
    end
    dx = min(dx, 0.99-x);
    dy = min(dy, 0.99-y);
    
    ps = [                ...
        x, y, z+dz;        ...
        x, y+dy, z+dz;     ...
        x+dx, y+dy, z+dz;  ...
        x+dx, y, z+dz;     ...
        x, y, z;          ...
        x, y+dy, z;       ...
        x+dx, y+dy, z;    ...
        x+dx, y, z        ...
    ];
    
    planes = [      ...
        1, 2, 3, 4; ... % top
%         1, 2, 6, 5; ... % left
        1, 4, 8, 5; ... % front
%         5, 6, 7, 8; ... % bottom
%         2, 3, 7, 6; ... % back
        3, 4, 8, 7; ... % right
    ];

    for(i = 1:size(planes, 1))
        box(i) = DrawPlane(ps(planes(i,1), :), ps(planes(i,2), :), ps(planes(i,3), :), ps(planes(i,4), :), clr, lineclr);
    end
end

function plane = DrawPlane(p1, p2, p3, p4, clr, lineclr)
    if(nargin < 6)
        lineclr = [0 0 0];
    end
    if(nargin < 5)
        clr = rand(1, 3);
    end
    x = [p1(1) p2(1) p3(1) p4(1)];
    y = [p1(2) p2(2) p3(2) p4(2)];
    z = [p1(3) p2(3) p3(3) p4(3)];
    
    plane = fill3(x, y, z, clr, 'EdgeColor', lineclr);
    plane.FaceLighting = 'gouraud';
    plane.BackFaceLighting = 'lit';
end