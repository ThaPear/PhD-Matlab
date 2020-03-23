function cropfigure(hFig, padding)
    % hFig: Figure to crop.
    % padding: Space to leave around the figure contents. [x, y]
    if(nargin < 1)
        hFig = gcf;
    end
    if(nargin < 2)
        padding = [0, 0];
    end
    padding = padding + [1 1];
    % Get figure color data.
    frame = getframe(hFig);
    cdata = frame.cdata;
    cdata = sum(cdata, 3);
    cdata = flipud(cdata);
    
    % Default background color, assume first value is empty.
    bgclr = cdata(1,1);
    
    % Figure size.
    szx = hFig.Position(3);
    szy = hFig.Position(4);
    
    % Detect the first column which has anything other than the default
    % background color.
    cdatarow = sum(cdata, 1);
    firstx = find(cdatarow ~= bgclr * szy, 1) - 1 - 1;
    lastx  = find(cdatarow ~= bgclr * szy, 1, 'last') + 1;
    width = lastx - firstx;
    
    % Detect the first row which has anything other than the default
    % background color.
    cdatacol = sum(cdata, 2);
%     figure; plot(cdatacol/szx);
    firsty = find(cdatacol ~= bgclr * szx, 1) - 1;
    lasty  = find(cdatacol ~= bgclr * szx, 1, 'last') + 1;
    height = lasty - firsty;
    
    elements = hFig.Children;
    for(i = 1:length(elements))
        el = elements(i);
        % If it has the 'Location' property it'll move along with the axes.
        if(isprop(el, 'Location') || ~isprop(el, 'Units'))
            continue;
        end
        el.Units = 'pixels';
        el.Position(1:2) = el.Position(1:2) - [firstx, firsty] + padding;
    end
    
    hFig.Position(3:4) = [lastx - firstx, lasty - firsty] + 2.*padding;
end