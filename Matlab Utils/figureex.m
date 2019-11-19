% ---- Simulation
% ---- Util
% ---- figureex.m

function hFig = figureex(iFig)
    if(nargin == 0 || isempty(iFig))
        hFig = figure;
    else
        hFig = figure(iFig);
        toolbars = findall(hFig,'type','uitoolbar');
        if(length(toolbars) == 2)
            return;
        end
    end
    %% Settings
    limstringsX = {'-inf inf', '0 20', '0 15', '1 10'};
    limstringsY = {'-inf inf', '0 1', '-30 0', '-40 0', '-200 200'};
    plotindices = {'1', '2', 'end-1', 'end', 'all'};
    
    
    %% Hide the figure while it's being built.
    hFig.Visible = 'off';
    grid on; box on; hold on;
%     toolbars = findall(hFig,'type','uitoolbar');
%     if(~isempty(toolbars))
%         toolbar = toolbars(1);
%     else
        toolbar = uitoolbar(hFig);
%     end
    %% Toggle dB scale
    icon = ones(16,16,3) .* 240/255;
    icon = insertText(icon, [-2 -1], 'dB', 'BoxOpacity', 0);
    icon(icon == 240/255) = nan;
    dBButton = uipushtool(toolbar, 'TooltipString', 'Convert y scale to dB');
    dBButton.ClickedCallback = @(~, ~) SetYScaleTodB(hFig);
    dBButton.CData = icon;
    icon(10, 1:16, :) = 0;
    undBButton = uipushtool(toolbar, 'TooltipString', 'Convert y scale to linear');
    undBButton.ClickedCallback = @(~, ~) SetYScaleToLinear(hFig);
    undBButton.CData = icon;
    
    %% Toggle grid
    icon(:) = 1;
    icon([6 11], 1:16, :) = 200/255;
    icon(1:16, [6 11], :) = 200/255;
    icon([1 16], 1:16, :) = 0;
    icon(1:16, [1 16], :) = 0;
    icon(icon == 240/255) = nan;
    gridButton = uipushtool(toolbar, 'TooltipString', 'Toggle grid', 'Separator', 'on');
    gridButton.ClickedCallback = @(~, ~) ToggleGrid(hFig);
    gridButton.CData = icon;
    
    %% Expand axes
    icon(:) = nan; % nan is transparent
    icon(1:16, 8:9, :) = 0;
    icon([3 14], 6:11, :) = 0;
    icon([2 15], 7:10, :) = 0;
    icon(8:9, 1:16, :) = 0;
    icon(6:11, [3 14], :) = 0;
    icon(7:10, [2 15], :) = 0;
    axesButton = uipushtool(toolbar, 'TooltipString', 'Expand axes');
    axesButton.ClickedCallback = @(~, ~) expandaxes(hFig);
    axesButton.CData = icon;
    
    %% Xlim
    icon(:) = 1;
    icon(11:12, 1:16, :) = 0;
    icon(9:14, [3 14], :) = 0;
    icon(10:13, [2 15], :) = 0;
    icon = insertText(icon, [1.5,-6], 'x', 'BoxOpacity', 0);
    icon(icon == 1) = nan;
%     axesButton = uipushtool(toolbar, 'TooltipString', 'Set X limits', 'Separator', 'on');
%     axesButton.ClickedCallback = @(~, ~) RequestLim(hFig, 'x');
%     axesButton.CData = icon;
    but = uisplittool(toolbar, 'TooltipString', 'Set X limits', 'Separator', 'on');
    but.CData = icon;
    but.ClickedCallback = @(~, ~) RequestLim(hFig, 'x');
    drawnow;
    jButton = get(but, 'JavaContainer');
    drawnow;
    jMenu = get(jButton, 'MenuComponent');
    for(i = 1:length(limstringsX))
        limstr = limstringsX{i};
        jActionItem = jMenu.add(limstr);
        jActionItem.setText(strrep(limstr, ' ', ' to '));
        set(jActionItem, 'ActionPerformedCallback', @(~, ~) RequestLim(hFig, 'x', limstr));
    end
    
    %% Ylim
    icon(:) = 1;
    icon(1:16, 11:12, :) = 0;
    icon([3 14], 9:14, :) = 0;
    icon([2 15], 10:13, :) = 0;
    icon = insertText(icon, [-2,-1], 'Y', 'BoxOpacity', 0);
    icon(icon == 1) = nan;
    but = uisplittool(toolbar, 'TooltipString', 'Set Y limits');
    but.CData = icon;
    but.ClickedCallback = @(~, ~) RequestLim(hFig, 'y');
    drawnow;
    jButton = get(but, 'JavaContainer');
    drawnow;
    jMenu = get(jButton, 'MenuComponent');
    for(i = 1:length(limstringsY))
        limstr = limstringsY{i};
        jActionItem = jMenu.add(limstr);
        jActionItem.setText(strrep(limstr, ' ', ' to '));
        set(jActionItem, 'ActionPerformedCallback', @(~, ~) RequestLim(hFig, 'y', limstr));
    end
    
    %% Move legend
    % Doesn't work properly due to the drawnows messing with imagesc.
    positions = {'sw', 'nw', 'ne', 'se', 'w', 'n', 'e', 's'};
    texts = {'Bottom Left', 'Top Left', 'Top Right', 'Bottom Right', 'Left', 'Top', 'Right', 'Bottom'};
    x0s = [1, 1, 11, 11, 1, 6, 11, 6];
    y0s = [9, 1, 1, 9, 5, 1, 5, 9];
    for(i = 1:length(positions))
        pos = positions{i};
        x0 = x0s(i); y0 = y0s(i);
        icon(:) = 1;
        icon([1 16], 1:16, :) = 0;
        icon(1:16, [1 16], :) = 0;
        icon([y0 y0+7], x0:x0+5, :) = 0;
        icon(y0:y0+7, [x0 x0+5], :) = 0;
        icon(icon == 240/255) = nan;
        legendButton = uipushtool(toolbar, 'TooltipString', ['Move legend to ', lower(texts{i})]);
        legendButton.ClickedCallback = @(~, ~) movelegend(pos, hFig.CurrentAxes);
        legendButton.CData = icon;
        if(i == 1 || i == 5)
            legendButton.Separator = 'on';
        end
    end
    %% Copy figure to clipboard
    icons = load('icons.mat');
    icon = icons.clipboard;
    clipboardButton = uipushtool(toolbar, 'TooltipString', 'Copy figure to clipboard', 'Separator', 'on');
    clipboardButton.ClickedCallback = @(~,~) CopyToClipboard(hFig);
%     clipboardButton.ClickedCallback = @(~,~) hgexport(hFig,'-clipboard');
%     clipboardButton.ClickedCallback = 'print -clipboard -dmeta';
    clipboardButton.CData = icon;
    
    %% Copy EMF to clipboard
    icon = icons.powerpnt;
    clipboardButton = uipushtool(toolbar, 'TooltipString', 'Copy figure to clipboard as emf');
%     clipboardButton.ClickedCallback = @(~,~) CopyToClipboard(hFig);
%     clipboardButton.ClickedCallback = @(~,~) hgexport(hFig,'-clipboard');
    clipboardButton.ClickedCallback = 'print -clipboard -dmeta';
    clipboardButton.CData = icon;
    %print -clipboard -dbitmap
    
    %% Remove plots
    icon(:) = 1;
    for(xy = 3:14)
        icon(xy, xy, :) = [0.8 0 0];
        icon(xy-1, xy, :) = [0.8 0 0];
        icon(xy+1, xy, :) = [0.8 0 0];
        icon(18-xy, xy, :) = [0.8 0 0];
        icon(17-xy, xy, :) = [0.8 0 0];
        icon(16-xy, xy, :) = [0.8 0 0];
    end
    icon = insertText(icon, [-1,-5], '1', 'BoxOpacity', 0, 'FontSize', 16);
    icon(icon == 1) = nan;
    but = uisplittool(toolbar, 'TooltipString', 'Remove plot', 'Separator', 'on');
    but.CData = icon;
    but.ClickedCallback = @(~, ~) RemoveLine(hFig, '1');
    drawnow;
    jButton = get(but, 'JavaContainer');
    drawnow;
    jMenu = get(jButton, 'MenuComponent');
    for(i = 1:length(plotindices))
        idxstr = plotindices{i};
        jActionItem = jMenu.add(idxstr);
        jActionItem.setText(idxstr);
        set(jActionItem, 'ActionPerformedCallback', @(~, ~) RemoveLine(hFig, idxstr));
    end
    
    %% Restore the pre-2018 figure toolbar.
    % From https://nl.mathworks.com/matlabcentral/answers/419036-what-happened-to-the-figure-toolbar-in-r2018b-why-is-it-an-axes-toolbar-how-can-i-put-the-buttons
    % To do it by default:
    % set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
    % set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
    ver = version('-release');
    if(str2double(ver(1:4)) >= 2018)
        addToolbarExplorationButtons(gcf); % Adds buttons to figure toolbar
        ax = gca;
        ax.Toolbar.Visible = 'off'; % Turns off the axes toolbar
    end
    %% Show the figure
    hFig.Visible = 'on';
end

function SetYScaleTodB(hFig)
    hAx = hFig.CurrentAxes;
%     hAx.YData = db(hAx.YData);
    %% Check if the dB string is already present.
    k = strfind(hAx.YLabel.String, '[dB]');
    if(~isempty(k))
        return;
    end
    
    children = hAx.Children;
    for(i = 1:length(children))
        child = children(i);
        if(isfield(child, 'YData') || isprop(child, 'YData'))
            if(min(child.YData) < 0)
                warning('Cannot change y scale to dB, there are negative values.');
                return;
            end
        end
    end
    %% Change y data of each line to dB.
    for(i = 1:length(children))
        child = children(i);
        if(isfield(child, 'YData') || isprop(child, 'YData'))
            child.YData = 10*log10(child.YData);
        end
    end
    % Insert dB string.
    hAx.YLabel.String = [hAx.YLabel.String, '[dB]'];
end

function SetYScaleToLinear(hFig)
    hAx = hFig.CurrentAxes;
%     hAx.YData = db(hAx.YData);

    %% Remove dB string from y label, or if it's not present, return.
    k = strfind(hAx.YLabel.String, '[dB]');
    if(isempty(k))
        return;
    end
    k = k(1);
    hAx.YLabel.String = [hAx.YLabel.String(1:k-1), hAx.YLabel.String(k+5:end)];
    
    %% Change y data of each line to linear.
    children = hAx.Children;
    for(i = 1:length(children))
        child = children(i);
%         if(isa(child, 'matlab.graphics.chart.primitive.Line'))
        if(isfield(child, 'YData') || isprop(child, 'YData'))
            child.YData = 10.^(child.YData / 10);
        end
    end
end

function ToggleGrid(hFig)
    hAx = hFig.CurrentAxes;
    names = get(hAx,'DimensionNames');
    xgrid = [names{1} 'Grid'];
    xminorgrid = [names{1} 'MinorGrid'];
    if(strcmp(get(hAx, xminorgrid), 'on'))
        grid(hAx, 'off');
    elseif(strcmp(get(hAx, xgrid), 'on'))
        grid(hAx, 'minor');
    else
        grid(hAx, 'on');
    end
end

function RequestLim(hFig, lim, limstr)
    switch(lim)
        case 'x'
        case 'y'
        otherwise
            error('Wrong lim, use ''x'' or ''y''.');
    end
    
    if(nargin < 3)
        prompt = {['Enter space-separated ', lim, ' limits.']};
        title = 'Input';
        dims = [1];
        definput = {'-inf inf'};
        answerstr = inputdlgex(prompt,title,dims,definput);
        if(isempty(answerstr))
            return;
        end
    else
        answerstr = {limstr};
    end
    % Split into 2 strings.
    answer = strsplit(cell2mat(answerstr));
    if(length(answer) ~= 2)
        warning('Setting limits requires lower and upper value.');
        return;
    end
    % Extract lims, set to inf if invalid.
    lims = [0 0];
    for(i = 1:2)
        lims(i) = str2double(answer{i});
    end
    if(isempty(lims(1)) || isnan(lims(1)))
        warning(['First number in ''', cell2mat(answerstr), ''' is not a valid number.']);
        return;
    end
    if(isempty(lims(2)) || isnan(lims(2)))
        warning(['Second number in ''', cell2mat(answerstr), ''' is not a valid number.']);
        return;
    end
    
    if(lims(2) < lims(1))
        lims = lims([2 1]);
    end
    hAx = hFig.CurrentAxes;
    limfunc = str2func([lim, 'lim']);
    limfunc(hAx, lims);
end

function CopyToClipboard(hFig)
    % Store old color.
    clr = hFig.Color;
    % Make background white.
    hFig.Color = [1 1 1];
    drawnow;
    % Get image data.
    f = getframe(hFig);
    % Copy to clipboard.
    imclipboard('copy', f.cdata);
    % Restore background color.
    hFig.Color = clr;
    drawnow;
end

function RemoveLine(hFig, index)
    for(iax = 1:length(hFig.Children))
        ax = hFig.Children(iax);
        if(~isa(ax, 'matlab.graphics.axis.Axes'))
            continue;
        end
        if(~isempty(ax.Children))
            if(strcmp(index, 'all'))
                for(i = 1:length(ax.Children))
                    delete(ax.Children(1));
                end
            elseif(contains(index, 'end'))
                if(length(index) > 3)
                    idx = str2double(index(4:end));
                    delete(ax.Children(end+idx));
                else
                    delete(ax.Children(end));
                end
            else
                if(~isnumeric(index))
                    index = str2double(index);
                    if(isnan(index))
                        error('Invalid index %s', index);
                    end
                end
                delete(ax.Children(index));
            end
        end
    end
end