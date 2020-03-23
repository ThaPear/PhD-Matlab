function addlegendentry(varargin)
% Adds entries to the legend.
% ADDLEGENDENTRY(entries) adds N entries for the last N lines.
%
% ADDLEGENDENTRY(hLines, entries) adds N entries for the given N lines.
%
% ADDLEGENDENTRY(hAx, hLines, entries) adds N entries for the given N lines
% in the given axes.

    switch(nargin)
        case 1
            arg1 = varargin{1};

            entries = arg1;
            if(~iscell(entries))
                entries = {entries};
            end

            hAx = gca;
            hLines = hAx.Children(1:length(entries));
        case 2
            arg1 = varargin{1};
            arg2 = varargin{2};

            entries = arg2;
            if(~iscell(entries))
                entries = {entries};
            end

            if(ishandle(arg1) && strcmp(arg1.Type, 'axes'))
                hAx = arg1;
                hLines = hAx.Children(1:length(entries));
            elseif(ishandle(arg1) && strcmp(arg1.Type, 'figure'))
                hAx = arg1.CurrentAxes;
                hLines = hAx.Children(1:length(entries));
            else
                hAx = gca;
                hLines = arg1;
            end
        case 3
            arg1 = varargin{1};
            arg2 = varargin{2};
            arg3 = varargin{3};

            if(strcmp(arg1.Type, 'axes'))
                hAx = arg1;
            elseif(strcmp(arg1.Type, 'figure'))
                hAx = arg1.CurrentAxes;
            else
                error('Invalid first argument type. Should be axes or figure.');
            end
            hLines = arg2;
            entries = arg3;
            if(~iscell(entries))
                entries = {entries};
            end
    end
    
    leg = hAx.Legend;
    if(isempty(leg))
        % Don't have to add since there's no legend yet.
        legend(hAx, hLines, entries);
        leg = hAx.Legend;
        leg.AutoUpdate = 'off';
        leg.UserData.hLines = hLines;
        leg.UserData.entries = entries;
        return;
    elseif(isempty(leg.UserData))
        warning('Cannot add legend entry to legend that was not created by addlegendentry.');
        return;
    end
    
    entries = [leg.UserData.entries, entries];
    hLines = [leg.UserData.hLines; hLines];
    
    % Remove entries that refer to invalid lines.
    entries = entries(ishandle(hLines));
    hLines = hLines(ishandle(hLines));
    
    legend(hAx, hLines, entries);
    
    leg.UserData.hLines = hLines;
    leg.UserData.entries = entries;
end
