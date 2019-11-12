function settodb(hFig)
    if(~ishandle(hFig))
        warning('Invalid handle supplied.');
        return;
    elseif(~strcmp(hFig.Type, 'figure'))
        warning(['''', hFig.Type, ''' can not set  to db.']);
        return;
    end
    toolbars = findall(hFig,'type','uitoolbar');
    if(~isempty(toolbars))
        toolbars(1).Children(end).ClickedCallback(hFig, [])
    end
end