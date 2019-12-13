function isopen = isopenineditor(filename)
    path = which(filename);
    if(~isempty(path))
        isopen = matlab.desktop.editor.isOpen(path);
    else
        isopen = 0;
    end
end