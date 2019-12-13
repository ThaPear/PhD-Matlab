function closeineditor(filename)
    % Get editor.
    edtSvc  = com.mathworks.mlservices.MLEditorServices();
    % Get open files.
    edtList = edtSvc.getEditorApplication().getOpenEditors().toArray();

    path = which(filename);
    if(~isempty(path))
        % Find the right file to close.
        for k=1:length(edtList)
            if(strcmpi(path, edtList(k).getLongName().toString()))
                edtList(k).close();
            end
        end
    end
end

