function ExportResult(project, resultname, exportfilename)
    % NOTE: Do not include file extension in exportfilename.
    plot1d = project.Plot1D();
    asciiexport = project.ASCIIExport();
    
    % If the specified path is a relative path, make it a global path
    exportfilename = GetFullPath(exportfilename);
%     exportfilename = exportfilename(1:find(exportfilename == '.', 1, 'last'));
    if(contains(exportfilename, '.'))
        warning('Parameter exportfilename contains ''.'', do not specify extension.');
    end
    
    % Ensure that the correct slash is used in specifying the Navigation Tree address.
    if(contains(resultname, '/'))
        warning('Result names must be specified with ''\'', not with ''/''. Autocorrecting.');
        resultname = strrep(resultname, '/', '\');
    end

    success = project.SelectTreeItem(resultname);
    
    
    resulttree = project.ResultTree();
    if(~success)
        childname = resulttree.GetFirstChildName(resultname);
        project.SelectTreeItem(childname);
    end
    
%     res = resulttree.GetResultIDsFromTreeItem([resultname, '\SZmax(1),Zmax(1)']);
    
    resultsplit = strsplit(resultname, '\');
    
    if(strcmp(resultsplit{end}, 'S-Parameters') || strcmp(resultsplit{end-1}, 'S-Parameters'))
        % Loop over parametric results.
        ret = project.ResultNavigatorRequest('set selection', '0');
        th = project.RestoreParameter('aa_theta');
        ph = project.RestoreParameter('aa_phi');
        touchstone = project.Touchstone();
        touchstone.Reset();
        touchstone.FileName(exportfilename);
        touchstone.Impedance(120*pi);
        touchstone.FrequencyRange('Full');
        touchstone.Renormalize(0);
%         touchstone.UseARResults(0);
%         touchstone.SetNSamples(100);
        touchstone.Write();
        return;
    end

    plot1d.PlotView('real');
    asciiexport.Reset();
    asciiexport.FileName([exportfilename, '.real']);
    asciiexport.Execute();
    plot1d.PlotView('imaginary');
    asciiexport.Reset();
    asciiexport.FileName([exportfilename, '.imag']);
    asciiexport.Execute();
end