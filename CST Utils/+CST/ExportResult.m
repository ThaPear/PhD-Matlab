function ExportResult(project, resultname, exportfilename)
    plot1d = project.Plot1D();
    asciiexport = project.ASCIIExport();

    project.SelectTreeItem(resultname);

    plot1d.PlotView('real');
    asciiexport.Reset();
    asciiexport.FileName([exportfilename, '.real']);
    asciiexport.Execute();
    plot1d.PlotView('imaginary');
    asciiexport.Reset();
    asciiexport.FileName([exportfilename, '.imag']);
    asciiexport.Execute();
end