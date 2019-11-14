function XPolCancellation(project, port1, port2, frequency, theta, phi)
    %%
    if(frequency >= 1e9)
        frequency = frequency / 1e9;
    end

    plot1d = project.Plot1D();

    plot1d.DeleteAllMarker();
    plot1d.XMarker(1);
    plot1d.XMarkerPos(60);
%     plot1d.ShowMarkerAtMax();
    plot1d.Plot();

    resultname = sprintf('Farfields\\farfield (f=%i) [%i]\\Ludwig 3 Horizontal', frequency, port1);
    project.SelectTreeItem(resultname);
    p1_hor = inputdlgex('Enter value shown in CST.');
    if(isempty(p1_hor)); return; end

    resultname = sprintf('Farfields\\farfield (f=%i) [%i]\\Ludwig 3 Hor. Phase', frequency, port1);
    project.SelectTreeItem(resultname);
    p1_horphase = inputdlgex('Enter value shown in CST.');
    if(isempty(p1_horphase)); return; end

    resultname = sprintf('Farfields\\farfield (f=%i) [%i]\\Ludwig 3 Horizontal', frequency, port2);
    project.SelectTreeItem(resultname);
    p2_hor = inputdlgex('Enter value shown in CST.');
    if(isempty(p2_hor)); return; end

    resultname = sprintf('Farfields\\farfield (f=%i) [%i]\\Ludwig 3 Hor. Phase', frequency, port2);
    project.SelectTreeItem(resultname);
    p2_horphase = inputdlgex('Enter value shown in CST.');
    if(isempty(p2_horphase)); return; end

    p1_hor = str2double(p1_hor{:});
    p1_horphase = str2double(p1_horphase{:});
    p2_hor = str2double(p2_hor{:});
    p2_horphase = str2double(p2_horphase{:});

    fprintf('Amplitude: %.15g / %.15g = \n%.15g\nPhase: 180 - (%.15g - %.15g) = \n%.15g\n', p2_hor, p1_hor, p2_hor / p1_hor, p2_horphase, p1_horphase, mod(180-(p2_horphase - p1_horphase), 360));
end