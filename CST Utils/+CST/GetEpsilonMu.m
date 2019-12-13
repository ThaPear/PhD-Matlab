function [er, mu] = GetEpsilonMu(tline, f, Nlayers)
    project = CST.InitializeUnitCellProject();
    
    project.StoreParameter('dx', 4.4);
    project.StoreParameter('dy', 'dx');
    
    % Repeat the structure this number of times.
    if(nargin < 3)
        % Multilayer makes no sense as permittivity extraction only works for thin layers.
        Nlayers = 1;
    end
    
    wcs = project.WCS();
    wcs.Enable();
    for(i = 1:Nlayers)
        tline.BuildCST(project);
    end
    wcs.Disable();
    
    % 
    th0 = 0 * pi/180;
    th = 10 * pi/180;
    ph = 0 * pi/180;
    
    er = [];
    mu = [];
    
    %% Simulate only the given frequency.
    if(f >= 1e9)
        fCST = f / 1e9;
    else
        fCST = f;
    end
    fdsolver = project.FDSolver();
    fdsolver.ResetSampleIntervals('all');
    fdsolver.AddSampleInterval(fCST, '', 1, 'Single', 1);
    
    project.StoreParameter('fmesh', fCST);
    project.StoreParameter('fmin', fCST-1);
    project.StoreParameter('fmax', fCST+1);
    
    % If the ADL is in free space, de-embed the port to the correct distance.
    floquetport = project.FloquetPort();
    if(isprop(tline.elements{1}, 'er') && tline.elements{1}.er == 1)
        floquetport.StartBulkMode();
        floquetport.Port('Zmin');
        floquetport.SetDistanceToReferencePlane(['-openboundary_distance+', num2str(tline.elements{1}.L*1e3, '%.15g')]);
        floquetport.Port('Zmax');
        floquetport.SetDistanceToReferencePlane(['-openboundary_distance+', num2str(tline.elements{end}.L*1e3, '%.15g')]);
        floquetport.EndBulkMode();
    end
%     project.StoreParameter('openboundary_distance', [CST.Defaults.OpenBoundaryDistance, '+', num2str(tline.elements{end}.L*1e3, '%.15g')]);


    %% Simulate for broadside incidence.
    project.StoreParameter('aa_theta', th0 * 180/pi);
    if(~project.Rebuild())
        dispex('Project rebuild failed.\n');
        return;
    end
    if(~fdsolver.Start())
        dispex('Simulation failed.\n');
        return;
    end
    
    filename = 'Temp/S0';
    CST.ExportResult(project, '1D Results\S-Parameters', filename);
    [fs, parameters, Ste0, Stm0] = CST.LoadData([filename, '.s4p']);

    %% Simulate for 10-degree incidence.
    project.StoreParameter('aa_theta', th * 180/pi);
    if(~project.Rebuild())
        dispex('Project rebuild failed.\n');
        return;
    end
    if(~fdsolver.Start())
        fprintf('Simulation failed.\n');
        return;
    end
    
    filename = 'Temp/S';
    CST.ExportResult(project, '1D Results\S-Parameters', filename);
    [fs, parameters, Ste, Stm] = CST.LoadData([filename, '.s4p']);
    
    d = Nlayers * tline.GetHeight();
    [er, mu] = EpsilonMu(f, Ste, Stm, Ste0, Stm0, d, th0, th, ph);
    
%     project.ResetAll();
    project.Quit();
end




























