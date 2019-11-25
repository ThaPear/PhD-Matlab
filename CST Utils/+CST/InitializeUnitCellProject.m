function project = InitializeUnitCellProject()
    [project, dsproject] = CST.InitializeBasicProject();

    %% Set up boundaries.
    boundary = project.Boundary();
    boundary.AllBoundaries('unit cell', 'unit cell', ...   % x
                           'unit cell', 'unit cell', ...   % y
                           'expanded open', 'expanded open'); % z

    %% De-embed the floquet ports.
    floquetport = project.FloquetPort();
    floquetport.StartBulkMode();
    floquetport.Port('Zmin');
    floquetport.SetDistanceToReferencePlane('-openboundary_distance');
    floquetport.SetNumberOfModesConsidered(2);
    floquetport.Port('Zmax');
    floquetport.SetDistanceToReferencePlane('-openboundary_distance');
    floquetport.SetNumberOfModesConsidered(2);
    floquetport.EndBulkMode();
end