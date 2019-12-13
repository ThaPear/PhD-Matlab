function [project, dsproject] = InitializeBasicProject()
    project = CST.Application.NewMWS();
    dsproject = CST.Application.ActiveDS();

    project.AddToHistory('ChangeSolverType "HF Frequency Domain"');

    %% Set up units.
    units = project.Units();
    units.AllUnits('mm', 'GHz');

    project.StoreParameter('c0', Constants.c0*1e3); % mm

    %% Set up solver.
    fmin = CST.Defaults.FrequencyRange(1);
    fmax = CST.Defaults.FrequencyRange(2);
    fmesh = CST.Defaults.MeshFrequency;
    project.StoreParameter('fmin', fmin);
    project.StoreParameter('fmax', fmax);
    project.StoreParameter('fmesh', fmesh);
    project.StoreParameter('nsamplesperGHz', CST.Defaults.SamplesPerGHz);

    %% Frequency range and sampling
    solver = project.Solver();
    fdsolver = project.FDSolver();
    solver.FrequencyRange('fmin', 'fmax');
    fdsolver.StartBulkMode();
    fdsolver.ResetSampleIntervals('all');
    % Add a single mesh adaptation sample at fmesh.
    fdsolver.AddSampleInterval('fmesh', '', 1, 'Single', 1);
    fdsolver.AddSampleInterval('fmin', 'fmax', 'nsamplesperGHz * (fmax-fmin) + 1', 'Equidistant', 0);
    fdsolver.AddSampleInterval('', '', '', 'Automatic', 0);

    % PML on top.
    fdsolver.SetOpenBCTypeTet('PML');
    fdsolver.EndBulkMode();

    %% Set up mesh adaptation settings
    meshadaption3d = project.MeshAdaption3D();
    meshadaption3d.SetType('HighFrequencyTet');
    meshadaption3d.MinPasses(6);
    meshadaption3d.MaxPasses(50);

    %% Set up open boundary
    boundary = project.Boundary();
    boundary.AllBoundaries('open', 'open', ...   % x
                           'open', 'open', ...   % y
                           'open', 'open'); % z
%     boundary.MinimumDistanceType('Fraction');
%     boundary.MinimumDistancePerWavelengthNewMeshEngine(4);
%     boundary.MinimumDistanceReferenceFrequencyType('User');
%     boundary.FrequencyForMinimumDistance('fqwextraspace');
    project.StoreParameter('openboundary_distance', CST.Defaults.OpenBoundaryDistance);
    boundary.MinimumDistanceType('Absolute');
    boundary.SetAbsoluteDistance('openboundary_distance');

    %% Set up scanning settings
    project.StoreParameter(CST.Defaults.ThetaName, 0);
    project.StoreParameter(CST.Defaults.PhiName, 0);
    boundary.PeriodicUseConstantAngles(1);
    boundary.SetPeriodicBoundaryAngles(CST.Defaults.ThetaName, CST.Defaults.PhiName);
%     boundary.YPeriodicShift('phi');

    material = project.Material();
    material.Reset();
    material.Name('Transparent');
    material.Colour(1, 1, 1);
    material.Transparency(75);
    material.Create();

    % Set Background material to normal instead of PEC.
    material = project.Material();
    material.Reset();
    material.Type('Normal');
    material.ChangeBackgroundMaterial();

    try % Not supported by all CST versions.
        %% Enable bounding box rendering.
        plot = project.Plot();
        plot.DrawBox(1);
        plot.DrawWorkplane(0);

        %% Rotate the view.
        plot.RestoreView('Bottom');
        plot.RotationAngle(30); plot.Rotate('Right');
        plot.RotationAngle(15); plot.Rotate('Down');
        plot.Update();

        %% Set shape accuracy
        solid = project.Solid();
        solid.ShapeVisualizationAccuracy2('120')

        plot.ZoomToStructure();
    catch ex
    end

    %% Enable Y and Z matrix calculation.
    postprocess1d = project.PostProcess1D();
    postprocess1d.ActivateOperation('YZ-matrices', 1);
end