function [eps, mu] = GetEpsilonMu(tline, f)
    project = CST.InitializeUnitCellProject();
    tline.BuildCST(project);
    
    if(f >= 1e9)
        f = f / 1e9;
    end
    fdsolver = project.FDSolver();
    fdsolver.ResetSampleIntervals('all');
    fdsolver.AddSampleInterval(f, '', 1, 'Single', 1);
    
    project.StoreParameter('fmesh', f);
    project.StoreParameter('fmin', f-1);
    project.StoreParameter('fmax', f+1);
    
    floquetport = project.FloquetPort();
    floquetport.StartBulkMode();
    floquetport.Port('Zmin');
    floquetport.SetDistanceToReferencePlane(['-openboundary_distance+', num2str(tline.elements{end}.L*1e3, '%.15g')]);
    floquetport.Port('Zmax');
    floquetport.SetDistanceToReferencePlane(['-openboundary_distance+', num2str(tline.elements{end}.L*1e3, '%.15g')]);
    floquetport.EndBulkMode();
%     project.StoreParameter('openboundary_distance', [CST.Defaults.OpenBoundaryDistance, '+', num2str(tline.elements{end}.L*1e3, '%.15g')]);

    project.Rebuild();
    success = project.Rebuild();
    if(~success)
        fprintf('%s: Project rebuild failed.\n', mfilename);
        return;
    end
    
    parametersweep = project.ParameterSweep();
    parametersweep.StartBulkMode();
    parametersweep.SetSimulationType('Frequency');
    parametersweep.AddSequence('Sequence 1');
    parametersweep.AddParameter_ArbitraryPoints('Sequence 1', 'aa_theta', '0;10');
    parametersweep.AddParameter_ArbitraryPoints('Sequence 1', 'aa_phi', '0');
    parametersweep.EndBulkMode();
    parametersweep.Start();
    
%     success = fdsolver.Start();
%     if(~success)
%         fprintf('%s: Simulation failed.\n', mfilename);
%         return;
%     end
    
%     filename = 'Temp/testS.s4p';
%     CST.ExportResult(project, '1D Results/S-Parameters', filename);
%     [fs, parameters, S] = CST.LoadData(filename);
end