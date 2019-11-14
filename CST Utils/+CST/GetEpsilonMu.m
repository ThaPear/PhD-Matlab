function [eps, mu] = GetEpsilonMu(tline, f)
    project = CST.InitializeUnitCellProject();
    tline.BuildCST(project);
    
    if(f >= 1e9)
        f = f / 1e9;
    end
    fdsolver.ResetSampleIntervals('all');
    fdsolver.AddSampleInterval(f, '', 1, 'Single', 1);
    
    project.StoreParameter('fmin', f-0.1);
    project.StoreParameter('fmax', f+0.1);
    
    success = fdsolver.Start();
    if(~success)
        fprintf('%s: Simulation failed.\n', mfilename);
        return;
    end
    
    
end