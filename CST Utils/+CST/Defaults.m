classdef Defaults
    properties(Constant)
        
        FrequencyRange = [9, 33];
        MeshFrequency = 29;
        SamplesPerGHz = 0.25; % (fmax-fmin)*N should be integer.
        ThetaName = 'aa_theta';
        PhiName   = 'aa_phi';
        OpenBoundaryDistance = 'clight/(fmin*1e9)/4*1e3';
    end
end