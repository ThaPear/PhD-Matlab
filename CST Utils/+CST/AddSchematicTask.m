% Creates a schematic task in the given DS project.
function AddSchematicTask(dsproject)
    if(nargin == 0)
        ok = questdlg('No arguments provided, add to active project?','','Yes','No','Yes');
        if(strcmp(ok, 'No'))
            dispex('Cancelled building schematic.\n');
            return;
        end
        dsproject = CST.Application.ActiveDS();
    end
    
    simulationtask = dsproject.SimulationTask();
    simulationtask.Reset();
    simulationtask.Name('SPara1');
    if(~simulationtask.DoesExist())
        simulationtask.Type('s-parameters');
        simulationtask.Create();
    end
    
    try
        parametersweep = dsproject.ParameterSweep();
        parametersweep.AddSequence('Sequence 1');
        parametersweep.AddParameter_ArbitraryPoints('Sequence 1', 'aa_phi', '0');
        parametersweep.AddParameter_ArbitraryPoints('Sequence 1', 'aa_theta', '0;60');
    catch exception
        % Ignore it if the sweep was already defined.
       if(~strcmpi(exception.identifier, 'MATLAB:COM:E2147549183'))
           rethrow(exception);
       end
    end
    
    
end