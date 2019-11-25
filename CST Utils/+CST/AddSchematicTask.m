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
    simulationtask.Type('s-parameters');
    simulationtask.Name('SPara1');
    simulationtask.Create();
end