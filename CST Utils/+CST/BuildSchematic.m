% Builds the schematic view for N ports, adds capacitances and external ports at the specified
% values.
% If capacitance = inf or not provided, capacitors are not placed.
function BuildSchematic(project, dsproject, externalimpedance, capacitance)
    if(nargin < 4)
        capacitance = inf;
    end
    if(nargin == 0)
        ok = questdlg('No arguments provided, build in active project?','','Yes, with Capacitance','Yes, without Capacitance','No','Yes, with Capacitance');
        if(strcmp(ok, 'No'))
            dispex('Cancelled building schematic.\n');
            return;
        end
        if(contains(ok, 'without'))
            capacitance = inf;
        else
            capacitance = 1e-12;
        end
        
        project = CST.Application.Active3D();
        dsproject = CST.Application.ActiveDS();
        externalimpedance = 80;
    end
        

    gridsize = 25;
    
    block = dsproject.Block();
    link = dsproject.Link();
    units = dsproject.Units();
    
    capUnit = units.GetCapacitanceUnitToSI();

    %% Ensure the MWSSCHEM1 block exists and get its position.
    block.Reset();
    block.Name('MWSSCHEM1');
    if(~block.DoesExist())
        dispex('Block ''MWSSCHEM1'' not found.\n');
        return;
    end
    baseposX = block.GetPositionX();
    baseposY = block.GetPositionY();
    Nports = block.GetNumberOfPins();
    
    port = []
    
    if(isinf(capacitance))
        % No capacitances, so just link the ports to the block directly.
        for(i = 1:Nports)
            port(i).block = 'MWSSCHEM1';
            port2pin = num2str(i);
        end
    else
        for(i = 1:Nports)
            %% Create the capacitances.
            name = ['CAP', num2str(i)];
            block.Reset();
            block.Name(name);
            if(block.DoesExist())
                block.Delete();
                block.Name(name);
            end
            block.Type('CircuitBasic\Capacitor');
            block.SetDoubleProperty('Capacitance', capacitance/capUnit);
            block.Position(baseposX + (2*i-3)*5 * gridsize, baseposY + 5 * gridsize);
            block.Rotate(90);
            block.Create();
            
            %% Link the capacitances to the array.
            link.Reset();
            link.SetConnection(1, 'B', name, '1');
            link.SetConnection(2, 'B', 'MWSSCHEM1', i);
            link.Create();
            
            % Set the external ports to connect to the capacitances.
            port(i).block = name;
            port(i).pin = '2';
        end
    end
    
    for(i = 1:Nports)
        %% Create the external ports.
        externalport = dsproject.ExternalPort();
        externalport.Reset();
        externalport.Number(i);
        if(externalport.DoesExist())
            externalport.Delete();
            externalport.Number(i);
        end
        externalport.Position(baseposX + (2*i-3)*5 * gridsize, baseposY + 10 * gridsize);
        externalport.Create();
        externalport.SetImpedanceType(1);
        externalport.SetImpedance(externalimpedance);
        
        %% Link the ports to the appropriate blocks.
        link.Reset();
            link.SetConnection(1, 'B', port(i).block, port(i).pin);
            link.SetConnection(2, 'P', num2str(i), "");
        link.Create();
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end