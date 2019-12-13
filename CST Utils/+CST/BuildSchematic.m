% Builds the schematic view for N ports, adds capacitances and external ports at the specified
% values.
% If capacitance = inf or not provided, capacitors are not placed.
function BuildSchematic(dsproject, externalimpedance, capacitance)
    if(nargin < 3)
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
        
        dsproject = CST.Application.ActiveDS();
        externalimpedance = 80;
    end
        

    gridsize = 25;
    capdistance = 18*gridsize;
    portdistance = 8*gridsize;
    
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
    
    
    pinnames = cell(1,Nports);
    port = [];
    for(i = 1:Nports)
        pinnames{i} = block.GetPinName(i-1);
        port(i).relpos = struct('x', block.GetPinPositionX(i-1) - baseposX, ...
                                'y', block.GetPinPositionY(i-1) - baseposY);
        if(abs(port(i).relpos.y) > abs(port(i).relpos.x))
            if(port(i).relpos.y < 0)
                % Port is above block
                port(i).orientation = 'up';
            else
                % Port is below block
                port(i).orientation = 'down';
            end
        elseif(port(i).relpos.x < 0)
            % Port is left of block
            port(i).orientation = 'left';
        else
            % Port is right of block
            port(i).orientation = 'right';
        end
    end
    
    if(isinf(capacitance))
        % No capacitances, so just link the ports to the block directly.
        for(i = 1:Nports)
            port(i).block = 'MWSSCHEM1';
%             port(i).pin = num2str(i);
            port(i).pin = pinnames{i};
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
            
%             block.Position(baseposX + (2*i-3)*5 * gridsize, baseposY + 5 * gridsize);
            
            switch(port(i).orientation)
                case 'left'
                    block.Rotate(180);
                    block.Position(baseposX - capdistance, baseposY + port(i).relpos.y);
                case 'right'
                    block.Position(baseposX + capdistance, baseposY + port(i).relpos.y);
                case 'up'
                    block.Rotate(-90);
                    block.Position(baseposX + port(i).relpos.x, baseposY - capdistance);
                case 'down'
                    block.Rotate(90);
                    block.Position(baseposX + port(i).relpos.x, baseposY + capdistance);
                otherwise
                    error('Invalid port orientation');
            end
            
            block.Create();
            
            %% Link the capacitances to the array.
            link.Reset();
            link.SetConnection(1, 'B', name, '1');
            link.SetConnection(2, 'B', 'MWSSCHEM1', pinnames{i});
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
        %externalport.Position(baseposX + (2*i-3)*5 * gridsize, baseposY + 10 * gridsize);
        switch(port(i).orientation)
            case 'left'
                externalport.Position(baseposX - capdistance - portdistance, baseposY + port(i).relpos.y);
            case 'right'
                externalport.Position(baseposX + capdistance + portdistance, baseposY + port(i).relpos.y);
            case 'up'
                externalport.Position(baseposX + port(i).relpos.x, baseposY - capdistance - portdistance);
            case 'down'
                externalport.Position(baseposX + port(i).relpos.x, baseposY + capdistance + portdistance);
            otherwise
                error('Invalid port orientation');
        end
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






















