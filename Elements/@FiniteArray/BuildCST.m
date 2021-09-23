function BuildCST(this, project, parentcomponent)
    if(nargin < 3 || isempty(parentcomponent))
        parentcomponent = '';
    else
        if(~strcmp(parentcomponent(end), '/'))
            parentcomponent = [parentcomponent, '/'];
        end
    end
    componentname = [parentcomponent, 'Array'];

    this.unitcell.BuildCST(project, [componentname, '/UnitCell']);

    wcs = project.WCS();
    wcs.Enable(); 
    wcs.Store('Pre-Slot'); 
    wcs.RotateWCS('u', 180); 

    % Build down-stratification. 
    this.tlinedown.BuildCST(project, [componentname, '/Stratification']); 

    wcs.Restore('Pre-Slot'); 

    % Build up-stratification. 
    this.tlineup.BuildCST(project, [componentname, '/Stratification']); 

    wcs.Restore('Pre-Slot'); 
    wcs.Delete('Pre-Slot'); 
    wcs.Disable(); 

    project.StoreParameter('lambda_min', 'c0/fmin/1e9');
    project.StoreParameter('Nx', num2str(this.Nx, '%.15g'));
    project.StoreParameter('Ny', num2str(this.Ny, '%.15g'));
    project.StoreParameter('edge_distance', num2str(this.dedge*1e3, '%.15g'));
    project.StoreParameter('edge_length', '5/3 * sqr(slot_width * lambda_min)');
%     project.StoreParameter('padding_x', 'max(0, lambda_min/4-edge_length)');
%     project.StoreParameter('padding_y', 'max(0, lambda_min/4-(dy/2-slot_width/2))');
    % Round up the X padding to a multiple of dx, to ensure a good fit of the stratification
    % above the unit cells.
    %     goal_x = max(-Int(-(lambda_min/4+edge_distance)/dx)*dx, -Int(-(edge_length+edge_distance-dx/2)/dx)*dx)
    %     padding_x = goal_x - edge_length - edge_distance
    project.StoreParameter('padding_x', 'max(-Int(-(lambda_min/4)/dx)*dx, -Int(-(edge_length+edge_distance-dx/2)/dx)*dx) - edge_length - edge_distance + dx/2');
    % Round up the Y padding to a multiple of dy, to ensure a good fit of the stratification
    % above the unit cells.
    %     goal_y = max(-Int(-(lambda_min/4)/dx)*dx, -Int(-(dy/2-slot_width/2)/dx)*dx)
    %     padding_y = goal_y
    project.StoreParameter('padding_y', 'max(-Int(-(lambda_min/4)/dx)*dx, -Int(-(dy/2-slot_width/2)/dx)*dx)');

    transform = project.Transform();
    brick = project.Brick();
    solid = project.Solid();

    %% Copy the unit cell in the x-direction along the array
    % Only if there's more than 1 feed in X.
    project.NextCommandConditional('Nx > 1');
        transform.Reset();
        transform.Name([componentname, '/UnitCell']);
        transform.Vector('dx', 0, 0);
        transform.MultipleObjects(1);
        transform.GroupObjects(1);
        transform.Repetitions('Nx-1');
        transform.Transform('Shape', 'Translate');
    project.NextCommandConditional('Nx > 1');
        transform.Reset();
        transform.Name([componentname, '/Stratification']);
        transform.Vector('dx', 0, 0);
        transform.MultipleObjects(1);
        transform.GroupObjects(1);
        transform.Repetitions('Nx-1');
        transform.Transform('Shape', 'Translate');

    %% Copy the stratification above and below along the X-terminations and X-padding.
    transform.Reset();
    transform.Name([componentname, '/Stratification']);
    transform.Vector('-dx', 0, 0);
    transform.MultipleObjects(1);
    transform.GroupObjects(1);
    transform.Repetitions('-Int(0.01-(edge_distance+edge_length+padding_x-dx/2)/dx)'); % -Int(-x) = Ceiling(x)
    transform.Transform('Shape', 'Translate');

    transform.Reset();
    transform.Name([componentname, '/Stratification']);
    transform.Vector('dx', 0, 0);
    transform.MultipleObjects(1);
    transform.GroupObjects(1);
    transform.Repetitions('-Int(0.01-(edge_distance+edge_length+padding_x-dx/2)/dx)'); % -Int(-x) = Ceiling(x)
    transform.Transform('Shape', 'Translate');

    %% Create the terminations
    brick.Reset();
    brick.Name('Termination1');
    brick.Component([componentname, '/UnitCell/Slot']);
    brick.Xrange('-edge_distance-edge_length',  '-dx/2');
    brick.Yrange('-dy/2', 'dy/2');
    brick.Zrange('0', '0');
    brick.Material('PEC');
    brick.Create();

    brick.Reset();
    brick.Name('Termination1_Slot');
    brick.Component([componentname, '/UnitCell/Slot']);
    brick.Xrange('-edge_distance',  '-dx/2');
    brick.Yrange('-slot_width/2', 'slot_width/2');
    brick.Zrange('0', '0');
    brick.Material('PEC');
    brick.Create();

    brick.Reset();
    brick.Name('Termination2');
    brick.Component([componentname, '/UnitCell/Slot']);
    brick.Xrange('(Nx-1)*dx+dx/2',  '(Nx-1)*dx+edge_distance+edge_length');
    brick.Yrange('-dy/2', 'dy/2');
    brick.Zrange('0', '0');
    brick.Material('PEC');
    brick.Create();

    brick.Reset();
    brick.Name('Termination2_Slot');
    brick.Component([componentname, '/UnitCell/Slot']);
    brick.Xrange('(Nx-1)*dx+dx/2',  '(Nx-1)*dx+edge_distance');
    brick.Yrange('-slot_width/2', 'slot_width/2');
    brick.Zrange('0', '0');
    brick.Material('PEC');
    brick.Create();

    solid.Subtract([componentname, '/UnitCell/Slot:Termination1'], [componentname, '/UnitCell/Slot:Termination1_Slot']);
    solid.Subtract([componentname, '/UnitCell/Slot:Termination2'], [componentname, '/UnitCell/Slot:Termination2_Slot']);

    %% Copy the slot in the y-direction
    % Only if there's more than 1 slot in Y.
    project.NextCommandConditional('Ny > 1');
        transform.Reset();
        transform.Name([componentname, '/UnitCell']);
        transform.Vector(0, 'dy', 0);
        transform.MultipleObjects(1);
        transform.GroupObjects(1);
        transform.Repetitions('Ny-1');
        transform.Transform('Shape', 'Translate');
    project.NextCommandConditional('Ny > 1');
        transform.Reset();
        transform.Name([componentname, '/Stratification']);
        transform.Vector(0, 'dy', 0);
        transform.MultipleObjects(1);
        transform.GroupObjects(1);
        transform.Repetitions('Ny-1');
        transform.Transform('Shape', 'Translate');

    %% Copy the stratification above and below along the Y-terminations and Y-padding.
    transform.Reset();
    transform.Name([componentname, '/Stratification']);
    transform.Vector(0, '-dy', 0);
    transform.MultipleObjects(1);
    transform.GroupObjects(1);
    transform.Repetitions('-Int(-padding_y/dy)'); % -Int(-x) = Ceiling(x)
    transform.Transform('Shape', 'Translate');

    transform.Reset();
    transform.Name([componentname, '/Stratification']);
    transform.Vector(0, 'dy', 0);
    transform.MultipleObjects(1);
    transform.GroupObjects(1);
    transform.Repetitions('-Int(-padding_y/dy)'); % -Int(-x) = Ceiling(x)
    transform.Transform('Shape', 'Translate');

    %% Copy the port in the x-direction
    % Only if there's more than 1 feed in X.
    project.NextCommandConditional('Nx > 1');
        transform.Reset();
        transform.Name('port1 (SlotFeed)');
        transform.Vector('dx', 0, 0);
        transform.MultipleObjects(1);
        transform.Repetitions('Nx-1');
        transform.Transform('Port', 'Translate');

    %% Copy the port in the x-direction
    % Only if there's more than 1 slot in Y.
    project.NextCommandConditional('Ny > 1');
        % Copy the port of each nx.
        project.NextCommandLoop('nxi', '1', 'Nx')
            transform.Reset();
            transform.Name('port" & nxi & " (SlotFeed)');
            transform.Vector(0, 'dy', 0);
            transform.MultipleObjects(1);
            transform.Repetitions('Ny-1');
            transform.Transform('Port', 'Translate');

    %% Add the padding to ensure lambda/4 spacing from the slot to the open boundary
    project.NextCommandConditional('padding_x > 0');
        brick.Reset();
        brick.Name('Padding_x1');
        brick.Component(componentname);
        brick.Xrange('-edge_distance-edge_length-padding_x',  '-edge_distance-edge_length');
        brick.Yrange('-dy/2', '(Ny-1)*dy+dy/2');
        brick.Zrange('0', '0');
        brick.Material('PEC');
        brick.Create();
    project.NextCommandConditional('padding_x > 0');
        brick.Reset();
        brick.Name('Padding_x2');
        brick.Component(componentname);
        brick.Xrange('(Nx-1)*dx+edge_distance+edge_length',  '(Nx-1)*dx+edge_distance+edge_length+padding_x');
        brick.Yrange('-dy/2', '(Ny-1)*dy+dy/2');
        brick.Zrange('0', '0');
        brick.Material('PEC');
        brick.Create();
    project.NextCommandConditional('padding_y > 0');
        brick.Reset();
        brick.Name('Padding_y1');
        brick.Component(componentname);
        brick.Xrange('-edge_distance-edge_length-padding_x', '(Nx-1)*dx+edge_distance+edge_length+padding_x');
        brick.Yrange('-dy/2-padding_y',  '-dy/2');
        brick.Zrange('0', '0');
        brick.Material('PEC');
        brick.Create();
    project.NextCommandConditional('padding_y > 0');
        brick.Reset();
        brick.Name('Padding_y2');
        brick.Component(componentname);
        brick.Xrange('-edge_distance-edge_length-padding_x', '(Nx-1)*dx+edge_distance+edge_length+padding_x');
        brick.Yrange('(Ny-1)*dy+dy/2',  '(Ny-1)*dy+dy/2+padding_y');
        brick.Zrange('0', '0');
        brick.Material('PEC');
        brick.Create();

    % Create the vertical walls if applicable.
    if(this.unitcell.walled)
        brick.Reset()
        brick.Name('Wall')
        brick.Component(componentname);
        brick.Xrange('-edge_distance-edge_length-padding_x', '(Nx-1)*dx+edge_distance+edge_length+padding_x');
        brick.Yrange('-dy/2','-dy/2');
        brick.Zrange('-hback', '0');
        brick.Material('PEC');
        brick.Create();

        transform.Reset();
        transform.Name([componentname, ':Wall']);
        transform.Vector(0, 'dy', 0);
        transform.MultipleObjects(1);
        transform.GroupObjects(1);
        transform.Repetitions('Ny');
        transform.Transform('Shape', 'Translate');
    end
end