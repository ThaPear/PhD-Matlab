% ADS.BuildCST
function BuildCST(this, project)
    wcs = project.WCS();
    solid = project.Solid();
    component = project.Component();
    brick = project.Brick();
    material = project.Material();
    transform = project.Transform();

    number = component.GetNextFreeNameWithBase('adl');
    adli = ['adl', num2str(number)];
    componentname = adli;
    component.New(componentname);
    

    % If dx and dy don't exist, create them to be equal to p.
    project.MakeSureParameterExists('dx', this.p*1e3);
    project.MakeSureParameterExists('dy', this.p*1e3);
    project.MakeSureParameterExists('adl_p', this.p*1e3);
    project.MakeSureParameterExists('adl_s0', 'slot_s0 * (dx / adl_p) - int(slot_s0 * (dx / adl_p))');
    project.MakeSureParameterExists('adl_h0', 0);
    project.MakeSureParameterExists('slot_s0', Globals.slot_s0);

    % adli_N
    % adli_h
    % adli_w
    % adli_s0
    % adli_s
    
    % Set s0 as parameter.
%     adlis0 = [adli, '_s0'];
%     if(number == 1)
%         project.StoreParameter(adlis0, ['adl_s0']);
%     else
%         % adli_s0 = (adl(i-1)_s0 + adl(i-1)_s * adl(i-1)_N) Mod 1
%         project.StoreParameter(adlis0, ...
%             [   '(adl', num2str(number-1), '_s0 + adl', num2str(number-1), '_s * adl', num2str(number-1), '_N) - ', ...
%              'int(adl', num2str(number-1), '_s0 + adl', num2str(number-1), '_s * adl', num2str(number-1), '_N)']);
%     end
    % Set h0 as parameter.
%     adlih0 = ['adl', num2str(number), '_h0'];
%     if(number == 1)
%         project.StoreParameter(adlih0, ['adl_h0']);   
%     else
        % adli_h0 = adl(i-1)_h0 + adl(i-1)_h
%         project.StoreParameter(adlih0, ...
%             ['adl', num2str(number-1), '_h0 + adl', num2str(number-1), '_h']);
%     end

    % Set h0 directly in history list.
    adlis0 = 'adl_s0';
    for(n = 1:number-1)
        adlis0 = [adlis0, ' + adl', num2str(n), '_s * adl', num2str(n), '_N']; %#ok<AGROW>
    end
    adlis0 = ['((', adlis0, ') - int(', adlis0, '))'];
    
    % Set h0 directly in history list.
    adlih0 = 'adl_h0';
    for(n = 1:number-1)
        adlih0 = [adlih0, ' + adl', num2str(n), '_h']; %#ok<AGROW>
    end
    adlih0 = ['(', adlih0, ')'];

    adlip = 'adl_p';%[adli, '_p']; project.StoreParameter(adlip, 'adl_p');
    adliN = [adli, '_N']; project.StoreParameter(adliN, this.NADL);
    adlih = [adli, '_h']; project.StoreParameter(adlih, this.GetHeight() * 1e3);
    adliw = [adli, '_w']; project.StoreParameter(adliw, this.ws(1) * 1e3);
    adlis = [adli, '_s']; project.StoreParameter(adlis, this.ss(1) * 1e3);
    % Host medium is only adjustable if it's not an ADS_Real
    if(~contains(class(this), '_Real'))
        adlierhost = [adli, '_erhost']; project.StoreParameter(adlierhost, this.erhosts(1));
    end
    
    % Create unit cell for later.
    brick.Reset();
    brick.Component(componentname);
    brick.Name('UnitCell');
    brick.Xrange('-dx/2', 'dx/2');
    brick.Yrange('-dy/2', 'dy/2');
    brick.Zrange(adlih0, [adlih0, ' + ', adlih]);
    brick.Material('Vacuum');
    brick.Create();

    %% Create the metal patches.
    % Create a single patch.
    brick.Reset();
    brick.Component(componentname);
    brick.Name('Metal');
    brick.Xrange(['-dx/2 + ', adliw, '/2'], ['-dx/2 + ', adlip, ' - ', adliw, '/2']);
    brick.Yrange(['-dy/2 + ', adliw, '/2'], ['-dy/2 + ', adlip, ' - ', adliw, '/2']);
    brick.Zrange([adlih0, ' + ', adlih, ' / ', adliN, ' / 2'], ...
                 [adlih0, ' + ', adlih, ' / ', adliN, ' / 2']);
    brick.Material('PEC');
    brick.Create();
    
    % Transform it back far enough.
    transform.Reset();
    transform.Name([componentname, ':Metal']);
    transform.Material('PEC');
    % (s0 - 1 - floor(N * s)) * p
    transform.Vector(['(', adlis0, ' - 1 - int(', adliN, ' * ', adlis, ')) * ', adlip], ...
                     ['(', adlis0, ' - 1 - int(', adliN, ' * ', adlis, ')) * ', adlip], ...
                     0);
    transform.MultipleObjects(0);
    transform.GroupObjects(0);
    transform.Repetitions(1);
    transform.Transform('Shape', 'Translate');
    
    % Copy it in x-direction.
    transform.Reset();
    transform.Name([componentname, ':Metal']);
    transform.Material('PEC');
    transform.Vector([adlip], 0, 0);
    transform.MultipleObjects(1);
    transform.GroupObjects(1);
    transform.Repetitions(['dx / ', adlip, ' + int(', adlis, ' * ', adliN, ')']);
    transform.Transform('Shape', 'Translate');
    
    % Copy it in y-direction.
    transform.Reset();
    transform.Name([componentname, ':Metal']);
    transform.Material('PEC');
    transform.Vector(0, [adlip], 0);
    transform.MultipleObjects(1);
    transform.GroupObjects(1);
    transform.Repetitions(['dx / ', adlip, ' + int(', adlis, ' * ', adliN, '+ 0.99)']);
    transform.Transform('Shape', 'Translate');
    
    % Copy it in z-direction.
    transform.Reset();
    transform.Name([componentname, ':Metal']);
    transform.Material('PEC');
    transform.Vector([adlis, ' * ', adlip], [adlis, ' * ', adlip], [adlih, ' / ', adliN]);
    transform.MultipleObjects(1);
    transform.GroupObjects(1);
    transform.Repetitions([adliN, ' - 1']);
    transform.Transform('Shape', 'Translate');

    % Boolean the metal with the fullsized block.
    solid.Intersect([componentname, ':Metal'], [componentname, ':UnitCell']);
    
    %% Introduce the host medium.
    erhost = this.erhosts(1);
    % Create necessary material
    material.Reset();
    material.Name(['" & ', adlierhost, ' & "']);
    material.Folder('Generated');
    clr = ['" & 0.8-', adlierhost, '/10 & "'];
    material.Colour(clr, clr, clr);
    material.Epsilon(['" & ', adlierhost, ' & "']);
    material.Transparency(50);
    material.CreateConditional(['Not Material.Exists("Generated/" & ', adlierhost, ')']);
    project.AddToHistory('conditional material History Update', ['With Material', newline, '.Epsilon(', adlierhost, ')', newline, 'End With']);
    
    % Create the host medium brick.
    brick.Reset();
    brick.Component(componentname);
    brick.Name('Host');
    brick.Xrange('-dx/2', 'dx/2');
    brick.Yrange('-dy/2', 'dy/2');
    brick.Zrange(adlih0, [adlih0, ' + ', adlih]);
    brick.Material(['Generated/" & ', adlierhost, ' & "']);
    brick.Create();
end
