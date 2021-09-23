% ADS.BuildCST
function BuildCST(this, project, parentcomponent)
    if(nargin < 3 || isempty(parentcomponent))
        parentcomponent = '';
    else
        if(~strcmp(parentcomponent(end), '/'))
            parentcomponent = [parentcomponent, '/'];
        end
    end
    wcs = project.WCS();
    solid = project.Solid();
    component = project.Component();
    brick = project.Brick();
    material = project.Material();
    transform = project.Transform();

    number = component.GetNextFreeNameWithBase([parentcomponent, 'ADS']);
    componentname = [parentcomponent, 'ADS', num2str(number)];
    component.New(componentname);

    % If dx and dy don't exist, create them to be equal to p.
    project.MakeSureParameterExists('dx', this.p*1e3);
    project.MakeSureParameterExists('dy', this.p*1e3);
%     project.MakeSureParameterExists('p', this.p*1e3);
    project.MakeSureParameterExists('adl_s0', 0);

%     global s0;
%     if(isempty(s0))
%         s0 = 'adl_s0';
%     end
%     if(isempty(s0))
%         if(Globals.exists('slot_s0'))
%             s0 = 0;%Globals.slot_s0 * (dx / this.p) + 0.5;
%             while(s0 > 0.5)
%                 s0 = s0 - 1;
%             end
%         else
%             s0 = 0.5;
%         end
%     end
    project.MakeSureParameterExists('slot_s0', 0);

    wcs.Enable();

    h = this.GetHeight();
    h_total = '0';
    
    function [paramname, i] = GetParameterName(base)
        i = 0;
        while(project.DoesParameterExist(sprintf(base, i)))
            i = i + 1;
        end
        paramname = sprintf(base, i);
    end

    [param_adl_p, iADS] = GetParameterName('adl%i_p');
    project.StoreParameter(param_adl_p, this.p*1e3);
    param_adl = param_adl_p(1:end-1);
    
    %% Determine shift of previous layers.
    s0 = 'dx*slot_s0 + adl0_p*adl_s0';
    for(i = 0:iADS-1)
        s0 = [s0, ' + ', sprintf('adl%i_p * adl%i_stotal', i, i)]; %#ok<AGROW>
    end
    
    % Number of repetitions for the patches in x and y.
    param_Nx = ['dx / ', param_adl_p];
    param_Ny = ['dy / ', param_adl_p];

    stotal = '0';
    htotal = '';
    mademetal = 0;
    for(i = 1:length(this.elements))
        brick.Reset();
        brick.Component(componentname);
        element = this.elements{i};
        switch(class(element))
            case 'Line'
                line = element;
                
                param_adl_h = GetParameterName([param_adl, 'h%i']);
                project.StoreParameter(param_adl_h, line.L*1e3);
                
                h_total = [h_total, ' + ', param_adl_h]; %#ok<AGROW>
                h = param_adl_h;
                if(h == 0)
                    continue;
                end
%                 if(line.er ~= 1)
                    % Create necessary material
                    if(length(unique(line.er)) > 1 || length(unique(line.mu)) > 1)
                        materialname = [num2str(line.er(1), 5), '_', num2str(line.er(2), 5), '_', num2str(line.er(3), 5)];
                        material.Name(materialname);
                        material.Folder('Generated');
                        material.Colour(0, min(1, line.er(1)/20), 1);
                        material.Type('Anisotropic');
                        material.EpsilonX(line.er(1));
                        material.EpsilonY(line.er(2));
                        material.EpsilonZ(line.er(3));
                        material.MuX(line.mu(1));
                        material.MuY(line.mu(2));
                        material.MuZ(line.mu(3));
                        material.Transparency(0.5);
                        material.Create();
                    else
                        materialname = [num2str(line.er(1), 5)];
                        material.Name(materialname);
                        material.Folder('Generated');
                        material.Colour(0, min(1, line.er(1)/20), 1);
                        material.Type('Normal');
                        material.Epsilon(line.er(1));
                        material.Mu(line.mu(1));
                        material.Transparency(0.5);
                        material.Create();
                    end

                    % Create dielectric slab.
                    linename = ['Line ', num2str(i)];
                    brick.Name(linename);
                    brick.Xrange('-dx/2', 'dx/2');
                    brick.Yrange('-dy/2', 'dy/2');
                    brick.Zrange(0, h);
                    brick.Material(['Generated/', materialname]);
                    brick.Create();

                    solid.MergeMaterialsOfComponent([componentname, ':', linename]);
%                 end

                wcs.MoveWCS('local', 0, 0, h);

            case {'ADL', 'ADL_Real'}
                adl = element;
                if(mademetal)
                    name = 'Metal2';
                else
                    name = 'Metal';
                end

                param_adl_w = GetParameterName([param_adl, 'w%i']);
                project.StoreParameter(param_adl_w, adl.w*1e3);
                
                brick.Name(name);
                brick.Xrange(sprintf('-%s/2+%s/2', param_adl_p, param_adl_w), sprintf('%s/2-%s/2', param_adl_p, param_adl_w));
                brick.Yrange(sprintf('-%s/2+%s/2', param_adl_p, param_adl_w), sprintf('%s/2-%s/2', param_adl_p, param_adl_w));
                brick.Zrange(0, 0);
                brick.Material('PEC');
                brick.Create();

                % Move the plate to the corner of the ADS.
                transform.Reset();
                transform.Name([componentname, ':', name]);
                transform.Material('PEC');
%                                 % -p/2*(dx/p) + dx*slot_s0 + s0
%                 transform.Vector(['-', param_adl_p, '/2*(dx/', param_adl_p, ') + dx*slot_s0 + ', s0], ...
%                                  ['-', param_adl_p, '/2*(dy/', param_adl_p, ') + dy*slot_s0 + ', s0], 0);
                                % -p/2*(dx/p) + dx*slot_s0 + p*stotal + s0
                                
                % Calc: s = p*stotal + s0
                s = sprintf('(%s)*(%s) + %s', param_adl_p, stotal, s0);
                % Calc: s = (s - int(s/p)*p)
                s = sprintf('(%s - int((%s)/%s)*%s)', s, s, param_adl_p, param_adl_p);
                transform.Vector(sprintf('-%s/2*(dx/%s) + %s', ...
                                     param_adl_p, param_adl_p, s), ...
                                 sprintf('-%s/2*(dy/%s) + %s', ...
                                     param_adl_p, param_adl_p, s), ...
                                 0);
                transform.MultipleObjects(0);
                transform.GroupObjects(0);
                transform.Repetitions(1);
                transform.Transform('Shape', 'Translate');

                % Copy the plates for following translates
                transform.MultipleObjects(1);
                transform.GroupObjects(1);

                % Copy the plate Nx+1 times.
                transform.Vector(param_adl_p, 0, 0);
                transform.Repetitions([param_Nx, ' + 1']);
                transform.Transform('Shape', 'Translate');

                % Copy the plate Ny+1 times.
                transform.Vector(0, param_adl_p, 0);
                transform.Repetitions([param_Ny, ' + 1']);
                transform.Transform('Shape', 'Translate');

                param_adl_s = GetParameterName([param_adl, 's%i']);
                project.StoreParameter(param_adl_s, adl.snext/this.p);
                
                if(length(stotal) < 2)
                    stotal = param_adl_s;
                    htotal = param_adl_h;
                else
                    stotal = [stotal, ' + ', param_adl_s];
                    htotal = [htotal, ' + ', param_adl_h];
                end

                if(mademetal)
                    solid.Add([componentname, ':Metal'], [componentname, ':Metal2']);
                end

                mademetal = 1;
            otherwise
                error('Invalid element found in ADS');
        end
    end
    % Store the total shift in a parameter.
    param_stotal = [param_adl, 'stotal'];
    project.StoreParameter(param_stotal, stotal);
    
    % Store the total height in a parameter.
    % TODO: Make dynamic depending on adl%i_N
    param_htotal = [param_adl, 'htotal'];
    project.StoreParameter(param_htotal, htotal);

    brick.Reset();
    brick.Component(componentname);
    brick.Name('UnitCell');
    brick.Xrange('-dx/2', 'dx/2');
    brick.Yrange('-dy/2', 'dy/2');
    brick.Zrange(['-(', h_total, ')'], 0);
    brick.Material('Transparent');
    brick.Create();
    % Boolean the metal with the fullsized block.
    solid.Intersect([componentname, ':Metal'], [componentname, ':UnitCell']);
end
