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
    project.MakeSureParameterExists('p', this.p*1e3);
    project.MakeSureParameterExists('adl_s0', -1);

    dx = str2double(project.RestoreParameter('dx'))/1e3;
    dy = str2double(project.RestoreParameter('dy'))/1e3;
    Nx = round(dx / this.p);
    Ny = round(dy / this.p);

    global s0;
    if(isempty(s0))
        if(Globals.exists('slot_s0'))
            s0 = 0;%Globals.slot_s0 * (dx / this.p) + 0.5;
            while(s0 > 0.5)
                s0 = s0 - 1;
            end
        else
            s0 = 0.5;
        end
    end
    project.MakeSureParameterExists('slot_s0', s0/2);

    wcs.Enable();

    h = this.GetHeight();

    brick.Reset();
    brick.Component(componentname);
    brick.Name('UnitCell');
    brick.Xrange('-dx/2', 'dx/2');
    brick.Yrange('-dy/2', 'dy/2');
    brick.Zrange(0, h*1e3);
    brick.Material('Transparent');
    brick.Create();

    mademetal = 0;
    for(i = 1:length(this.elements))
        brick.Reset();
        brick.Component(componentname);
        element = this.elements{i};
        switch(class(element))
            case 'Line'
                line = element;
                h = line.L*1e3;
                if(h == 0)
                    continue;
                end
                if(line.er ~= 1)
                    % Create necessary material
                    if(length(line.er) > 1)
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
                        material.Colour(0, min(1, line.er/20), 1);
                        material.Type('Anisotropic');
                        material.Epsilon(line.er);
                        material.Mu(line.mu);
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
                end

                wcs.MoveWCS('local', 0, 0, h);

            case {'ADL', 'ADL_Real'}
                adl = element;
                if(mademetal)
                    name = 'Metal2';
                else
                    name = 'Metal';
                end

                brick.Name(name);
                brick.Xrange((-this.p/2+adl.w/2)*1e3, (this.p/2-adl.w/2)*1e3);
                brick.Yrange((-this.p/2+adl.w/2)*1e3, (this.p/2-adl.w/2)*1e3);
                brick.Zrange(0, 0);
                brick.Material('PEC');
                brick.Create();

                % Move the plate to the corner of the ADS.
                transform.Reset();
                transform.Name([componentname, ':', name]);
                transform.Material('PEC');
                transform.Vector(['-p/2*(dx/p)+dx*slot_s0 + p*(adl_s0 + ', num2str(s0), ' -', num2str(s0), '\1)'], ...
                                 ['-p/2*(dy/p)+dx*slot_s0 + p*(adl_s0 + ', num2str(s0), ' -', num2str(s0), '\1)'], 0);
%                         transform.Vector((-this.p/2*Nx+this.p*s0)*1e3, (-this.p/2*Nx+this.p*s0)*1e3, 0);
                transform.MultipleObjects(0);
                transform.GroupObjects(0);
                transform.Repetitions(1);
                transform.Transform('Shape', 'Translate');

                % Copy the plates for following translates
                transform.MultipleObjects(1);
                transform.GroupObjects(1);

                % Copy the plate Nx+1 times.
                transform.Vector(this.p*1e3, 0, 0);
                transform.Repetitions(Nx+1);
                transform.Transform('Shape', 'Translate');

                % Copy the plate Ny+1 times.
                transform.Vector(0, this.p*1e3, 0);
                transform.Repetitions(Ny+1);
                transform.Transform('Shape', 'Translate');

                s0 = s0 + adl.snext/this.p;
                if(s0 > 0.5)
                    s0 = s0 - 1;
                end

                if(mademetal)
                    solid.Add([componentname, ':Metal'], [componentname, ':Metal2']);
                end

                mademetal = 1;
            otherwise
                error('Invalid element found in ADS');
        end
    end

    % Boolean the metal with the fullsized block.
    solid.Intersect([componentname, ':Metal'], [componentname, ':UnitCell']);
end
