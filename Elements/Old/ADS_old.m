% Suppress warnings:
    % 'The variable '...' appears to change size on every loop iteration. 
    %  Consider preallocating for speed.'
    %#ok<*AGROW>

% Artificial Dielectric Slab
classdef ADS < TLine
    properties(SetAccess = protected)
        p       % Unit cell size.
        ds      % Distances between layer i and i+1.
        ss      % Shift between layer i and i+1.
        ws
        erhosts
    end
    methods
        function this = ADS(p, ds, ss, ws, erhosts)
            this.p = p;
            this.ds = ds;
            this.ss = ss;
            this.ws = ws;
            
            N = length(ws);
            
            if(nargin < 4)
                % If not specified, assume free space.
                warning('erhosts not specified, assuming free space');
                erhosts = ones(1, N+1);
            elseif(length(erhosts) == 1)
                % If only a single value was given, duplicate it.
                warning('Single-length erhost supplied, duplicating.');
                erhosts = ones(1, N+1) .* erhosts;
            end
            this.erhosts = erhosts;
            
            if(length(ss) < N-1)
                error('Insufficient ss specified.');
            end
            if(length(ds) < N+1)
                error('Insufficient ds specified.');
            end
            if(length(erhosts) < N+1)
                error('Insufficient erhosts specified.');
            end
            
            elements = {}; % Start with an empty transmission line.
            if(N == 1)
                elements = [elements, {Line(erhosts(1), ds(1))}];
                elements = [elements, {ADL(p,                   ... % p
                                    inf, inf,                   ... % d(n-1,n) , d(n,n+1)
                                    0, 0,                       ... % s(n-1,n), s(n,n+1)
                                    0, ws(1), 0,                ... % w(n-1), w(n), w(n+1)
                                    erhosts(1), erhosts(2))}];      % er(n-1,n), er(n,n+1)
                elements = [elements, {Line(erhosts(2), ds(2))}];
            else
%                 ADL(p, 
%                 dprev, dnext, 
%                 sprev, snext, 
%                 wprev, w, wnext,
%                 erhostdown, erhostup)
                %% First layer
                elements = [elements, {Line(erhosts(1), ds(1))}];
                elements = [elements, {ADL(p,                   ... % p
                                    inf, ds(2),                 ... % d(n-1,n) , d(n,n+1)
                                    0, ss(1),                   ... % s(n-1,n), s(n,n+1)
                                    0, ws(1), ws(2),            ... % w(n-1), w(n), w(n+1)
                                    erhosts(1), erhosts(2))}];      % er(n-1,n), er(n,n+1)
                elements = [elements, {Line(erhosts(2), ds(2))}]; 

                %% Middle layers
                for(n = 2:N-1)
                    elements = [elements, {ADL(p,               ... % p
                                        ds(n), ds(n+1),         ... % d(n-1,n), d(n,n+1)
                                        ss(n-1), ss(n),         ... % s(n-1,n), s(n,n+1)
                                        ws(n-1), ws(n), ws(n+1),... % w(n-1), w(n), w(n+1)
                                        erhosts(n), erhosts(n+1))}];% er(n-1,n), er(n,n+1)
                    elements = [elements, {Line(erhosts(n+1), ds(n+1))}]; 
                end
                
                %% Final layer
                elements = [elements, {ADL(p,                   ... % p
                                    ds(end-1), inf,             ... % d(n-1,n), d(n,n+1)
                                    ss(end), 0,                 ... % s(n-1,n), s(n,n+1)
                                    ws(end-1), ws(end), 0,      ... % w(n-1), w(n), w(n+1)
                                    erhosts(end-1), erhosts(end))}];
                % Half-width line on the end.
                elements = [elements, {Line(erhosts(end), ds(end))}];
            end
            
            this.elements = elements;
        end
        function [this, ads] = SetNextADS(this, ads)
            lastADL = this.elements{end-1};
            nextADL = ads.elements{2};
            
            dtotal = this.ds(end) + ads.ds(1);
            
            lastADL.dnext = dtotal;
            nextADL.dprev = dtotal;
            
            ermean = (this.ds(end)*this.erhosts(end) + ads.ds(1)*ads.erhosts(1)) / dtotal;
            
            lastADL.ernext = ermean;
            nextADL.erprev = ermean;
            
            this.elements{end-1} = lastADL;
            ads.elements{2} = nextADL;
        end
        function str = PrintParameters(this, f0)
            lambda0 = Constants.c0 / f0;
            [epsilon, ~] = this.GetEpsilonMu(f0, 0, 0);
            er = epsilon.x;
            lambda = lambda0 / sqrt(er);
            
            str = '';
            if(nargout == 0)
                delim = '\n';
            else
                delim = ' - ';
            end
            % TODO: Fix erhost when varying per-layer.
            str = [str, [char(949), 'rhost = ', num2str(this.erhosts(1))], ', '];
            str = [str, [char(949), 'r = ', num2str(er)], delim];
            str = [str, VarStr('p', this.p, lambda0, lambda, delim), delim];
            str = [str, VarStr('ds', this.ds, lambda0, lambda, delim), delim];
            str = [str, VarStr('ss', this.ss, lambda0, lambda, delim), delim];
            str = [str, VarStr('ws', this.ws, lambda0, lambda, delim), delim];
            str = [str, VarStr('h', this.GetHeight(), lambda0, lambda, delim), delim];
            if(nargout == 0)
                dispex(str);
            end
            function str = VarStr(name, values, lambda0, lambda, delim)
                if(all(values == values(1)))
                    value = values(1);
                    str = [name, ' = ', num2str(value*1e3, '%.5f'), 'mm = ', num2str(value/lambda0, '%.5f'), char(955), '0 = ', num2str(value/lambda, '%.5f'), char(955), 'eff'];
                else % Values are different for each layer.
                    mmstr = [name, ' = ['];
                    l0str = [repmat(' ', 1, length(name)*2+2), ' = ['];
                    lestr = [repmat(' ', 1, length(name)*2+2), ' = ['];
                    for(i = 1:length(values))
                        value = values(i);
                        mmstr = [mmstr, num2str(value*1e3, '%.5f'), ', '];
                        l0str = [l0str, num2str(value/lambda0, '%.5f'), ', '];
                        lestr = [lestr, num2str(value/lambda, '%.5f'), ', '];
                    end
                    mmstr = [mmstr, '] mm'];
                    l0str = [l0str, '] ', char(955), '0'];
                    lestr = [lestr, '] ', char(955), 'eff'];
                    str = [mmstr, delim, l0str, delim, lestr];
                end
            end
        end
        function PlotEpsilons(this, fs, hAx)
            if(nargin < 3 || isempty(hAx))
                hFig = figureex;
                hAx = hFig.CurrentAxes;
            end
            
            th = eps; ph = eps;
            ers = zeros(size(fs));
            parfor(fi = 1:length(fs))
                f = fs(fi);
                [epsilon, ~] = this.GetEpsilonMu(f, th, ph); %#ok<PFBNS>
                ers(fi) = epsilon.x;
            end
            plot(hAx, fs/1e9, ers);
        end
        function BuildCST(this, project)
            wcs = project.WCS();
            solid = project.Solid();
            component = project.Component();
            brick = project.Brick();
            material = project.Material();
            transform = project.Transform();
            
            global s0;
            if(isempty(s0))
                s0 = 0.5;
            end
            
            number = component.GetNextFreeNameWithBase('ADS');
            componentname = ['ADS', num2str(number)];
            component.New(componentname);
            
            % If dx and dy don't exist, create them to be equal to p.
            project.MakeSureParameterExists('dx', this.p*1e3);
            project.MakeSureParameterExists('dy', this.p*1e3);
            
            dx = str2double(project.RestoreParameter('dx'))/1e3;
            dy = str2double(project.RestoreParameter('dy'))/1e3;
            Nx = round(dx / this.p);
            Ny = round(dy / this.p);
            
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
                            materialname = num2str(line.er, 5);
                            material.Name(materialname);
                            material.Folder('Generated');
                            material.Colour(0, min(1, line.er/20), 1);
                            material.Epsilon(line.er);
                            material.Transparency(0.5);
                            material.Create();

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
                        transform.Vector((-this.p/2*Nx+this.p*s0)*1e3, (-this.p/2*Nx+this.p*s0)*1e3, 0);
                        transform.MultipleObjects(0);
                        transform.GroupObjects(0);
                        transform.Repetitions(1);
                        transform.Transform('Shape', 'Translate');
                        
                        % Copy the plates for following translates
                        transform.MultipleObjects(1);
                        transform.GroupObjects(1);
                        
                        % Copy the plate Nx times.
                        transform.Vector(this.p*1e3, 0, 0);
                        transform.Repetitions(Nx);
                        transform.Transform('Shape', 'Translate');
                        
                        % Copy the plate Ny times.
                        transform.Vector(0, this.p*1e3, 0);
                        transform.Repetitions(Ny);
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
        function h = GetEffectiveHeight(this, f)
            [epsilon, ~] = this.GetEpsilonMu(f, 0, 0);
            h = this.GetHeight() .* sqrt(epsilon.x);
        end
        function [ds, ss, ws] = CondenseParams(this, ds, ss, ws)
            if(nargin < 2); ds = [0]; end
            if(nargin < 3); ss = [];  end
            if(nargin < 4); ws = [];  end
            
            ds(end) = ds(end) + this.ds(1);
            ds = [ds, this.ds(2:end)];
            ss = [ss, this.ss];
            ws = [ws, this.ws];
            
            % Assume the last shift specified is used between the ADSs.
            % This is a correct assumption.
            if(length(this.ws) > 1)
                ss = [ss, this.ss(end)];
            end
        end
        function [ds, ss, ws] = CondenseParamsSeparate(this, ds, ss, ws)
            if(nargin < 2); ds = []; end
            if(nargin < 3); ss = [];  end
            if(nargin < 4); ws = [];  end
            
            ds = [ds, NaN, this.ds];
            ss = [ss, NaN, this.ss];
            ws = [ws, NaN, this.ws];
            
            % Assume the last shift specified is used between the ADSs.
            % This is a correct assumption.
            if(length(this.ws) > 1)
                ss = [ss, this.ss(end)];
            end
        end
        %%
%        function er = GetSimpleEpsilon(obj, isTE, f, th, ph)
%{
            [k0, ~, ~, kz0] = k(f, 1, th, ph);
            % elements{1} is a Line object.
            [kd, ~, ~, kz] = k(f, obj.elements{1}.er, th, ph);
            
            kr = kd .* sin(th);
            
            dz = obj.elements{1}.L * 2;
            % elements{2} is an ADL object.
            Bn = obj.elements{2}.GetLayerSusceptance(f);
            
            if(isTE)
                Zadl = 1 ./ (1j .* (Bn .* (1-sin(th).^2./2)));
            else
                Zadl = 1 ./ (1j .* Bn);
            end
            
            % C = 1/Zr sinh(kz dz)
%             line = TLine({obj.elements{1}, obj.elements{2}, obj.elements{1}});
%             ABCD = line.GetABCD(isTE, f, th, ph, kz0);
%             Zr = 1 ./ (ABCD.C ./ sinh(kz .* dz));
            Zr = Constants.z0 ./ sqrt(obj.elements{1}.er);
            
            kzeff = 1 ./ dz .* acos(cos(kz .* dz) + 1j .* Zr ./ Zadl .* sin(kz .* dz));
            
            er = (kzeff.^2 + kr.^2) ./  k0.^2;
            
        end
        function [fs, desc] = GetResonances(obj, fmin, fmax)
            fs = {};
            desc = {};
            [epsilon, mu] = obj.GetEpsilonMu((fmin + fmax) / 2, 0, 0);
            er = epsilon.x;
            for(i = 1:length(obj.elements))
                element = obj.elements{i};
                if(isa(element, 'ADL'))
                    [fsel, descel] = element.GetResonances(fmin, fmax, er);
                else
                    [fsel, descel] = element.GetResonances(fmin, fmax);
                end
                descel = strcat({[class(element), ':']}, descel);
                fs = [fs, fsel];
                desc = [desc, descel];
            end
        end
%}
    end
end