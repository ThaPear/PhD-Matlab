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
    properties(Dependent)
        NADL
    end
    methods
        function this = ADS(p, ds, ss, ws, erhosts)
            if(nargin == 0)
                return;
            end
            this.p = p;
            this.ds = ds;
            this.ss = ss;
            this.ws = ws;
            
            NADL = this.NADL;
            
            if(nargin < 4)
                % If not specified, assume free space.
                warning('erhosts not specified, assuming free space');
                erhosts = ones(1, NADL+1);
            elseif(length(erhosts) == 1)
                % If only a single value was given, duplicate it.
                warning('Single-length erhost supplied, duplicating.');
                erhosts = ones(1, NADL+1) .* erhosts;
            end
            this.erhosts = erhosts;
            
            if(length(ss) < NADL)
                error('Insufficient ss specified.');
            end
            if(length(ds) < NADL+1)
                error('Insufficient ds specified.');
            end
            if(length(erhosts) < NADL+1)
                error('Insufficient erhosts specified.');
            end
            
            elements = {}; % Start with an empty transmission line.
            if(NADL == 1)
                elements = [elements, {Line(ds(1), erhosts(1))}];
                elements = [elements, {ADL(p,                   ... % p
                                    inf, inf,                   ... % d(n-1,n) , d(n,n+1)
                                    0, 0,                       ... % s(n-1,n), s(n,n+1)
                                    0, ws(1), 0,                ... % w(n-1), w(n), w(n+1)
                                    erhosts(1), erhosts(2))}];      % er(n-1,n), er(n,n+1)
                elements = [elements, {Line(ds(2), erhosts(2))}];
            else
%                 ADL(p, 
%                 dprev, dnext, 
%                 sprev, snext, 
%                 wprev, w, wnext,
%                 erhostdown, erhostup)
                %% First layer
                elements = [elements, {Line(ds(1), erhosts(1))}];
                elements = [elements, {ADL(p,                   ... % p
                                    inf, ds(2),                 ... % d(n-1,n) , d(n,n+1)
                                    0, ss(1),                   ... % s(n-1,n), s(n,n+1)
                                    0, ws(1), ws(2),            ... % w(n-1), w(n), w(n+1)
                                    erhosts(1), erhosts(2))}];      % er(n-1,n), er(n,n+1)
                elements = [elements, {Line(ds(2), erhosts(2))}]; 

                %% Middle layers
                for(n = 2:NADL-1)
                    elements = [elements, {ADL(p,               ... % p
                                        ds(n), ds(n+1),         ... % d(n-1,n), d(n,n+1)
                                        ss(n-1), ss(n),         ... % s(n-1,n), s(n,n+1)
                                        ws(n-1), ws(n), ws(n+1),... % w(n-1), w(n), w(n+1)
                                        erhosts(n), erhosts(n+1))}];% er(n-1,n), er(n,n+1)
                    elements = [elements, {Line(ds(n+1), erhosts(n+1))}]; 
                end
                
                %% Final layer
                elements = [elements, {ADL(p,                   ... % p
                                    ds(end-1), inf,             ... % d(n-1,n), d(n,n+1)
                                    ss(end-1), 0,                 ... % s(n-1,n), s(n,n+1)
                                    ws(end-1), ws(end), 0,      ... % w(n-1), w(n), w(n+1)
                                    erhosts(end-1), erhosts(end))}];
                % Half-width line on the end.
                elements = [elements, {Line(ds(end), erhosts(end))}];
            end
            
            this.elements = elements;
        end
        function n = get.NADL(this)
            n = length(this.ws);
        end
        function [varargout] = SetNeighbours(this, prevADS, nextADS)
            % [...] = SetNeighbours(this, prevads, nextads)
            % Sets the parameters on the edge elements on the interface
            % between this and prevads, and between this and nextads
            %
            % [...] = SetNeighbours(this, prevads)
            % Sets the parameters on the edge elements on the interface
            % between this and prevads
            %
            % [this, prevads, nextads] = SetNeighbours(...)
            % Sets the parameters on the edge elements of the previous and
            % next ADS as well.
            %
            % [this, prevads] = SetNeighbours(...)
            % Sets the parameters on the edge elements of the previous ADS
            % as well.
            %
            % [this] = SetNeighbours(...)
            % Only sets the parameters on the current ads.
            %
            dprev = 0;
            dnext = 0;
            iprev = 0;
            ifirst = 0;
            ilast = 0;
            inext = 0;
            
            if(~isempty(prevADS))
                % Find last ADL in previous ADS
                for(i = prevADS.N:-1:1)
                    if(contains(class(prevADS.elements{i}), 'ADL')) % Also catches ADL_Real
                        iprev = i;
                        break;
                    else
                        dprev = dprev + prevADS.elements{i}.GetHeight();
                    end
                end
                % Find first ADL in this ADS
                for(i = 1:this.N)
                    if(contains(class(this.elements{i}), 'ADL'))
                        ifirst = i;
                        break;
                    else
                        dprev = dprev + this.elements{i}.GetHeight();
                    end
                end
                % Last ADL of prev ADS
                prevADS.elements{iprev}.dnext = dprev;
                prevADS.elements{iprev}.wnext = this.elements{ifirst}.w;
                prevADS.elements{iprev}.snext = prevADS.ss(end);
                % First ADL of this ADS
                this.elements{ifirst}.dprev = dprev;
                this.elements{ifirst}.wprev = prevADS.elements{iprev}.w;
                this.elements{ifirst}.sprev = prevADS.ss(end);
            end
            
            if(nargin == 3 && ~isempty(nextADS))
                % Find last ADL in this ADS
                for(i = this.N:-1:1)
                    if(contains(class(this.elements{i}), 'ADL'))
                        ilast = i;
                        break;
                    else
                        dnext = dnext + this.elements{i}.GetHeight();
                    end
                end
                % Find first ADL in next ADS
                for(i = 1:nextADS.N)
                    if(contains(class(nextADS.elements{i}), 'ADL'))
                        inext = i;
                        break;
                    else
                        dnext = dnext + nextADS.elements{i}.GetHeight();
                    end
                end
                % Last ADL of this ADS
                this.elements{ilast}.dnext = dnext;
                this.elements{ilast}.wnext = nextADS.elements{inext}.w;
                this.elements{ilast}.snext = this.ss(end);
                % First ADL of next ADS
                nextADS.elements{inext}.dprev = dnext;
                nextADS.elements{inext}.wprev = this.elements{ilast}.w;
                nextADS.elements{inext}.sprev = this.ss(end);
            end
            
            
            if(nargout == 1)
                varargout = {this};
            elseif(nargout == 2)
                varargout = {this, prevADS};
            else
                varargout = {this, prevADS, nextADS};
            end
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
            if(~isempty(this.erhosts))
                str = [str, [char(949), 'rhost = ', num2str(this.erhosts(1))], ', '];
            end
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
            
            global s0;
            if(isempty(s0))
                if(Globals.exists('slot_s0'))
                    s0 = Globals.slot_s0 * (dx / this.p) + 0.5;
                    while(s0 > 0.5)
                        s0 = s0 - 1;
                    end
                else
                    s0 = 0.5;
                end
            end
            
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
    end
end