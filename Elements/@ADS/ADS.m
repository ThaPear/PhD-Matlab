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
        BuildCST(this, project, parentcomponent)
        BuildCSTMaterial(this, project)
        [ds, ss, ws] = CondenseParams(this, ds, ss, ws)
        [ds, ss, ws] = CondenseParamsSeparate(this, ds, ss, ws)
        h = GetEffectiveHeight(this, f)
        PlotEpsilons(this, fs, hAx)
        str = PrintParameters(this, f0)
        [varargout] = SetNeighbours(this, prevADS, nextADS)
    end
    methods
        function n = get.NADL(this)
            n = length(this.ws);
        end
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
                                    ss(end-1), ss(end),         ... % s(n-1,n), s(n,n+1)
                                    ws(end-1), ws(end), 0,      ... % w(n-1), w(n), w(n+1)
                                    erhosts(end-1), erhosts(end))}];
                % Half-width line on the end.
                elements = [elements, {Line(ds(end), erhosts(end))}];
            end
            
            this.elements = elements;
        end
    end
end