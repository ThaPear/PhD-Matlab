% Suppress warnings:
    % 'The variable '...' appears to change size on every loop iteration. 
    %  Consider preallocating for speed.'
    %#ok<*AGROW>

% Artificial Dielectric Slab with real layers
classdef ADS_Real < ADS
    properties(SetAccess = protected)
%         p       % Unit cell size.
%         ds      % Distances between layer i and i+1.
%         ss      % Shift between layer i and i+1.
%         ws
    end
    properties(Dependent)
%         NADL
    end
    methods
        BuildCST(this, project)
    end
    methods
        function this = ADS_Real(p, ds, ss, ws, ~)
            this.p = p;
            this.ds = ds;
            this.ss = ss;
            this.ws = ws;
            
            NADL = this.NADL;
            
            if(length(ss) < NADL)
                error('Insufficient ss specified.');
            end
            if(length(ds) < NADL+1)
                error('Insufficient ds specified.');
            end
            
            elements = {}; % Start with an empty transmission line.
            if(NADL == 1)
                % Foam-Glue-Metal-Substrate-Glue-Foam
                % Foam
                elements = [elements, {Line(Materials.Foam.permittivity, ds(1) - Materials.Glue.thickness)}];
                % Glue
                elements = [elements, {Line(Materials.Glue.permittivity, Materials.Glue.thickness)}];
                % Metal
                elements = [elements, {ADL_Real(p,                   ... % p
                                    inf, inf,                   ... % d(n-1,n) , d(n,n+1)
                                    0, 0,                       ... % s(n-1,n), s(n,n+1)
                                    0, ws(1), 0,                ... % w(n-1), w(n), w(n+1)
                                    Materials.Foam.permittivity, Materials.Foam.permittivity)}];      % er(n-1,n), er(n,n+1)
                % Substrate
                elements = [elements, {Line(Materials.ADLSubstrate.permittivity, Materials.ADLSubstrate.thickness)}];
                % Glue
                elements = [elements, {Line(Materials.Glue.permittivity, Materials.Glue.thickness)}];
                % Foam
                elements = [elements, {Line(Materials.Foam.permittivity, ds(2) - Materials.Glue.thickness  - Materials.ADLSubstrate.thickness)}];
            else
%                 ADL_Real(p, 
%                 dprev, dnext, 
%                 sprev, snext, 
%                 wprev, w, wnext,
%                 erhostdown, erhostup)
                % First layer
                % Foam
                elements = [elements, {Line(Materials.Foam.permittivity, ds(1) - Materials.Glue.thickness)}];
                % Glue
                elements = [elements, {Line(Materials.Glue.permittivity, Materials.Glue.thickness)}];
                % Metal
                elements = [elements, {ADL_Real(p,                   ... % p
                                    inf, ds(2),                 ... % d(n-1,n) , d(n,n+1)
                                    0, ss(1),                   ... % s(n-1,n), s(n,n+1)
                                    0, ws(1), ws(2),            ... % w(n-1), w(n), w(n+1)
                                    Materials.Foam.permittivity, Materials.Foam.permittivity)}];      % er(n-1,n), er(n,n+1)
                % Substrate
                elements = [elements, {Line(Materials.ADLSubstrate.permittivity, Materials.ADLSubstrate.thickness)}];
                % Glue
                elements = [elements, {Line(Materials.Glue.permittivity, Materials.Glue.thickness)}];
                % Foam
                elements = [elements, {Line(Materials.Foam.permittivity, ds(2) - Materials.Glue.thickness - Materials.ADLSubstrate.thickness)}]; 

                % Middle layers
                for(n = 2:NADL-1)
                    % Glue
                    elements = [elements, {Line(Materials.Glue.permittivity, Materials.Glue.thickness)}];
                    % Metal
                    elements = [elements, {ADL_Real(p,               ... % p
                                        ds(n), ds(n+1),         ... % d(n-1,n), d(n,n+1)
                                        ss(n-1), ss(n),         ... % s(n-1,n), s(n,n+1)
                                        ws(n-1), ws(n), ws(n+1),... % w(n-1), w(n), w(n+1)
                                        Materials.Foam.permittivity, Materials.Foam.permittivity)}];% er(n-1,n), er(n,n+1)
                    % Substrate
                    elements = [elements, {Line(Materials.ADLSubstrate.permittivity, Materials.ADLSubstrate.thickness)}];
                    % Glue
                    elements = [elements, {Line(Materials.Glue.permittivity, Materials.Glue.thickness)}];
                    % Foam
                    elements = [elements, {Line(Materials.Foam.permittivity, ds(n+1) - 2*Materials.Glue.thickness - Materials.ADLSubstrate.thickness)}]; 
                end
                
                % Glue
                elements = [elements, {Line(Materials.Glue.permittivity, Materials.Glue.thickness)}];
                % Final layer
                elements = [elements, {ADL_Real(p,              ... % p
                                    ds(end-1), inf,             ... % d(n-1,n), d(n,n+1)
                                    ss(end-1), ss(end),         ... % s(n-1,n), s(n,n+1)
                                    ws(end-1), ws(end), 0,      ... % w(n-1), w(n), w(n+1)
                                    Materials.Foam.permittivity, Materials.Foam.permittivity)}];
                % Substrate
                elements = [elements, {Line(Materials.ADLSubstrate.permittivity, Materials.ADLSubstrate.thickness)}];
                % Glue
                elements = [elements, {Line(Materials.Glue.permittivity, Materials.Glue.thickness)}];
                % Foam - Half-width line on the end.
                elements = [elements, {Line(Materials.Foam.permittivity, ds(end) - 2*Materials.Glue.thickness - Materials.ADLSubstrate.thickness)}];
            end
            
            this.elements = elements;
        end
    end
end