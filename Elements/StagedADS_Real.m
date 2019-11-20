% Artificial Dielectric Slab
classdef StagedADS_Real < StagedADS
    properties(SetAccess = protected)
%         p       % Unit cell size
%         heights % Heights of the sections
%         ers     % Permittivities of stages
%         isADLs  % Is this stage an ADL or a Line
    end
    properties(Dependent)
%         NADS       % Number of layers
    end
    methods
        % L is total length, N is number of stages.
        function this = StagedADS_Real(p, heights, ers, isADLs, f0)
            if(nargin == 0)
                return;
            end
            
            global hcut;
            heights(1) = heights(1) - hcut;
            
            tc = tic;
            this.p = p;
            this.heights = heights;
            this.ers = ers;
            this.isADLs = isADLs;
            
            NADS = length(ers);
            
            adss = cell(1, NADS);
            for(rep = 1:2)
%                 dispex('Repetition %s.\n' rep);
                for(n = 1:NADS)
                    er = ers(n);
                    h = heights(n);
                    if(isADLs(n))
                        nextlayer = [];
                        prevlayer = [];
                        if(n > 1)
                            if(isa(adss{n-1}, 'ADS_Real'))
                                prevlayer = adss{n-1};
                            end
                        end
                        if(length(adss) > n) % Second repetition.
                            if(isa(adss{n+1}, 'ADS_Real'))
                                nextlayer = adss{n+1};
                            end
                        end
                        doublespaceinfirstlayer = isempty(prevlayer);
                        adss{n} = DesignADS_Real(f0, p, h, er, prevlayer, nextlayer, doublespaceinfirstlayer);
                    else
                        adss{n} = Line(ers(n), h);
                    end
                    
%                     dispex('Designed layer %i in %fs, used %i layers.\n', n, dt, (length(adss{n}.elements)-1)/2);
                end
            end
            
            this.elements = adss;
            this = this.ApplyNeighbours();
            
            % Count the number of ADLs in the entire ADS
            nlayers = 0;
            for(i = 1:NADS)
                el = this.elements{i};
                nlayers = nlayers + el.NADL;
            end
            
            dt = toc(tc);
            dispex('Designed %i-stage %i-layer %.2fmm %s in %.3fs.\n', ...
                NADS, nlayers, this.GetHeight()*1e3, class(this), dt);
        end
    end
end