% Artificial Dielectric Slab
classdef StagedADS < TLine
    properties(SetAccess = protected)
        p       % Unit cell size
        heights % Heights of the sections
        ers     % Permittivities of stages
        isADLs  % Is this stage an ADL or a Line
    end
    properties(Dependent)
        NADS       % Number of layers
    end
    methods
        % L is total length, N is number of stages.
        function this = StagedADS(p, heights, ers, isADLs, f0)
            if(nargin == 0)
                return;
            end
            tc = tic;
            this.p = p;
            this.heights = heights;
            this.ers = ers;
            this.isADLs = isADLs;
            
            NADS = length(ers);
            
            adss = cell(1, NADS);
            for(rep = 1:2)
%                 disp(['Repetition ', num2str(rep)]);
                for(n = 1:NADS)
                    er = ers(n);
                    h = heights(n);
                    if(isADLs(n))
                        nextlayer = [];
                        prevlayer = [];
                        if(n > 1)
                            if(isa(adss{n-1}, 'ADS'))
                                prevlayer = adss{n-1};
                            end
                        end
                        if(length(adss) > n) % Second repetition.
                            if(isa(adss{n+1}, 'ADS'))
                                nextlayer = adss{n+1};
                            end
                        end
                        adss{n} = DesignADS(f0, p, h, er, prevlayer, nextlayer);
                    else
                        adss{n} = Line(ers(n), h);
                    end
                    
%                     disp(['Designed layer ', num2str(n), ' in ', num2str(dt), 's, used ', num2str((length(elements{n}.elements)-1)/2), ' layers.']);
                end
            end
            
            this.elements = adss;
            this = this.ApplyNeighbours();
            
            % Count the number of ADLs in the entire ADS
            nlayers = 0;
            for(i = 1:NADS)
                if(~this.isADLs(i))
                    continue;
                end
                el = this.elements{i};
                nlayers = nlayers + el.NADL;
            end
            
            dt = toc(tc);
            fprintf('Designed %i-stage %i-layer %.2fmm %s in %.3fs.\n', ...
                NADS, nlayers, this.GetHeight()*1e3, class(this), dt);
        end
        function this = ApplyNeighbours(this)
            Nads = this.NADS;
            if(Nads > 1)
                % Correct the distances seen from each layer of the ADLs.
                % If the number of layers did not change in the 2nd iteration
                % of the design above, this step does not change anything, as
                % the distances will have been correctly set in DesignADS.
                for(n = 2:Nads-1)
                    if(this.isADLs(n-1) && this.isADLs(n) && this.isADLs(n+1))
                        [this.elements{n}, this.elements{n-1}, this.elements{n+1}] = ...
                            this.elements{n}.SetNeighbours(this.elements{n-1}, this.elements{n+1});
                    end
                end
                if(this.isADLs(Nads-1) && this.isADLs(Nads))
                    % Correct the last layer as well.
                    [this.elements{Nads}, this.elements{Nads-1}] = ...
                        this.elements{Nads}.SetNeighbours(this.elements{Nads-1});
                end
            end
        end
        function obj = Copy(this)
            %TODO: Change this if it's a handle.
            %obj = copy(this);
            obj = this;
        end
        function n = get.NADS(this)
            n = length(this.elements);
        end
        function PlotEpsilons(this, fs, hAx)
            if(nargin < 3)
                PlotEpsilons@TLine(this, fs);
            else
                PlotEpsilons@TLine(this, fs, hAx);
            end
            hAx = gca;
        
            % Plot the desired ers.
            for(n = 1:this.NADS)
                er = this.ers(n);
%                 plot(hAx, [min(fs) max(fs)]/1e9, [er er], 'k--');
%                 text(min(fs), er, ['Layer ', num2str(n)])
            end
        end
        function [ds, ss, ws] = CondenseParams(this, ds, ss, ws)
            if(nargin < 2); ds = [0]; end
            if(nargin < 3); ss = [];  end
            if(nargin < 4); ws = [];  end
            
            for(i = 1:this.NADS)
                [ds, ss, ws] = this.elements{i}.CondenseParams(ds, ss, ws);
            end
            % Chop off the last shift, since there's no layer after it.
            ss = ss(1:end-1);
        end
        function [ds, ss, ws] = CondenseParamsSeparate(this, ds, ss, ws)
            if(nargin < 2); ds = []; end
            if(nargin < 3); ss = [];  end
            if(nargin < 4); ws = [];  end
            
            for(i = 1:this.NADS)
                [ds, ss, ws] = this.elements{i}.CondenseParamsSeparate(ds, ss, ws);
            end
            % Chop off the last shift, since there's no layer after it.
            ss = ss(1:end-1);
        end
    end
end