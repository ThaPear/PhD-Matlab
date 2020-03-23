% Artificial Dielectric Slab
classdef ChebyshevADS < StagedADS
    properties(SetAccess = protected)
        z1      % Impedance 1
        z2      % Impedance 2
    end
    methods
        % L is total length, N is number of stages.
        function this = ChebyshevADS(p, gamma, z1, z2, N, f0, f0design, useADLs)
            lambda0 = 3e8/f0;
            len = 0.25*lambda0;
            
            [Zs] = Chebyshev.GetImpedances(gamma, z1, z2, N);
            
            ers = (z2 ./ Zs).^2;
            
            heights = len ./ sqrt(ers);
            
            if(nargin < 8)
                useADLs = 1;
            end
            
            this@StagedADS(p, heights, ers, useADLs*ones(1, N), f0design);
            
            this.z1 = z1;
            this.z2 = z2;
        end
    end
end