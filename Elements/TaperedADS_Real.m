% Artificial Dielectric Slab
classdef TaperedADS_Real < StagedADS_Real
    properties(SetAccess = protected)
        z1      % Impedance 1
        z2      % Impedance 2
    end
    methods
        % L is total length, N is number of stages.
        function this = TaperedADS_Real(p, z1, z2, L, N, f0)
            % Determine the required ers.
            a = 1 / L * log(z2 / z1);
            Z = @(z) z1 * exp(a * z);
            len = L/N;
            ers = zeros(1,N);
            for(n = 1:N)
                ers(n) = (Constants.z0/Z(n*len - len/2))^2;
            end
            
            heights = len ./ sqrt(ers);
            
            this@StagedADS_Real(p, heights, ers, ones(1, N), f0);
            
            this.z1 = z1;
            this.z2 = z2;
        end
    end
end