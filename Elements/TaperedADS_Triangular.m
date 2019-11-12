% Artificial Dielectric Slab
classdef TaperedADS_Triangular < StagedADS
    properties(SetAccess = protected)
        z1      % Impedance 1
        z2      % Impedance 2
    end
    methods
        % L is total length, N is number of layers.
        function this = TaperedADS_Triangular(p, z1, z2, L, N, f0, useADLs)
            % Determine the required ers.
            a = @(z) log(z2 / z1) * ((z <  L/2) .* (          2*(z/L).^2    ) + ...
                                     (z >= L/2) .* (4 * z/L - 2*(z/L).^2 - 1) ...
                                    );
            Z = @(z) z1 * exp(a(z));
            len = L/N;
            ers = zeros(1,N);
            for(n = 1:N)
                ers(n) = (Constants.z0/Z(n*len - len/2))^2;
            end
            heights = len ./ sqrt(ers);
            
            if(nargin < 7)
                useADLs = 1;
            end
            
            this@StagedADS(p, heights, ers, useADLs*ones(1, N), f0);
            
            this.z1 = z1;
            this.z2 = z2;
        end
    end
end