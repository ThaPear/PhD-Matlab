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
            
            [Zs] = ChebyshevADS.GetImpedances(gamma, z1, z2, N);
            
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
    methods(Static)
        function [Zs] = GetImpedances(gamma, z1, z2, N)
            secthm = cosh(1/N*acosh(abs(log(z2/z1)/(2*gamma))));
%             secthm = cosh(1/N*acosh(1/gamma*abs((z2-z1)/(z2+z1))))
            
            polys = zeros(1, max(5, N+1));
                          % 1 x1 x2 x3 x4
            % T(1,x) = x
            polys(1, :) = [ 0  1  0  0  0 zeros(1, N-4)];
            % T(2,x) = 2x^2 - 1
            polys(2, :) = [-1  0  2  0  0 zeros(1, N-4)];
            % T(3,x) = 4x^3 - 3x
            polys(3, :) = [ 0 -3  0  4  0 zeros(1, N-4)];
            % T(4,x) = 8x^4 - 8x^2 + 1
            polys(4, :) = [ 1  0 -8  0  8 zeros(1, N-4)];
            
            for(n = 5:N)
                % T(n,x) = 2xT(n-1,x) - T(n-2,x)
                polys(n, :) = 2 .* [0, polys(n-1, 1:end-1)] - polys(n-2, :);
            end
            
%             A = (z2-z1)/(z2+z1)*1/sum(polys(N,:) .* secthm.^(0:max(4,N)));
            
            chebypol = zeros(1, N+1);
            for(i = 1:N+1)
                cosinecoeff = ChebyshevADS.ExpandCosinePower(i-1);
                cosinecoeff = [cosinecoeff, zeros(1, N+1-i)];
                chebypol = chebypol + polys(N, i) .* secthm^(i-1) .* cosinecoeff;
            end
            chebypol = chebypol * gamma;
            gammas = 1/2 * fliplr(chebypol(chebypol ~= 0));
            
            if(mod(N, 2))
                gammas = [gammas, fliplr(gammas)];
            else
                gammas = [gammas, fliplr(gammas(1:end-1))];
                gammas(ceil((N+1)/2)) = 2 * gammas(ceil((N+1)/2));
            end
            
            Zs(1) = exp(log(z1)+2*gammas(1));
            for(i = 2:N)
                Zs(i) = exp(2*gammas(i)) * Zs(i-1);
            end
        end
        function [mat] = ExpandCosinePower(n)
            % https://math.stackexchange.com/questions/117061/expansion-of-cosk-theta
            if(mod(n, 2)) % n is odd
                mat = zeros(1,n+1);
                i = n+1;
                for(k = 0:(n-1)/2)
                    % nchoosek doesn't accept matrices
                    mat(i) = 2^-(n-1) * nchoosek(n, k);
                    i = i - 2;
                end
            else % n is even
                mat = zeros(1,n+1);
                i = n+1;
                for(k = 0:(n/2-1))
                    % nchoosek doesn't accept matrices, so loop
                    mat(i) = 2^-(n-1) * nchoosek(n, k);
                    i = i - 2;
                end
                mat(1) = 2^-n * nchoosek(n, n/2);
            end
        end
    end
end