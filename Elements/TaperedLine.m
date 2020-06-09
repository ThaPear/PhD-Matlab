classdef TaperedLine < TLine
    properties(SetAccess = protected)
        z1      % Impedance 1
        z2      % Impedance 2
        L       % Total length
        N       % Number of layers
    end
    methods
        % L is total length, N is number of steps.
        function this = TaperedLine(z1, z2, L, N)
            this.z1 = z1;
            this.z2 = z2;
            this.N = N;
            
            % Determine the required ers.
            a = 1 / L * log(z2 / z1);
            Z = @(z) z1 * exp(a * z);
            len = L/N;
            ers = zeros(1,N);
            for(n = 1:N)
                ers(n) = (Constants.z0/Z(n*len - len/2))^2;
            end
            % Example ers for N=5.
            % ers = [37.9533 16.9163 7.5398 3.3606 1.4979];

            elements = cell(1,N);
            for(n = 1:N)
                er = ers(n);
                elements{n} = Line(L/N ./ sqrt(er), er);
            end
            this.elements = elements;
        end
        function OutputCST(this, f0)
            lambda0 = Constants.c0 / f0;
            for(i = 1:length(this.elements))
                element = this.elements{i};
                
                dispex('MakeHomogeneousLayer("Layer X", %.10f, p, h0, %.10f*l0, Nx, Ny)\n', ...
                    element.er, element.L/lambda0);
            end
        end
    end
end