classdef Transformer < Element
    properties
        N
    end
    methods
        function this = Transformer(N)
            this.N = N;
        end
        function Z = GetInputImpedance(this, isTE, f, k0, kr)
            error('%s::GetInputImpedance:\n\Input impedance is not valid on a transformer.\n\tUse a shunt element or an open/shorted line.', mfilename);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
            n = this.N; 
            ABCD = ABCDMatrix(n, 0, ...
                              0, 1./n);
        end
    end
end