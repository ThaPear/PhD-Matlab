classdef Impedance
    properties
        Z
    end
    methods
        function this = Impedance(Z)
            if(nargin < 1)
                Z = nan;
            end
            this.Z = Z;
        end
        function zin = GetInputImpedance(this, isTE, f, k0, kr)
            zin = this.Z;
        end
    end
end