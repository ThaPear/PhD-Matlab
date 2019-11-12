classdef Impedance
    properties
        Z
    end
    methods
        function obj = Impedance(Z)
            if(nargin < 1)
                Z = nan;
            end
            obj.Z = Z;
        end
        function zin = GetInputImpedance(obj, isTE, f, k0, kr)
            zin = obj.Z;
        end
    end
end