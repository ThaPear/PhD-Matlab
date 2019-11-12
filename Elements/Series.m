classdef Series < Element
    properties
        impedance
    end
    methods
        function this = Series(impedance)
            this.impedance = impedance;
        end
        function zin = GetInputImpedance(this, isTE, f, k0, kr)
            error('%s::GetInputImpedance:\n\Input impedance is not valid on a series element.\n\tUse a shunt element or an open/shorted line.', mfilename);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
            Z = this.impedance.GetInputImpedance(isTE, f, k0, kr);
            ABCD = ABCDMatrix(1, Z, ...
                              0, 1);
        end
    end
end