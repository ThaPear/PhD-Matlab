classdef Shunt < Element
    properties
        impedance
    end
    methods
        function this = Shunt(impedance)
            this.impedance = impedance;
        end
        function zin = GetInputImpedance(this, isTE, f, k0, kr)
            zin = this.impedance.GetInputImpedance(isTE, f, k0, kr);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
            Y = 1 ./ this.impedance.GetInputImpedance(isTE, f, k0, kr);
            ABCD = ABCDMatrix(1, 0, ...
                              Y, 1);
        end
    end
end