classdef Capacitance < Impedance
    properties
        C
    end
    methods
        function obj = Capacitance(C)
            obj.C = C;
        end
        function zin = GetInputImpedance(obj, isTE, f, k0, kr)
            zin = 1 ./ (1j.*2.*pi.*f.*obj.C); % 1 / (jwC)
        end
    end
end