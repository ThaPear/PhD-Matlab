classdef Inductance < Element
    properties
        L
    end
    methods
        function obj = Inductance(L)
            obj.L = L;
        end
        function zin = GetInputImpedance(obj, isTE, f, k0, kr)
            zin = (1j.*2.*pi.*f.*obj.L); % (jwL)
        end
%         function ABCD = GetABCD(obj, isTE, f, k0, kr)
%             Y = 1 ./ (1j.*2.*pi.*f.*obj.L); % 1 / (jwL)
%             ABCD = ABCDMatrix(1, 0, ...
%                               Y, 1);
%         end
    end
end