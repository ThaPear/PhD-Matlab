classdef Nothing < Element
    properties
    end
    methods
        function obj = Nothing()
        end
        function Z = GetInputImpedance(this, isTE, f, k0, kr)
            error(['%s::GetInputImpedance:\n',...
                   '\tInput impedance is not valid on a Nothing element.\n', ...
                   '\tUse a shunt element or an open/shorted/terminated line.'], mfilename);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
            ABCD = ABCDMatrix(1, 0, ...
                              0, 1);
        end
        function h = GetHeight(this)
            h = 0;
        end
    end
end