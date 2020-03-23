classdef Element
    properties
    end
    methods
        function Zin = GetInputImpedance(this, isTE, f, k0, kr)
            error('%s::GetInputImpedance:\n\tShould not be called on an Element object.\n\tUse a derived class.', mfilename);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
            error('%s::GetABCD:\n\tShould not be called on an Element object.\n\tUse a derived class.', mfilename);
        end
        function h = GetHeight(this)
            h = 0;
        end
        function h = GetEffectiveHeight(this, f)
            h = 0;
        end
        function BuildCST(this, project)
            error('BuildCST not defined for objects of class ''%s''.', class(this));
        end
    end
end