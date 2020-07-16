classdef FreeSpace < Element
    properties
    end
    methods
        function zin = GetInputImpedance(obj, isTE, f, k0, kr)
            z0 = Constants.z0;
            
            kz0 = -1j .* sqrt(-(k0.^2 - kr.^2));
            
            if(isTE)
                zin = z0 .* k0 ./ kz0;
            else
                zin = z0 .* kz0 ./ k0;
            end
        end
        function ABCD = GetABCD(obj, isTE, f, k0, kr)
            error('%s::GetABCD:\n\tShould not be called on FreeSpace.\n\tABCD matrix is not valid.', mfilename);
        end
        function BuildCST(obj, project, parentcomponent)
            % Nothing to do.
        end
    end
end