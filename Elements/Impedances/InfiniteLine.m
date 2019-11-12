classdef InfiniteLine < Element
    properties
        er
    end
    methods
        function obj = InfiniteLine(er)
            obj.er = er;
        end
        function zin = GetInputImpedance(obj, isTE, f, k0, kr)
            kd = k0 .* sqrt(obj.er);
            krdsq = krdsqfunc(obj.er);
            
            kzd = -1j .* sqrt(-(kd.^2 - krdsq));
            
            [~, zcte, zctm] = z(obj.er, kd, kzd);
            
            if(isTE)
                zin = zcte;
            else
                zin = zctm;
            end
        end
        function ABCD = GetABCD(obj, isTE, f, k0, kr)
            error('%s::GetABCD:\n\tShould not be called on InfiniteLine.\n\tABCD matrix is not valid.', mfilename);
        end
    end
end