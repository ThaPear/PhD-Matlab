classdef OpenLine < Element
    properties
        er
        L
    end
    methods
        function this = OpenLine(er, L)
            this.er = er;
            this.L = L;
        end
        function zin = GetInputImpedance(this, isTE, f, k0, kr)
            kd = k0 ./ sqrt(this.er);
            krdsq = krdsqfunc(this.er);
            
            kzd = -1j .* sqrt(-(kd.^2 - krdsq));
            
            [~, zcte, zctm] = z(this.er, kd, kzd);
            
            if(isTE)
                zc = zcte;
            else
                zc = zctm;
            end
            
            zin = -1j .* zc .* cot(kzd .* this.L);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
            error('%s::GetABCD:\n\tShould not be called on OpenLine.\n\tABCD matrix is not valid.', mfilename);
        end
        function h = GetHeight(this)
            h = this.L;
        end
    end
end