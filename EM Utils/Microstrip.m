classdef Microstrip
    methods(Static)
        function WoverD = GetWoverD(z0des, er)
            A = z0des / 60 * sqrt((er+1)/2) + (er-1)/(er+1) * (0.23+0.11/er);
            WoverD = 8*exp(A) / (exp(2*A)-2);
            
            if(WoverD > 2)
                B = 377*pi/(2*z0des*sqrt(er));
                WoverD = 2/pi*(B-1-log(2*B-1)+(er-1)/(2*er)*(log(B-1)+0.39-0.61/er));
            end
        end
    end
end