classdef Microstrip
    methods(Static)
        % Determines the width of a microstrip given an impedance, permittivity and distance.
        function w = GetWoverD(z0des, er, d)
            A = z0des ./ 60 .* sqrt((er+1)./2) + (er-1)./(er+1) .* (0.23+0.11./er);
            B = 60.*pi.^2./(z0des.*sqrt(er));
                        
            WoverD1 = 8.*exp(A) ./ (exp(2.*A)-2);
            WoverD2 = 2./pi.*(B-1-log(2.*B-1)+(er-1)./(2.*er).*(log(B-1)+0.39-0.61./er));
            
            % Pozar, page 148, eq 3.197
            % Hong, page 80-81, eq 4.10-4.11
            % Expression is identical in both.
            WoverD = (WoverD1 <= 2) .* WoverD1 + ...
                     (WoverD1 >  2) .* WoverD2;
            
            w = WoverD .* d;
        end
        % Calculates the effective permittivity of the line.
        % When two arguments are provided, the second is expected to be the desired line impedance.
        function ee = EpsilonEffective(er, w_or_Z0, d)
            if(nargin < 3)
                z0 = w_or_Z0;
                d = 1;
                w = Microstrip.GetWoverD(z0, er, d);
            else
                w = w_or_Z0;
            end
            % Pozar, page 148, eq 3.195
            % ee = (er+1)./2+(er-1)./2 .* 1./sqrt(1+12.*d./w);
            % Hong, page 79, eq 4.4
            u = w./d;
            a = 1 + 1/49 .* log((u.^4+(u./52).^2)./(u.^4+0.432)) + 1/18.7 .* log(1+(u./18.1).^3);
            b = 0.564 .* ((er-0.9) ./ (er+3)).^ 0.053;
            ee = (er + 1) ./ 2 + (er - 1) ./ 2 .* (1 + 10./u).^-(a .* b);
            
        end
        % Calculates the attenuation found across a piece of line with the specified parameters.
        function alpha = Attenuation(er, w, d, f, tand)
            ee = Microstrip.EpsilonEffective(er, w, d);
            k0 = 2 .* pi .* f ./ Constants.c0;
            alpha = k0 .* er .* (ee - 1) .* tand ./ (2 .* sqrt(ee) .* (er - 1));
        end
        % Calculates the characteristic impedance of a microstrip.
        function zc = Impedance(er, w, d)
            ee = Microstrip.EpsilonEffective(er, w, d);
            % Pozar, page 148, eq 3.196
            % zc = (w/d <= 1) .* 60./sqrt(ee) .* log(8.*d./w+w./(4.*d)) + ...
            %      (w/d >  1) .* 120*pi./(sqrt(ee) .* (w./d + 1.393 + 0.667 .* log(w./d + 1.444)));
            % Hong, page 79, eq 4.5
            z0 = 120*pi;
            u = w./d;
            F = 6 + (2*pi - 6) .* exp(-(30.666 ./ u) .^ 0.7528);
            zc = z0 ./ (2*pi*sqrt(ee)) .* log(F ./ u + sqrt(1 + (2 ./ u).^2));
        end
    end
end