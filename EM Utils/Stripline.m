classdef Stripline
    methods(Static)
        function Woverb = GetWoverb(z0, er)
            x = 30.*pi./(sqrt(er).*z0) - 0.441;
            Woverb = (sqrt(er).*z0 <  120) .* x + ...
                     (sqrt(er).*z0 >= 120) .* (0.85 - sqrt(0.6-x));
        end
        function Z0 = CharacteristicImpedance(er, W, b)
            We = W - (W./b < 0.35) .* b.*(W./b - (0.35-W./b).^2);
            Z0 = 30.*pi./sqrt(er) .* b./(We+0.441.*b);
        end
        function alpha = Attenuation(er, W, b, t, Rs, f, tand)
            alphac = Stripline.ConductorAttenuation(er, W, b, t, Rs);
            alphad = Stripline.DielectricAttenuation(er, f, tand);
            alpha = alphac + alphad;
        end
        function alpha = ConductorAttenuation(er, W, b, t, Rs)
            z0 = Stripline.CharacteristicImpedance(er, W, b);
            
            A = 1 + 2.*W ./ (b-t) + 1/pi .* (b+t)./(b-t) .* log((2.*b - t) ./ t);
            B = 1 + b ./ (0.5.*W + 0.7.*t) .* (0.5 + 0.414.*t./W + 1./(2.*pi) .* log(4.*pi.*W ./ t));
            
            alpha = (sqrt(er).*z0 <  120) .* (2.7e-3 .* Rs .* er .* z0) ./ (30.*pi .* (b-t)) .* A + ...
                    (sqrt(er).*z0 >= 120) .* (0.16 .* Rs) ./ (z0 .* b) .* B;
        end
        function alpha = DielectricAttenuation(er, f, tand)
            k = Stripline.PropagationConstant(er, f);
            alpha = k .* tand ./ 2;
        end
        function k = PropagationConstant(er, f)
            k = 2.*pi.*f.*sqrt(er)./Constants.c0;
        end
    end
end