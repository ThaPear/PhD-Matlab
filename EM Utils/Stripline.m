classdef Stripline
    % -------------------------------------
    %                                / \
    %                                 |
    %             -------------       | b
    %             <----------->       |
    %                   w            \ /
    % -------------------------------------
    methods(Static)
        function w = GetWidth(z0des, er, b)
            x = 30.*pi./(sqrt(er).*z0des) - 0.441;
            u = (sqrt(er).*z0des <  120) .* x + ...
                (sqrt(er).*z0des >= 120) .* (0.85 - sqrt(0.6-x));
            w = b .* u;
        end
        function Z0 = CharacteristicImpedance(er, w, b)
            We = w - ((w./b >= 0.35) .* 0 + ...
                      (w./b <  0.35) .* b.*(w./b - (0.35-w./b).^2));
            Z0 = 30.*pi./sqrt(er) .* b./(We+0.441.*b);
        end
        function alpha = Attenuation(er, w, b, t, Rs, f, tand)
            alphac = Stripline.ConductorAttenuation(er, w, b, t, Rs);
            alphad = Stripline.DielectricAttenuation(er, f, tand);
            alpha = alphac + alphad;
        end
        function alpha = ConductorAttenuation(er, w, b, t, Rs)
            z0 = Stripline.CharacteristicImpedance(er, w, b);
            
            A = 1 + 2.*w ./ (b-t) + 1/pi .* (b+t)./(b-t) .* log((2.*b - t) ./ t);
            B = 1 + b ./ (0.5.*w + 0.7.*t) .* (0.5 + 0.414.*t./w + 1./(2.*pi) .* log(4.*pi.*w ./ t));
            
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