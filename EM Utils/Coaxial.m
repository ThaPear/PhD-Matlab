classdef Coaxial
    % Equations based on: https://www3.nd.edu/~wzech/CoaxialTransmissionLine.pdf
    methods(Static)
        function z = Impedance(outerR, innerR, er, mur)
            if(nargin < 4)
                mur = 1;
            end
%             C = Coaxial.CapacitancePerUnitLength(outerR, innerR, er, mur);
%             L = Coaxial.InductancePerUnitLength(outerR, innerR, er, mur);
%             z = sqrt(L/C);
            % Expanded the above.
            mu = Constants.mu0 * mur;
            ep = Constants.ep0 * er;
            z = 1 / (2*pi) * sqrt(mu / ep) * log(outerR / innerR);
        end
        function ratio = GetRadiusRatio(z0, er, mur)
            if(nargin < 4)
                mur = 1;
            end
            % Inverted the above equation.
            mu = Constants.mu0 * mur;
            ep = Constants.ep0 * er;
            ratio = exp(2*pi*z0*sqrt(ep/mu));
        end
        function C = CapacitancePerUnitLength(outerR, innerR, er, mur)
%             if(nargin < 4)
%                 mur = 1;
%             end
            C = 2*pi*Constants.ep0 * er / log(outerR/innerR);
        end
        function L = InductancePerUnitLength(outerR, innerR, er, mur)
            if(nargin < 4)
                mur = 1;
            end
            L = Constants.mu0 * mur / (2*pi) * log(outerR/innerR);
        end
        function f = Cutoff(outerR, innerR, er, mur)
%             if(nargin < 4)
%                 mur = 1;
%             end
            kc = 2 / (outerR + innerR);
            f = Constants.c0 * kc / (2*pi*sqrt(er));
        end
    end
end