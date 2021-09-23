
%      o       o
%      |       |
% -   ---     ---
% |   | |     | |
% | L | |     | | er
% |   | |     | |
% -   ---     ---
%      |       |
%      o       o

classdef MicrostripLine < Element
    properties
        L
        z0
        er
        d
        t
        w
        ee
    end
    methods
        % Microstrip line with given length and target impedance.
        % The microstrip is on a dielectric with density er and thickness d.
        % The metal of the microstrip has thickness t.
        function this = MicrostripLine(L, z0, er, d, t)
            if(nargin < 5)
                t = 0;
            end
            this.L = L;
            this.z0 = z0;
            this.er = er;
            this.d = d;
            this.t = t;
            this.w = Microstrip.GetWidth(z0, er, d);
            this.ee = Microstrip.EpsilonEffective(er, this.w, d, t);
        end
        function Z = GetInputImpedance(this, isTE, f, k0, kr)
            error('%s::GetInputImpedance:\n\tInput impedance is not valid on a Line element.\n\tUse a shunt element or an open/shorted/terminated line.', mfilename);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
            % ABCD = GetABCD(this, isTE, f, k0, kr)
%             [kd, ~, ~, kzd] = k(f, obj.er, 0, 0);
            
            kd = k0 .* sqrt(this.ee);
            kzd = -1j .* sqrt(-(kd.^2 - kr.^2));
            zd = this.z0;%Constants.z0 ./ sqrt(this.ee);
            
            % Copied from Elements/Line.m
            % Since kzd can be 0, the sin(kzd*L)/(kzd*L) are combined into
            % sinc such that they result in 1 instead of NaN.
            % Replaced code was kept for clarity.
            % [~, zcte, zctm] = z(this.er, kd, kzd);
            % 
            % if(isTE)
            %     zc = zcte;
            % else
            %     zc = zctm;
            % end
            % 
            % ABCD = ABCDMatrix(                                          ...
            %     cos(kzd.*this.L),                                        ...
            %     1j.*zc.*sin(kzd.*this.L),                                ...
            %     1j./zc.*sin(kzd.*this.L),                                ...
            %     cos(kzd.*this.L)                                         ...
            % );
            % Replacement code below.
            if(isTE)
                ABCD = ABCDMatrix(                                      ...
                    cos(kzd.*this.L),                                    ...
                    1j .* zd .* kd .* this.L .* sinc(kzd.*this.L./pi),    ... % Divide by pi for matlab
                    1j ./ (zd .* kd ./ kzd) .* sin(kzd.*this.L),         ...
                    cos(kzd.*this.L)                                     ...
                );
            else
                ABCD = ABCDMatrix(                                      ...
                    cos(kzd.*this.L),                                    ...
                    1j .* (zd .* kzd ./ kd) .* sin(kzd.*this.L),         ...
                    1j .* kd ./ zd .* this.L .* sinc(kzd.*this.L./pi),    ... % Divide by pi for matlab
                    cos(kzd.*this.L)                                     ...
                );
            end
        end
        function h = GetHeight(this)
            h = this.L;
        end
        function h = GetEffectiveHeight(this, f)
            h = this.L .* sqrt(this.ee);
        end
    end
end