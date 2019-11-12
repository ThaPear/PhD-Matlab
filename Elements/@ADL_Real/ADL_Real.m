classdef ADL_Real < ADL_Base
    properties
    end
    methods
    end
    methods
        function this = ADL_Real(p, dprev, dnext, sprev, snext, wprev, w, wnext, erhostdown, erhostup)
            %  |--     ------     -----     -----     -----  |  -
            %  | |wnext|          :                          |  |
            %  |                  |-|<snext         erhostup |  | dnext
            %  |                    :                        |  |
            %  |-----     -----     -----     -----     -----|  -
            %  |    |- w -|         :                        |  |
            %  |               |----|<sprev       erhostdown |  | dprev
            %  |               :                             |  |
            %  |    ------     -----     -----     -----     |  -
            %            |wprev|
            
            % Function is a fit through the points determined using the
            % results obtained in fab_adlglue.m.
% 127 um thick layers
%                    a =      -4.739215763446925; % (-26.57, 17.09)
%                    b =       0.116810429067856; % (-0.8618, 1.095)
%                    c =       1.607814456799277; % (-32.73, 35.95)
%                    d =       0.267814456799277; % (-0.9875, 1.524)
%                    e =       4.423492475147106; % (-8.115, 16.96)

% 25.4 um thick layers
                   a =     0.306018752897363;%  (-0.3667, 0.9787)
                   b =    -0.350282031462144;%  (-0.6335, -0.06705)
                   c =     0.545673344874803;%  (-5.821, 6.912)
                   d =     0.084263850068818;%  (-0.6258, 0.7943)
                   e =     0.312746876442275;%  (-6.726, 7.351)
            erfunc = @(w) a.*w.^b+c.*w.^d+e;
            
            ernext = erfunc(w/p);
            erprev = erfunc(w/p);
            this@ADL_Base(p, dprev, dnext, sprev, snext, wprev, w, wnext, ernext, erprev);
        end
    end
end