classdef ADL_Base < Element
    properties
        p
        dprev
        dnext
        sprev
        snext
        wprev
        w
        wnext
        
        ernext
        erprev
    end
    methods
        ABCD = GetABCD(this, isTE, f, k0, kr)
        Bn = GetLayerSusceptance(this, f)
    end
    methods
        function this = ADL_Base(p, dprev, dnext, sprev, snext, wprev, w, wnext, erhostdown, erhostup)
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
            
            this.p = p;
            this.dprev = dprev;
            this.dnext = dnext;
            this.sprev = sprev;
            this.snext = snext;
            this.wprev = wprev;
            this.w = w;
            this.wnext = wnext;
            
            this.ernext = erhostup;
            this.erprev = erhostdown;
        end
        function flippedadl = Flip(this)
            % Flip all references to 'next' and 'prev'.
            flippedadl = ADL_Base(this.p, this.dnext, this.dprev, ...
                                          this.snext, this.sprev, ...
                                          this.wnext, this.w, this.wprev, ...
                                          this.erhostup, this.erhostdown);
        end
    end
end