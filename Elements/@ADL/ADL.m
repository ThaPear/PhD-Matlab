classdef ADL < ADL_Base
    properties
    end
    methods
    end
    methods
        function this = ADL(p, dprev, dnext, sprev, snext, wprev, w, wnext, erhostdown, erhostup, varargin)
            this@ADL_Base(p, dprev, dnext, sprev, snext, wprev, w, wnext, erhostdown, erhostup, varargin{:})
        end
    end
end