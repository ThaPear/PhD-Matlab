classdef Cavity
    properties(SetAccess = protected)
        wcavity
        wcavitydiag
    end
    methods
        function this = Cavity(wcavity, wcavitydiag)
            this.wcavity = wcavity;
            this.wcavitydiag = wcavitydiag;
        end
        
    end
end