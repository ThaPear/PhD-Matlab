classdef Globals
    properties(Constant)
        slot_s0 = 0.25;
    end
    methods(Static)
        function bool = exists(name)
            bool = isprop(Globals, name);
        end
    end
end
% function glob = Globals
%     glob = globalsclass;
% end
