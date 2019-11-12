classdef globalsclass
    methods
        function value = subsref(this, S)
            global globalstruct;
            if(ischar(S.subs))
                % If the struct is not initialized, the requested value doesn't exist.
                % If the index doesn't exist in the struct, the requested value doesn't exist.                
                if(~isstruct(globalstruct) || ~isfield(globalstruct, S.subs))
                    error('Attempting to index nonexistent global %s.', S.subs);
                end
                value = globalstruct.(S.subs);
            else
                error('Attempting to index global with nonstring identifier.');
            end
        end
        function this = subsasgn(this, S, value)
            global globalstruct;
            if(ischar(S.subs))
                if(~isstruct(globalstruct))
                    globalstruct = struct();
                end
                globalsstruct.(S.subs) = value;
            else
                error('Attempting to index global with nonstring identifier.');
            end
        end
        function bool = exists(~, name)
            global globalstruct;
            if(~isstruct(globalstruct) || ~isfield(globalstruct, name))
                bool = 0;
                return;
            end
            bool = 1;
            return;
        end
    end
end