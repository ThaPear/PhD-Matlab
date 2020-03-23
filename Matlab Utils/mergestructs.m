function [struct1] = mergestructs(struct1, struct2)
    fname = fieldnames(struct2);
    for i = 1:length(fname)
        % Check if field already exists in struct1
        if(isfield(struct1, fname(i)))
            % Check if the value is different
            if(struct1.(fname{i}) ~= struct2.(fname{i}))
                error('Field exists in both structs with different values');
            end
        end
        
        struct1.(fname{i}) = struct2.(fname{i});
    end
end