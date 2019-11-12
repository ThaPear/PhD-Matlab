function printvar(varargin)
% varargin allows the writing of 'printvar c / f' with spaces.
    varname = cell2mat(varargin);
    if(length(varname) < 1)
        return;
    end
    val = evalin('caller', varname);
    if(isnumeric(val))
        if(length(val) > 1)
            for(i = 1:length(val))
                if(isreal(val(i)))
                    disp([varname, '(', num2str(i), ') = ', num2str(val(i))]);
                else
                    disp([varname, '(', num2str(i), ') = ', num2str(val(i)), ' = ', num2str(abs(val(i))), char(8736), num2str(rad2deg(angle(val(i))))]);
                end
%                 varname = repmat(' ', 1, length(varname));
            end
        else
            if(isreal(val))
                disp([varname, ' = ', num2str(val)]);
            else
                disp([varname, ' = ', num2str(val), ' = ', num2str(abs(val)), char(8736), num2str(rad2deg(angle(val)))]);
            end
        end
    elseif(isstruct(val))
        props = fieldnames(val);
        for(i = 1:length(props))
            propname = props{i};
            propval = val.(propname);
            if(isnumeric(propval))
                disp([varname, '.', propname, ' = ', num2str(propval)]);
            else
                disp([mfilename, ': Nonnumeric prop values not supported (', varname, '.', propname, ')']);
            end
        end
    elseif(iscell(val))
        if(isnumeric(val{1}))
            for(i = 1:length(val))
                if(isreal(val{i}))
                    disp([varname, '{', num2str(i), '} = ', num2str(val{i})]);
                else
                    disp([varname, '{', num2str(i), '} = ', num2str(val{i}), ' = ', num2str(abs(val{i})), char(8736), num2str(rad2deg(angle(val{i})))]);
                end
                varname = repmat(' ', length(varname));
            end
        end
    end
end
