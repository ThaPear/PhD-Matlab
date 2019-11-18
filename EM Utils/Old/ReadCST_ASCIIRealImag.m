function [f, values, parameters] = ReadCST_ASCIIRealImag(filename)
    fileID = fopen(filename);
    if(fileID == -1)
        error(['File ''', filename, ''' not found.']);
    end
    
    i = 0;
    set = 1;
    reading = 0;
    
    f = [];
    while(1)
        line = fgetl(fileID);
        if(line == -1)
            break;
        end
        if(line(1) == '#' || isnan(str2double(line(1:4))))
            if(reading)
                reading = 0;
                set = set + 1;
            end
            if(contains(line, '#Parameters = {'))
                parameters{set} = GetParameters(line);
            end
        else
            if(~reading)
                i = 0;
                reading = 1;
            end
            i = i + 1;
            dat = str2double(strsplit(line));
            if(isnan(dat(1)))
                dat = dat(2:end);
            end
            f(i, set) = dat(1);
            if(length(dat) == 2)     % Frequency & Real value
                values(i, set) = dat(2);
            elseif(length(dat) >= 3) % Frequency & Complex value
                values(i, set) = dat(2) + 1j * dat(3);
            end
        end
    end
    fclose(fileID);
end

function p = GetParameters(line)
    % Cut off all except for the parameters.
    line = line(16:end-1);
    
    % Ensure each parameter will be called 'p_parametername' to avoid
    % conflicts with existing variables in this scope.
    p = struct();
    line = strrep(line, '; ', '; p.');
    line = ['p.', line, ';'];
    
    eval(line);
end