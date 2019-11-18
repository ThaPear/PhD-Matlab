function [f, parameters, values] = Txt(filename)
    fileID = fopen(filename);
    if(fileID == -1)
        error(['File ''', filename, ''' not found.']);
    end
    
    i = 0;
    set = 1;
    reading = 0;
    % Format 1 is a+j*b, format 2 is a*exp(j*b)
    format = 1;
    
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
            if(contains(line, '#"Frequency / GHz"'))
                % Determine format of following data.
                if(strcmp(line(end-4:end), '[Im]"'))
                    format = 1;
                else
                    format = 2;
                end
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
                switch(format)
                    case 1
                        values(i, set) = dat(2) + 1j * dat(3);
                    case 2
                        values(i, set) = dat(2) *exp(1j * dat(3));
                end
            end
        end
    end
    fclose(fileID);
end

function parameters = GetParameters(line)
    % Cut off all except for the parameters.
    line = line(16:end-1);
    
    % Ensure each parameter will be called 'p_parametername' to avoid
    % conflicts with existing variables in this scope.
    parameters = struct();
    line = strrep(line, '; ', '; parameters.');
    line = ['parameters.', line, ';'];
    
    eval(line);
end