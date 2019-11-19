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
    firstline = 1;
    
    f = [];
    while(1)
        line = fgetl(fileID);
        if(line == -1)
            break;
        end
        if(line(1) == '#' || isnan(str2double(line(1:4))))
            if(firstline && line(1) ~= '#')
                % The first line does not start with a '#', so this is an automatically exported
                % file.
                [f, parameters, values] = ReadAutomaticExport(fileID);
            end
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
        % We're no longer on the first line.
        firstline = 0;
    end
    fclose(fileID);
end

function parameters = GetParameters(line)
    % Cut off all except for the parameters.
    line = line(16:end-1);
    
    parameters = [];
    % Rename each parameter to 'params.parametername' to load them into a struct.
    line = strrep(line, '; ', '; parameters.');
    line = ['parameters.', line, ';'];
    
    eval(line);
end

function [f, parameters, values] = ReadAutomaticExport(fileID)
    % Restart the file.
    frewind(fileID);
    
    i = 1;
    f = [];
    values = [];
    while(1)
        % Skip the header.
        line = fgetl(fileID);
        while(line(1) ~= -1 && (isempty(line) || ~strcmp(line(1:10), '----------')))
            line = fgetl(fileID);
        end
        if(line == -1); break; end
        
        % Read the data from the file.
        j = 1;
        line = fgetl(fileID);
        while(~isempty(line))
            % Split line into its separate numbers.
            split = strsplit(line);
            % Delete empty entries.
            split([cellfun(@isempty, split)]) = [];
            
            dat = str2double(split);
            
            f(i, j) = dat(1);
            values(i, j) = dat(2);
            if(length(dat) > 2)
                error('Unknown number of columns in export.');
            end
            
            j = j + 1;
            line = fgetl(fileID);
        end
        i = i + 1;
    end
    
    % Parameters are not exported in the automatic export.
    parameters = [];
end





























