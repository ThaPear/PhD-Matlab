function [parameters, out] = Txt(filename)
    fileID = fopen(filename);
    if(fileID == -1)
        error(['File ''', filename, ''' not found.']);
    end
    
    format = DetermineFormat(fileID);
    
    switch(format)
        case 1
            [parameters, out] = ReadAutomaticExport(fileID);
        case 2
            [f, parameters, out] = ReadManualExport(fileID);
        case 3
            [f, parameters, out] = ReadCopyPaste(fileID);
    end
    % Wrap up the values in a cell since they're returned as varargout in LoadData.
    out = {out};
end

% Formats:
% 1 - Automatically exported
% 2 - Manually exported
% 3 - Copy-pasted
function format = DetermineFormat(fileID)
    line = fgetl(fileID);
    line2 = fgetl(fileID);
    format = 0;
    
    split = strsplit(line);
    split([cellfun(@isempty, split)]) = []; % Remove any empty values.
    
    if(~isnan(str2double(split{1})) || length(line2) > 9 && strcmp(line2(1:10), '----------'))
        format = 1;
    end
    if(contains(line, '#Parameters = {'))
        format = 2;
    end
    if(contains(line, 'CST XY Data Exchange Format V2'))
        format = 3;
    end
    if(format == 0)
        error('Unrecognized format provided.');
    end
    
    % Restart the file.
    frewind(fileID);
end

function [parameters, out] = ReadAutomaticExport(fileID)
    % Check if the file has a header, or if it starts with numbers right away.
    line = fgetl(fileID);
    split = strsplit(line);
    split([cellfun(@isempty, split)]) = [];
    if(isnan(str2double(split{1})))
        % There is a header.
        % Determine axis labels.
        line = replace(line, '   ', '  ');      % Replace any triple space with a double one.
        split = strsplit(line, '  ');           % Split the string on double spaces.
        split([cellfun(@isempty, split)]) = []; % Remove any empty values.
                                                % The string is now split into sections separated by at
                                                % least 2 spaces.
        parameters{1}.labels = split;
        line = fgetl(fileID); % Get the ------------ line.
        line = fgetl(fileID); % Get the first data line.
    else
        parameters{1}.labels = {};
    end
    
    % Determine the number of data values per line.
    split = strsplit(line);
    split([cellfun(@isempty, split)]) = [];
    nValsPerLine = length(split);
    
    % Go back 1 line to ensure all data is read.
    fseek(fileID, -(length(line)+2), 0);
    [dat, count] = fscanf(fileID, '%g', [nValsPerLine, inf]);
    
    out{1} = dat;
    
    
%     i = 1;
%     f = [];
%     out = [];
%     while(1)
%         % Skip the header.
%         line = fgetl(fileID);
%         while(line(1) ~= -1 && (isempty(line) || ~strcmp(line(1:10), '----------')))
%             line = fgetl(fileID);
%         end
%         if(~ischar(line)); break; end
%         
%         % Read the data from the file.
%         j = 1;
%         line = fgetl(fileID);
%         while(~isempty(line))
%             % Split line into its separate numbers.
%             split = strsplit(line);
%             % Delete empty entries.
%             split([cellfun(@isempty, split)]) = [];
%             
%             dat = str2double(split);
%             
%             f{i}(j) = dat(1);
%             out{i}(j) = dat(2);
%             if(length(dat) > 2)
%                 error('Unknown number of columns in export.');
%             end
%             
%             j = j + 1;
%             line = fgetl(fileID);
%         end
%         i = i + 1;
%     end
%     
%     % Parameters are not exported in the automatic export.
%     parameters = [];
end

function [f, parameters, values] = ReadManualExport(fileID)
    i = 0;
    reading = 0;
    format = 1;
    set = 1;
    
    f = [];
    while(1)
        line = fgetl(fileID);
        if(~ischar(line))
            break;
        end
        if(line(1) == '#' || isnan(str2double(line(1:4))))
            if(reading)
                reading = 0;
                set = set + 1;
            end
            if(contains(line, '#Parameters = {'))
                parameters{set} = GetParametersForManualExport(line);
            end
            if(length(line >= 2) && strcmp(line(1:2), '#"'))%contains(line, '#"Frequency / GHz"'))
                % Determine format of following data.
                if(strcmp(line(end-4:end), '[Im]"'))
                    format = 1;
                else
                    format = 2;
                end
                split = strsplit(line, '"');
                split = split(2:2:end);
                parameters{set}.labels = split;
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
            % Assume the first column is frequency.
            f{set, i} = dat(1);
            if(length(dat) == 2)     % Frequency & Real value
                values{set}(i) = dat(2);
            elseif(length(dat) >= 3) % Frequency & Complex value
                % TODO: Support multiple columns of data.
                % Such as #"Frequency / GHz"	"S1,1 (12) [Re]"	"S1,1 (12) [Im]"	"Ref.Imp. [Re]"	"Ref.Imp. [Im]"
                switch(format)
                    case 1
                        values{set}(i) = dat(2) + 1j * dat(3);
                    case 2
                        values{set}(i) = dat(2) *exp(1j * dat(3));
                end
            end
        end
    end
end

function parameters = GetParametersForManualExport(line)
    % Cut off all except for the parameters.
    line = line(16:end-1);
    
    parameters = [];
    % Rename each parameter to 'params.parametername' to load them into a struct.
    line = strrep(line, '; ', '; parameters.');
    line = ['parameters.', line, ';'];
    
    eval(line);
end

function [f, parameters, values] = ReadCopyPaste(fileID)
    set = 1;
    reading = 0;
    i = 0;
    while(1)
        line = fgetl(fileID);
        if(~ischar(line))
            break;
        end
        if(~reading)
            [parameters{set}, format] = ReadCopyPasteHeader(fileID); %#ok<AGROW>
            reading = 1;
        else
            %% Read data
            split = strsplit(line, char(9));
            split(cellfun(@isempty,split)) = [];
            if(isempty(split))
                reading = 0;
                i = 0;
                set = set + 1;
                continue;
            end
            i = i + 1;
            
            data = str2double(split);
            f{set}(i) = data(1);
%             if(length(data) > 2)
%                 error('Unknown number of elements (%i) in data.', length(data));
%             end
            values{set}(i, :) = data(2:end);
        end
    end
end

function [parameters, format] = ReadCopyPasteHeader(fileID)
    while(1)
        %% Read header
        line = fgetl(fileID);
        if(~ischar(line))
            break;
        end
        % Replace tabs with spaces.
%         linerep = strrep(line, char(9), ' ');
        if(~isnan(str2double(strsplit(line, char(9)))))
            % Reading numbers, so go back 1 line and return.
            fseek(fileID, -(length(line)+2), 0);
            return;
        end
        
        split = strtrim(strsplit(line, '='));
        split(cellfun(@isempty,split)) = [];
        
        if(isempty(split))
            continue;
        end
        
        parameters = struct('labels', '');
        switch(split{1})
            case 'Title'
                parameters.title = [split{2:end}];
            case 'Xlabel'
                parameters.labels = [parameters.labels, [split{2:end}]];
            case 'Ylabel'
                parameters.labels = [parameters.labels, [split{2:end}]];
            case 'Parameters'
                params = GetParametersForManualExport(line);
                parameters = mergestructs(parameters, params);
            case 'Result type'
                switch(split{2})
                    case 'real:cartesian'
                    case 'complex:none;impedance'
                end
            case 'View type'
                switch(split{2})
                    case 'Cartesian'
                        format = 1;
                    case 'Polar'
                        format = 2;
                    otherwise
                        error('Unknown format in copy-paste file.');
                end
            otherwise
        end

%         breakpoint;
    end
    a=1;
end








































