function [parameters, out] = Touchstone(filename)
    switch(exist(filename, 'file'))
        case 0 % Doesn't exist
            error('File ''%s'' does not exist.\n', filename);
        case 2 % It's a file
%             fprintf('Reading Touchstone file ''%s''.\n', filename);
        otherwise
            error('Given filename ''%s'' refers to unknown type.\n', filename);
    end
    
    % Determine number of ports in file.
    % Extract it from the extension specified in the filename.
    extension = strsplit(filename, '.');
    extension = extension{end};
    numPorts = str2double(extension(2:end-1));
    if(isnan(numPorts) || numPorts < 1)
        error('Invalid number of ports extracted.');
    end
        
    
    fileID = fopen(filename);
    line = '!';
    % Extract information from CST header.
    while(line(1) == '!')
        % CST project parameters.
        if(contains(line, '! Parameters = {'))
            % Cut off all except for the parameters.
            line = line(17:end-1);

            parameters = [];
            % Rename each parameter to 'params.parametername' to load them into a struct.
            line = strrep(line, '; ', '; parameters.');
            line = ['parameters.', line, ';'];

            eval(line);
        end
        % Port names.
        if(contains(line, '! Touchstone port ') && ~ contains(line, '! Touchstone port assignment:'))
            % ! Touchstone port 1 = CST MWS port <number> ("<name>")
            split = strsplit(line);
            portI = str2double(split{4});
            portN = str2double(split{9});
            portName = [split{10:end}];
            portName = portName(3:end-2);
            parameters.ports(portI) = struct('number', portN, 'name', portName);
        end
        
        line = fgetl(fileID);
        if(~ischar(line))
            error('Invalid CST TOUCHSTONE file, no data found.');
        end
    end
    % The next line should be the TOUCHSTONE header.
    if(line(1) ~= '#')
        error('Invalid CST TOUCHSTONE file, no TOUCHSTONE header found.');
    end
    
    % Extract information from header line.
    split = strsplit(line);
    header = [];
    header.frequencyunit = split{2};
    header.matrixtype = split{3};
    header.dataformat = split{4};
    header.impedance = str2double(split{6});
    
    % Read data.
    [dat, count] = fscanf(fileID, '%g', [1+2*numPorts^2, inf]);
    parameters.frequencies = dat(1,:);
    % Convert to complex value according to format.
    switch(lower(header.dataformat))
        case 'db' % "DB for dB-angle (dB = 20*log10|magnitude|)"
            dat = 10.^(dat(2:2:end, :)./20) .* exp(1j .* pi/180 .* (dat(3:2:end, :)));
        case 'ma' % "MA for magnitude-angle"
            dat = dat(2:2:end, :) .* exp(1j .* pi/180 .* (dat(3:2:end, :)));
        case 'ri' % "RI for real-imaginary"
            dat = dat(2:2:end, :) + 1j .* dat(3:2:end, :);
        otherwise
            error('Invalid data format %s.', header.dataformat);
    end
    % Reshape into the right shape for output.
    out = reshape(dat, numPorts, numPorts, []);
    
    % Scale frequency values to Hz.
    switch(lower(header.frequencyunit))
        case 'hz'
            parameters.f = parameters.f * 1e0;
        case 'khz'
            parameters.f = parameters.f * 1e3;
        case 'mhz'
            parameters.f = parameters.f * 1e6;
        case 'ghz'
            parameters.f = parameters.f * 1e9;
        otherwise
            error('Unknown frequency unit specified in TOUCHSTONE header.');
    end
end