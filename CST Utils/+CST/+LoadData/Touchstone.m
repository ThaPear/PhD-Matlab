function [f, parameters, out] = Touchstone(filename)
    switch(exist(filename, 'file'))
        case 0 % Doesn't exist
            disp(['File "', filename, '" does not exist']);
            f = [];
            parameters = [];
            for(i = 1:nargout)
                out{i} = [];
            end
            return;
        case 2 % It's a file
            disp(['Reading Touchstone file ''', filename, '''.']);
        otherwise
            error(['Given filename ("', filename, '") refers to unknown type']);
    end
    
    if(filename(end-2) == 's')
        [f, out] = ReadTouchstone_S(filename);
    else
        [f, out] = ReadTouchstone_Z(filename);
    end
    
    parameters = ReadParameters(filename);
end

function [f, out] = ReadTouchstone_S(filename)
    data = sparameters(filename);
    f = data.Frequencies;

    switch(data.NumPorts)
        case 1
            S = [];
            S.s11 = squeeze(data.Parameters(1,1,:));

            out{1} = S;
        case 2
            S = [];
            S.s11 = squeeze(data.Parameters(1,1,:));
            S.s12 = squeeze(data.Parameters(1,2,:));
            S.s21 = squeeze(data.Parameters(2,1,:));
            S.s22 = squeeze(data.Parameters(2,2,:));
            
            out{1} = S;
        case 4
            % TE
            Ste = [];
            Ste.s11 = squeeze(data.Parameters(1,1,:));
            Ste.s12 = squeeze(data.Parameters(1,3,:));
            Ste.s21 = squeeze(data.Parameters(3,1,:));
            Ste.s22 = squeeze(data.Parameters(3,3,:));

            % TM
            Stm = [];
            Stm.s11 = squeeze(data.Parameters(2,2,:));
            Stm.s12 = squeeze(data.Parameters(2,4,:));
            Stm.s21 = squeeze(data.Parameters(4,2,:));
            Stm.s22 = squeeze(data.Parameters(4,4,:));

            out{1} = Ste;
            out{2} = Stm;
        otherwise
            out{1} = data.Parameters;
%             error(['Unknown number of ports in S-parameters: ', num2str(data.NumPorts)]);
    end
end

function [f, out] = ReadTouchstone_Z(filename)
    data = zparameters(filename);
    f = data.Frequencies;

    switch(data.NumPorts)
        case 1
            out{1} = squeeze(data.Parameters(1,1,:));
        otherwise
            out{1} = data.Parameters;
%             error(['Unknown number of ports in Z-parameters: ', num2str(data.NumPorts)]);
    end
end

function parameters = ReadParameters(filename, parametername)
    file = fopen(filename);
    line = '';
    % Find the start of the Parameters line.
    while(~contains(line, '! Parameters = {'))
        line = fgetl(file);
    end
    % Close the file.
    fclose(file);
    
    % Cut off all except for the parameters.
    line = line(17:end-1);
    
    parameters = [];
    % Rename each parameter to 'params.parametername' to load them into a struct.
    line = strrep(line, '; ', '; parameters.');
    line = ['parameters.', line, ';'];
    
    eval(line);
end