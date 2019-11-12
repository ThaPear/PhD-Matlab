function [value] = GetParameterFromTouchstone(filename, parametername)
    file = fopen(filename);
    line = '';
    % Find the start of the Parameters line.
    while(~contains(line, '! Parameters = {'))
        line = fgetl(file);
    end
    % Cut off all except for the parameters.
    line = line(17:end-1);
    
    % Ensure each parameter will be called 'p_parametername' to avoid
    % conflicts with existing variables in this scope.
    line = strrep(line, '; ', '; p_');
    line = ['p_', line, ';'];
    
    eval(line);
    
    % Get the requested parameter.
    try
        value = eval(['p_', parametername]);
    catch exception
        error('Requested variable ''%s'' not found in ''%s''.', parametername, filename);
    end
end