function [parameters, out] = Grasp(filename)
    fileID = fopen(filename);
    if(fileID == -1)
        error(['File ''', filename, ''' not found.']);
    end
    
    %% Determine the required size of the output array.
    Nlines = 0;
    line = fgetl(fileID);
    % Get the header to determine the number of points in the first dimension.
    format = fscanf(fileID, '%g', [1 7]);
%     V_INI = format(1); % Initial value.
%     V_INC = format(2); % Increment.
    V_NUM = format(3); % Number of values in cut.
%     C     = format(4); % Constant.
%     ICOMP = format(5); % Polarisation control parameter.
    ICUT  = format(6); % Control parameter of cut.
    NCOMP = format(7); % Number of field components.
    
    % Count the number of lines in the file.
    while(ischar(line))
        line = fgetl(fileID);
        Nlines = Nlines+1;
    end
    frewind(fileID);
    
    N(1) = V_NUM;
    N(2) = Nlines / (N(1)+2);
    
    out = zeros(NCOMP, N(1), N(2));
    
    parameters = [];
    parameters.thetas = zeros(1, N(3-ICUT));
    parameters.phis = zeros(1, N(ICUT));
    
    %% Read the data.
    for(i = 1:N(2)) % Loop over the number of blocks in the file.
        header = fgetl(fileID);
        if(~ischar(header))
            break;
        end  % End of file
        
        % First read the header
        %'CST MWS Results: <name> [<number>], Polarization: <Polarization>, FF Origin: (x=<number>, y=<number>, z=<number>)'
        
        % Cut out the origin.
        i1 = strfind(header, '('); i1 = i1(end);
        i2 = strfind(header, ')'); i2 = i2(end);
        str = header(i1+1:i2-1);

        % Parse the origin.
        str = strrep(str, ', ', '; parameters.origin.');
        str = ['parameters.origin.', str, ';'];
        parameters.origin = [];
        try
            eval(str);
        catch(exception)
            warning('GRASP CST header could not be read.');
        end

        % Get the filename, monitor name and polarization.
        str = header(1:i1);
        split = strsplit(str, ',');
        parameters.filename = split{1}(18:end);
        parameters.monitorname = split{2}(2:end);
        parameters.polarization = split{3}(16:end);
        
        % Read header as described in F2.1 in the GRASP manual.
        format = fscanf(fileID, '%g', [1 7]);
        V_INI = format(1); % Initial value.
        V_INC = format(2); % Increment.
        V_NUM = format(3); % Number of values in cut.
        C     = format(4); % Constant.
        ICOMP = format(5); % Polarisation control parameter.
        ICUT  = format(6); % Control parameter of cut.
        NCOMP = format(7); % Number of field components.
        
        if(NCOMP > 2 && ICOMP ~= 1)
            error('Unknown number of components (%i) for component index %i.\nOnly ICOMP=1 supports NCOMP=3.', NCOMP, ICOMP);
        end
        % Determine the appropriate component labels.
        switch(abs(ICOMP))
            case 1
                if(NCOMP == 3)
                    parameters.components = {'Etheta', 'Ephi', 'Erho'};
                else
                    parameters.components = {'Etheta', 'Ephi'};
                end
            case 2
                parameters.components = {'Erhc', 'Elhc'};
            case 3
                parameters.components = {'Eco', 'Ecx'};
            case 4
                parameters.components = {'Emaj', 'Emin'};
            case 5
                parameters.components = {'Etheta/Ephi', 'Ephi/Etheta'};
            case 6
                parameters.components = {'Erhc/Elhc', 'Elhc/Erhc'};
            case 7
                parameters.components = {'Eco/Ecx', 'Ecx/Eco'};
            case 8
                parameters.components = {'Emaj/Emin', 'Emin/Emaj'};
            case 9
                parameters.components = {'Eabs', 'sqrt(Erhc/Elhc)'};
        end
        
        % Determine the correct theta and phi axes.
        switch(ICUT)
            case 1 % A standard polar cut where phi is fixed (C) and theta is varying (V_)
                theta = V_INI + (0:V_NUM-1) * V_INC;
                phi = C;
                parameters.thetas = theta;
                parameters.phis(i) = phi;
            case 2 % A conical cut where theta is fixed (C) and phi is varying (V_).
                phi = V_INI + (0:V_NUM-1) * V_INC;
                theta = C;
                parameters.phis = phi;
                parameters.thetas(i) = theta;
        end
        % Read the next block of data from the file.
        [dat, count] = fscanf(fileID, '%g', [2*NCOMP inf]);
        if(count/(2*NCOMP) ~= V_NUM)
            error('Read incorrect number of data elements (%i/%i)', count, V_NUM);
        end

        % Store data in correct polarization & angle.
        for(k = 1:NCOMP)
            out(k, :, i) = dat(2*k-1, :) + 1j .* dat(2*k, :);
        end
    end
end