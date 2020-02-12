function [f, parameters, values] = LoadData(filename)
    if(~exist(filename, 'file'))
        error(['File ''', filename, ''' not found.']);
    end

    if(~contains(filename, '.'))
        error('Please specify file extension.');
    end
    dotindices = find(filename == '.');
    extension = filename(dotindices(end)+1:end);
    switch extension
        case 'csv'
            % readtable method - Gives warnings.
%             data = readtable(filename);
%             f = table2array(data(:,1));
%             values = table2array(data(:,2:end));
%             parameters = [];
%             parameters.labels = data.Properties.VariableNames;
            % csvread method
            data = csvread(filename, 1, 0);
            f = data(:,1);
            values = data(:,2:end);
            
            fileID = fopen(filename);
            titles = fgetl(fileID);
            fclose(fileID);
            
            parameters = [];
            parameters.labels = strsplit(titles, ',');
        otherwise
            error('Unknown extension ''%s''.', extension);
    end
end