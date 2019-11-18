function [f, S] = ReadCST(filename)
    switch(exist(filename, 'file'))
        case 0 % Doesn't exist
            warning(['File "', filename, '" does not exist']);
            f = [];
            S = [];
        case 2 % It's a file
            T = readtable(filename);
            f = table2array(T(:,1));
            f = 1e9 * f.';
            S = table2array(T(:,2)) + 1j * table2array(T(:,3));
            S = S.';
        case 7 % It's a folder
            S = [];
            f = [];
            if(isfile([filename, '/s11.txt']))
                T = readtable([filename, '/s11.txt']);
                f = table2array(T(:,1));
                S.s11 = table2array(T(:,2)) + 1j * table2array(T(:,3));
                S.s11 = S.s11.';
            end
            if(isfile([filename, '/s12.txt']))
                T = readtable([filename, '/s12.txt']);
                f = table2array(T(:,1));
                S.s12 = table2array(T(:,2)) + 1j * table2array(T(:,3));
                S.s12 = S.s12.';
            end
            if(isfile([filename, '/s21.txt']))
                T = readtable([filename, '/s21.txt']);
                f = table2array(T(:,1));
                S.s21 = table2array(T(:,2)) + 1j * table2array(T(:,3));
                S.s21 = S.s21.';
            end
            if(isfile([filename, '/s22.txt']))
                T = readtable([filename, '/s22.txt']);
                f = table2array(T(:,1));
                S.s22 = table2array(T(:,2)) + 1j * table2array(T(:,3));
                S.s22 = S.s22.';
            end
            % Make the array 1xN instead of Nx1.
            f = 1e9 * f.';
        otherwise
            warning(['Argument filename ("', filename, '") refers to unknown type']);
            f = [];
            S = [];
    end
end