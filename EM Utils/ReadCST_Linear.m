function [f, z] = ReadCST_RealImag(filename)
    T = readtable([filename, '.real'], 'FileType', 'text');
    f = table2array(T(:,1));
    z = table2array(T(:,2));
    
end