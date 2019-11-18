function [f, z] = ReadCST_RealImag(filename)
    T = readtable([filename, '.real'], 'FileType', 'text');
    zreal = table2array(T(:,2));
    
    T = readtable([filename, '.imag'], 'FileType', 'text');
    zimag = table2array(T(:,2));
    
    f = table2array(T(:,1));
    
    z = zreal + 1j.*zimag;
end