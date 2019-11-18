function [f, z] = ReadCST_CopyPaste(filename)
    T = readtable([filename], 'FileType', 'text');
    zreal = table2array(T(:,2));
    zimag = table2array(T(:,3));
    
    f = table2array(T(:,1));
    
    z = zreal + 1j.*zimag;
end