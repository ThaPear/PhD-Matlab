function [f, value] = RealImag(filename)
    dotindex = strfind(filename, '.');
    filename = filename(1:dotindex(end)-1);
    
    T = readtable([filename, '.real'], 'FileType', 'text');
    valreal = table2array(T(:,2));
    
    T = readtable([filename, '.imag'], 'FileType', 'text');
    valimag = table2array(T(:,2));
    
    f = table2array(T(:,1));
    
    value = valreal + 1j.*valimag;
end