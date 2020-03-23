function [parameters, varargout] = LoadData(filename)
    split = strsplit(filename, '.');
    extension = split{end};
    
    switch(extension)
        case {'s1p', 's2p', 's3p', 's4p', 's5p', 's6p', 's7p', 's8p', 's9p', 's10p', ...
              'z1p', 'z2p', 'z3p', 'z4p', 'z5p', 'z6p', 'z7p', 'z8p', 'z9p', 'z10p'}
            [parameters, varargout] = CST.LoadData.Touchstone(filename);
        case {'real', 'imag', 'realimag'}
            [parameters, varargout] = CST.LoadData.RealImag(filename);
        case {'txt'}
            [parameters, varargout] = CST.LoadData.Txt(filename);
        case 'mat'
            dat = load(filename);
            parameters = dat;
            varargout = {};
        case 'cut'
            [parameters, varargout] = CST.LoadData.Grasp(filename);
        otherwise
            error('Unknown file format ''%s'' specified.', extension);
    end
    
    if(iscell(varargout) && nargout ~= length(varargout)+1)
        spl = strsplit(filename, '.');
        error(['Invalid number of output arguments for file type ''', spl{end}, ''', ', num2str(length(varargout)+2), ' expected.']);
    elseif(~iscell(varargout))
        error('Non-cell value returned.');
        varargout = {varargout};
    end
end