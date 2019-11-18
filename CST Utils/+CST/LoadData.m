function [f, parameters, varargout] = LoadData(filename)
    split = strsplit(filename, '.');
    extension = split{end};
    
    switch(extension)
        case {'s1p', 's2p', 's3p', 's4p', 's5p', 's6p', 's7p', 's8p', 's9p', 's10p', ...
              'z1p', 'z2p', 'z3p', 'z4p', 'z5p', 'z6p', 'z7p', 'z8p', 'z9p', 'z10p'}
            [f, parameters, varargout] = CST.LoadData.Touchstone(filename);
        case {'real', 'imag', 'realimag'}
            parameters = [];
            [f, varargout] = CST.LoadData.RealImag(filename);
        case {'txt'}
            [f, parameters, varargout] = CST.LoadData.Txt(filename);
        case 'mat'
            dat = load(filename);
            if(~isfield(dat, 'f'))
                error('No frequency found in ''%s''.', filename);
            end
            f = dat.f;
            parameters = dat;
            parameters.f = [];
            varargout = {};
        otherwise
            error('Unknown file format ''%s'' specified.', extension);
    end
    
    if(iscell(varargout) && nargout ~= length(varargout)+2)
        spl = strsplit(filename, '.');
        error(['Invalid number of output arguments for file type ''', spl{end}, ''', ', num2str(length(varargout)+2), ' expected.']);
    else
        varargout = {varargout};
    end
end