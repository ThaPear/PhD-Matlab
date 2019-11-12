function varargout = plotex(varargin)
    error('Not written');
    
    argout = builtin('plot', varargin{:});
    if(~iscell(argout))
        varargout = {argout};
    else
        varargout = argout;
    end
end