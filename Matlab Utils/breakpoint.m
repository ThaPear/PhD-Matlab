function breakpoint(~, expression)
    ST = dbstack('-completenames');
    if length(ST) < 2; return; end
    if(nargin > 0)
        dbstop('in', ST(2).file, 'at', num2str(ST(2).line+1), 'if', expression);
    else
        dbstop('in', ST(2).file, 'at', num2str(ST(2).line+1));
    end
end