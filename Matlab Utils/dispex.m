function dispex(varargin)
    % Get the caller function.
    ST = dbstack('-completenames');
    if(length(ST) > 1 && ~isempty(ST(2).file))
        file = ST(2).file;
        file = strrep(file, '\', '/');
        line = ST(2).line;
        href = ['<a href="matlab: opentoline(''', file, ''',', num2str(line), ')"><<</a> '];
%         href = ['<a href="matlab: opentoline(''', file, ''',', num2str(line), ')">', char(55357), char(56599), '</a> ', char(8197)];
    else
        href = '<< ';
    end
    
    % Prepend the link
    msg = varargin{1};
    msg = [href, msg];
    
%     if(nargin == 1)
%         % We're used as a normal disp.
%         disp(msg);
%     else
        % Assume fprintf style.
        fprintf(msg, varargin{2:end});
%     end
end

