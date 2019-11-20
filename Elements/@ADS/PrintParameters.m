function str = PrintParameters(this, f0)
    lambda0 = Constants.c0 / f0;
    [epsilon, ~] = this.GetEpsilonMu(f0, 0, 0);
    er = epsilon.x;
    lambda = lambda0 / sqrt(er);

    str = '';
    if(nargout == 0)
        delim = '\n';
    else
        delim = ' - ';
    end
    % TODO: Fix erhost when varying per-layer.
    if(~isempty(this.erhosts))
        str = [str, [char(949), 'rhost = ', num2str(this.erhosts(1))], ', '];
    end
    str = [str, [char(949), 'r = ', num2str(er)], delim];
    str = [str, VarStr('p', this.p, lambda0, lambda, delim), delim];
    str = [str, VarStr('ds', this.ds, lambda0, lambda, delim), delim];
    str = [str, VarStr('ss', this.ss, lambda0, lambda, delim), delim];
    str = [str, VarStr('ws', this.ws, lambda0, lambda, delim), delim];
    str = [str, VarStr('h', this.GetHeight(), lambda0, lambda, delim), delim];
    if(nargout == 0)
        dispex(str);
    end
    function str = VarStr(name, values, lambda0, lambda, delim)
        if(all(values == values(1)))
            value = values(1);
            str = [name, ' = ', num2str(value*1e3, '%.5f'), 'mm = ', num2str(value/lambda0, '%.5f'), char(955), '0 = ', num2str(value/lambda, '%.5f'), char(955), 'eff'];
        else % Values are different for each layer.
            mmstr = [name, ' = ['];
            l0str = [repmat(' ', 1, length(name)*2+2), ' = ['];
            lestr = [repmat(' ', 1, length(name)*2+2), ' = ['];
            for(i = 1:length(values))
                value = values(i);
                mmstr = [mmstr, num2str(value*1e3, '%.5f'), ', '];
                l0str = [l0str, num2str(value/lambda0, '%.5f'), ', '];
                lestr = [lestr, num2str(value/lambda, '%.5f'), ', '];
            end
            mmstr = [mmstr, '] mm'];
            l0str = [l0str, '] ', char(955), '0'];
            lestr = [lestr, '] ', char(955), 'eff'];
            str = [mmstr, delim, l0str, delim, lestr];
        end
    end
end
