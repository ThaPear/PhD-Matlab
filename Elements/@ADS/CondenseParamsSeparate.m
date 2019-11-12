function [ds, ss, ws] = CondenseParamsSeparate(this, ds, ss, ws)
    if(nargin < 2); ds = []; end
    if(nargin < 3); ss = [];  end
    if(nargin < 4); ws = [];  end

    ds = [ds, NaN, this.ds];
    ss = [ss, NaN, this.ss];
    ws = [ws, NaN, this.ws];

    % Assume the last shift specified is used between the ADSs.
    % This is a correct assumption.
    if(length(this.ws) > 1)
        ss = [ss, this.ss(end)];
    end
end
