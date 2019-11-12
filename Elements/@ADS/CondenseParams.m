function [ds, ss, ws] = CondenseParams(this, ds, ss, ws)
    if(nargin < 2); ds = [0]; end
    if(nargin < 3); ss = [];  end
    if(nargin < 4); ws = [];  end

    ds(end) = ds(end) + this.ds(1);
    ds = [ds, this.ds(2:end)];
    ss = [ss, this.ss];
    ws = [ws, this.ws];

    % Assume the last shift specified is used between the ADSs.
    % This is a correct assumption.
    if(length(this.ws) > 1)
        ss = [ss, this.ss(end)];
    end
end
