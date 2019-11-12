function [varargout] = SetNeighbours(this, prevADS, nextADS)
    % [...] = SetNeighbours(this, prevads, nextads)
    % Sets the parameters on the edge elements on the interface
    % between this and prevads, and between this and nextads
    %
    % [...] = SetNeighbours(this, prevads)
    % Sets the parameters on the edge elements on the interface
    % between this and prevads
    %
    % [this, prevads, nextads] = SetNeighbours(...)
    % Sets the parameters on the edge elements of the previous and
    % next ADS as well.
    %
    % [this, prevads] = SetNeighbours(...)
    % Sets the parameters on the edge elements of the previous ADS
    % as well.
    %
    % [this] = SetNeighbours(...)
    % Only sets the parameters on the current ads.
    %
    dprev = 0;
    dnext = 0;
    iprev = 0;
    ifirst = 0;
    ilast = 0;
    inext = 0;

    if(~isempty(prevADS))
        % Find last ADL in previous ADS
        for(i = prevADS.N:-1:1)
            if(contains(class(prevADS.elements{i}), 'ADL')) % Also catches ADL_Real
                iprev = i;
                break;
            else
                dprev = dprev + prevADS.elements{i}.GetHeight();
            end
        end
        % Find first ADL in this ADS
        for(i = 1:this.N)
            if(contains(class(this.elements{i}), 'ADL'))
                ifirst = i;
                break;
            else
                dprev = dprev + this.elements{i}.GetHeight();
            end
        end
        % Last ADL of prev ADS
        prevADS.elements{iprev}.dnext = dprev;
        prevADS.elements{iprev}.wnext = this.elements{ifirst}.w;
        prevADS.elements{iprev}.snext = prevADS.ss(end);
        % First ADL of this ADS
        this.elements{ifirst}.dprev = dprev;
        this.elements{ifirst}.wprev = prevADS.elements{iprev}.w;
        this.elements{ifirst}.sprev = prevADS.ss(end);
    end

    if(nargin == 3 && ~isempty(nextADS))
        % Find last ADL in this ADS
        for(i = this.N:-1:1)
            if(contains(class(this.elements{i}), 'ADL'))
                ilast = i;
                break;
            else
                dnext = dnext + this.elements{i}.GetHeight();
            end
        end
        % Find first ADL in next ADS
        for(i = 1:nextADS.N)
            if(contains(class(nextADS.elements{i}), 'ADL'))
                inext = i;
                break;
            else
                dnext = dnext + nextADS.elements{i}.GetHeight();
            end
        end
        % Last ADL of this ADS
        this.elements{ilast}.dnext = dnext;
        this.elements{ilast}.wnext = nextADS.elements{inext}.w;
        this.elements{ilast}.snext = this.ss(end);
        % First ADL of next ADS
        nextADS.elements{inext}.dprev = dnext;
        nextADS.elements{inext}.wprev = this.elements{ilast}.w;
        nextADS.elements{inext}.sprev = this.ss(end);
    end


    if(nargout == 1)
        varargout = {this};
    elseif(nargout == 2)
        varargout = {this, prevADS};
    else
        varargout = {this, prevADS, nextADS};
    end
end
