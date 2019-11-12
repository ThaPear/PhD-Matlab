% ---- resizefig.m

% RESIZEFIG(SzX, SzY, PX, PY, hFig)
% Default is 560 by 420 px.
function resizefig(SzX, SzY, PX, PY, hFig)
    switch(nargin)
        case 2
            hFig = gcf;
            PX = hFig.Position(1);
            PY = hFig.Position(2);
        case 4
            hFig = gcf;
        case 5
            if(isempty(PX))
                PX = hFig.Position(1);
                PY = hFig.Position(2);
            end
    end
    hFig.Position(1) = PX;
    hFig.Position(2) = PY;
    hFig.Position(3) = SzX;
    hFig.Position(4) = SzY;
end