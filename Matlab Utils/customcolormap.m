% Custom colormap following the given waypoints.
function [clrs] = customcolormap(waypoints, N)
    Nwp = size(waypoints, 1);
    n = [1 1+(1:N)/(N+1)*(Nwp-1), Nwp];
    clrs = [interp1(waypoints(:,1), n); interp1(waypoints(:,2), n); interp1(waypoints(:,3), n)].';
end