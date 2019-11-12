% Custom jet colormap with white in the center.
function [clrs] = customjet(N)
    if(nargin == 0)
        N = 99;
    end
    waypoints = zeros(1,3);
    waypoints( 1,:) = [0.0 0.0 0.5];
    waypoints( 2,:) = [0.0 0.0 1.0];
    waypoints( 3,:) = [0.5 0.5 1.0];
    waypoints( 4,:) = [1.0 1.0 1.0];
    waypoints( 5,:) = [1.0 0.5 0.5];
    waypoints( 6,:) = [1.0 0.0 0.0];
    waypoints( 7,:) = [0.5 0.0 0.0];
    
    clrs = customcolormap(waypoints, N);
end