classdef Limits
    properties(Constant)
        wmin = 0.1e-3;   % Minimum w is 100 micron.
        wmaxoverp = 0.7; % Maximum w is 0.7 of the unit cell.
        dmin = 0.5e-3;%@(lambda) lambda/10; % Minimum dz
        dmax = @(lambda) lambda/3;  % Maximum dz
        
        % Drill radius, used for the curved corners of the trench.
        drillradius = 0.5e-3 / 2;
        % Minimum distance from trench to slot metal.
        trench2metal = 150e-6;
        
    end
end