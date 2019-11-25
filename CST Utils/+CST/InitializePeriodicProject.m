function [project, dsproject] = InitializePeriodicProject()
    [project, dsproject] = CST.InitializeBasicProject();

    %% Set up boundaries.
    boundary = project.Boundary();
    boundary.AllBoundaries('periodic', 'periodic', ...   % x
                           'periodic', 'periodic', ...   % y
                           'electric', 'expanded open'); % z
end