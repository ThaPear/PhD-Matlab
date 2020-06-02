% Set up the path.
if(~exist('IsPathSetup', 'file') || ~IsPathSetup())
    disp('Setting up path.');
    restoredefaultpath;
    matlabrc;
    addpath(                                ...
        'CST-Interface',                    ...
        'Matlab-Utils',                     ...
        'Matlab-Utils/Downloaded',          ...
        'Matlab-Utils/EM Utils',            ...
        'Matlab-Utils/EM Utils/Matrices',   ...
        'PhD-Matlab',                       ...
        'PhD-Matlab/Elements',              ...
        'PhD-Matlab/Elements/Impedances',   ...
        'PhD-Matlab/ Reports/WP_2100&2400', ...
        '');
end

% Parallel pool.
hParallelpool = gcp('nocreate');
if(~isempty(hParallelpool))
    hParallelpool.IdleTimeout = inf;%24*60;
end
clear hParallelpool;

% Pause is error occurs.
dbstop if error;

% Disable new-style zoom stuff.
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'));

% Compact formatting in Command Window.
format compact;

% Close SetupPath if it's open.
try % Not worth erroring over.
    if(isopenineditor('SetupPath.m'))
        closeineditor('SetupPath.m');
    end
catch exception
end