% Set up the path.
if(~exist('IsPathSetup', 'file') || ~IsPathSetup())
    disp('Setting up path.');
    restoredefaultpath;
    matlabrc;
    addpath(                        ...
        'CST Interface',            ...
        'CST Utils',                ...
        'Elements',                 ...
        'Elements/Impedances',      ...
        'EM Utils',                 ...
        'EM Utils/Matrices',        ...
        'HFSS Interface',           ...
        'HFSS Utils',               ...
        'Matlab Utils',             ...
        'Matlab Utils/Downloaded',  ...
        '');
end

% Parallel pool.
hParallelpool = gcp('nocreate');
if(~isempty(hParallelpool))
    hParallelpool.IdleTimeout = 24*60;
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