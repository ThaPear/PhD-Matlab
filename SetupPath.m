% Set up the path.
if(~exist('IsPathSetup', 'file') || ~IsPathSetup())
    restoredefaultpath;
    matlabrc;
    addpath(                        ...
        'CST Interface',            ...
        'CST Utils',                ...
        'Elements',                 ...
        'Elements/Impedances',      ...
        'EM Utils',                 ...
        'Matlab Utils',             ...
        'Matlab Utils/Downloaded',  ...
        '');
    disp('Setting up path.');
end

% Parallel pool.
hParallelpool = gcp;
hParallelpool.IdleTimeout = 24*60;
clear hParallelpool;

% Pause is error occurs.
dbstop if error;

% Disable new-style zoom stuff.
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'));
