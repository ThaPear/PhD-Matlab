% Set up the path.
if(~exist('IsPathSetup', 'file') || ~IsPathSetup())
    disp('Setting up path.');
    restoredefaultpath;
    matlabrc;
    addpath(                                ...
        'CST-Interface',                    ...
        'Matlab-Utils',                     ...
        'Matlab-Utils/fastintegral',        ...
        'Matlab-Utils/Downloaded',          ...
        'Matlab-Utils/EM Utils',            ...
        'Matlab-Utils/EM Utils/Matrices',   ...
        'PhD-Matlab',                       ...
        'PhD-Matlab/Elements',              ...
        'PhD-Matlab/Elements/Impedances',   ...
        ...'PhD-Matlab/ Reports/WP_2100&2400', ...
        ...'PhD-Matlab/ Reports/WP_2400',      ...
        'PhD-Matlab/ Reports/WP_3100',      ...
        'PhD-Matlab/Tests',                 ...
        'PhD-Matlab/Tests/finitearray',     ...
        'PhD-Matlab/Validations',           ...
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
set(groot, 'DefaultFigureCreateFcn', @(hFig,~) SetFigureProperties(hFig));
set(groot, 'DefaultAxesCreateFcn',   @(hAx,~)  SetAxesProperties(hAx));

% Set some figure/plot defaults.
set(groot, 'DefaultFigureColor', [1 1 1]);
set(groot, 'DefaultAxesLineWidth', 1);
set(groot, 'DefaultLineLineWidth', 1.5);
if ~verLessThan('matlab', '9.10') % 2021a
    set(groot, 'DefaultAxesXLimitMethod', 'tight');
end

% Compact formatting in Command Window.
format compact;

% Close SetupPath if it's open.
try % Not worth erroring over.
    if(isopenineditor('SetupPath.m'))
        closeineditor('SetupPath.m');
    end
catch exception
end

function SetFigureProperties(hFig)
    addToolbarExplorationButtons(hFig);
end

function SetAxesProperties(hAx)
    if(isa(hAx, 'matlab.ui.control.UIAxes'))
        return;
    end
    set(hAx.Toolbar,'Visible','off')
end