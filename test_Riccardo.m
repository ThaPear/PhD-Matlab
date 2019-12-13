

f0=10e9;
l0=3e8/f0;

fs = (0.2:0.01:1)*f0;

dx = 0.5*l0;
dy=dx;
wslot = 0.1*l0;
dslot = 0.1*l0;
slot = Slot(dx, dy, wslot, dslot, FreeSpace(), FreeSpace(), 0);

th = eps * pi/180;
ph = 0;

tc = tic;
Zas = slot.GetInputImpedance(fs, th, ph);

figureex;
    plot(fs, real(Zas));
    plot(fs, imag(Zas), '--');

% Calculate reflection coefficient.
% Gamma = (ZasC - zfeed) ./ (ZasC + zfeed);
% VSWR = (1 + abs(Gamma)) ./ (1 - abs(Gamma));


% Feed optimization
% Algorithm: Trust Region Framework
% Number of evaluations: 921
%               (solver: 920, reloaded: 1)
% Initial goal function value = 8.77176092279  (reloaded)
% Best goal function value    = 8.63950548728
% Last goal function value    = 8.63950548728
% 
% Last solver evaluation time = 00:00:01 h
% 
% 
% Best parameters so far:
% 
%  feed_ms_width = 0.216492
%  feed_shield_distance = 0.7
%  feed_shield_totalangle = 70
%  feed_ms_buried_width = feed_ms_width

% Slot optimization
% Algorithm: Trust Region Framework
% Number of evaluations: 126
%               (solver: 126)
% Initial goal function value = 1.79377870952
% Best goal function value    = 1.16819273587
% Last goal function value    = 1.16819273587
% 
% Last solver evaluation time = 00:00:01 h
% 
% 
% Best parameters so far:
% 
%  slot_feedlength = 0.273897
%  slot_feedwidth = 0.458051
%  slot_width = 0.710773