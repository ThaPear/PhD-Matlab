clear;
SetupPath;
close all;

% tc = tic;
% for(i = 1:100000)
%     val = rand();
%     fprintf('val(%i) = %g\n', i, val);
% end
% dt = toc(tc);
% 
% tc = tic;
% for(i = 1:100000)
%     val = rand();
%     dispex('val(%i) = %g\n', i, val);
% end
% dt2 = toc(tc);
% clc;
% disp(['fprintf = ', num2str(dt), 's\n']);
% disp(['dispex = ', num2str(dt2), 's\n']);

z0=50:200;
Z11_00=79.427756+1j*5.5073314;

gamma = @(z0, Z) abs((Z-z0)./(Z+z0));

Z11_31=23.550842+1j*-18.277396;
Z22_31=131.81712+1j*24.106141;
G11_31 = gamma(z0, Z11_31);
G22_31 = gamma(z0, Z22_31);

Z11_28=31.563439+1j*-43.602029;
Z22_28=59.67009+1j*47.430883;
G11_28 = gamma(z0, Z11_28);
G22_28 = gamma(z0, Z22_28);

mx = max([20*log10(G11_31);20*log10(G22_31);20*log10(G11_28);20*log10(G22_28)],[],1);
[m, idx] = min(mx);
m
z0(idx)
G00 = 20*log10(abs((Z11_00-z0(idx))./(Z11_00+z0(idx))))

% words = fileread('words.txt');
% split = strsplit(words, newline);
% split(contains(split, '''')) = [];

aa = tic;

bb = tic;

dt = toc(aa);

dt2 = toc(bb);

% fid = fopen('words_stripped.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for(i = 1:5)
%     tc = tic;
%     project = CST.Application.Active3D();
%     project.Rebuild();
%     fdsolver = project.FDSolver();
%     fdsolver.Start();
%     dt = toc(tc);
%     dispex('dt = %.5g\n', dt);
% end

%% On E-drive
% All 6 cores
% << dt = 110.09
% << dt = 109.03
% << dt = 106.48
% << dt = 106.65
% << dt = 107.71
% 5 cores
% << dt = 107.55
% << dt = 109.75
% << dt = 115.43
% << dt = 111.51
% << dt = 111.59
% 4 cores
% << dt = 114.46
% << dt = 116.65
% << dt = 113.69
% << dt = 115.18
% << dt = 116.75

%% On C-drive
% All 6 cores
% << dt = 108.93
% << dt = 109.71
% << dt = 109.45
% << dt = 108.35
% << dt = 108.8