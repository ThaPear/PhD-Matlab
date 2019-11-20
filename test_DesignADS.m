close all;
clear;
clear global;
SetupPath;
clc;

f0 = 22e9;
p = 4.4e-3/2;
L = 1.7e-3;
Nlayers = 2;

for(erdes = 3:12)
%     dispex('Testing erdes = %f.\n', erdes);
    ads = DesignADS(f0, p, L, erdes);
    [er, mu] = CST.GetEpsilonMu(ads, f0-1e9, Nlayers);
    if(isempty(er) || isempty(mu)); break; end
    [er2, mu2] = CST.GetEpsilonMu(ads, f0, Nlayers);
    if(isempty(er2) || isempty(mu2)); break; end
    [er3, mu3] = CST.GetEpsilonMu(ads, f0+1e9, Nlayers);
    if(isempty(er3) || isempty(mu3)); break; end
    dispex('erdes = %f, er.x = [%g, %g, %g] at [f0-1, f0 = %g, f0+1]\n', erdes, er.x, er2.x, er3.x, f0);
end

%% Output on 19/11/2019, in vacuum, where dx=p, at f0
% Testing erdes = 3.000000.
% er.x = 2.5841
% er.y = 2.5833
% er.z = 0.94273
% Testing erdes = 4.000000.
% er.x = 3.676
% er.y = 3.6743
% er.z = 0.95735
% Testing erdes = 5.000000.
% er.x = 4.8014
% er.y = 4.8033
% er.z = 1.0177
% Testing erdes = 6.000000.
% er.x = 5.4233
% er.y = 5.4214
% er.z = 1.0317
% Testing erdes = 7.000000.
% er.x = 6.5953
% er.y = 6.5931
% er.z = 1.0618
% Testing erdes = 8.000000.
% er.x = 7.7697
% er.y = 7.7663
% er.z = 1.086
% Testing erdes = 9.000000.
% er.x = 8.9679
% er.y = 8.9695
% er.z = 1.1283
% Testing erdes = 10.000000.
% er.x = 9.1341
% er.y = 9.1346
% er.z = 1.177
% Testing erdes = 11.000000.
% er.x = 10.3823
% er.y = 10.3786
% er.z = 1.2035
% Testing erdes = 12.000000.
% er.x = 11.574
% er.y = 11.5702
% er.z = 1.22

%% Output on 20/11/2019, in vacuum, where dx = 2p, at f0
% Testing erdes = 3.000000.
% er.x = 2.6573
% er.y = 2.6566
% er.z = 0.94805
% Testing erdes = 4.000000.
% er.x = 3.8115
% er.y = 3.8133
% er.z = 1.0033
% Testing erdes = 5.000000.
% er.x = 4.9921
% er.y = 4.9871
% er.z = 1.0606
% Testing erdes = 6.000000.
% er.x = 5.6513
% er.y = 5.6491
% er.z = 1.0615
% Testing erdes = 7.000000.
% er.x = 6.8932
% er.y = 6.9032
% er.z = 1.0558
% Testing erdes = 8.000000.
% er.x = 8.1265
% er.y = 8.1313
% er.z = 1.07
% Testing erdes = 9.000000.
% er.x = 9.3583
% er.y = 9.375
% er.z = 1.188
% Testing erdes = 10.000000.
% er.x = 9.575
% er.y = 9.5758
% er.z = 1.1998
% Testing erdes = 11.000000.
% er.x = 10.9395
% er.y = 10.9399
% er.z = 1.192
% Testing erdes = 12.000000.
% er.x = 12.1917
% er.y = 12.1774
% er.z = 1.2526

%% Output on 20/11/2019, with glue, with dx=2p at f0+-1
% << erdes = 3.000000, er.x = [2.72049, 2.75, 2.76998] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 4.000000, er.x = [3.88604, 3.94318, 3.99547] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 5.000000, er.x = [5.16977, 5.20358, 5.30553] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 6.000000, er.x = [6.34544, 6.53569, 6.67046] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 7.000000, er.x = [7.44542, 8.11682, 7.87675] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 8.000000, er.x = [8.60626, 9.01098, 9.23929] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 9.000000, er.x = [9.83021, 10.0928, 9.96982] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 10.000000, er.x = [11.0803, 11.0598, 9.94992] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 11.000000, er.x = [11.6235, 12.0685, 11.9865] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 12.000000, er.x = [13.0506, 12.8972, 11.7655] at [f0-1, f0 = 2.2e+10, f0+1]