close all;
clear;
clear global;
SetupPath;
clc;

f0 = 22e9;
p = 4.4e-3/2;
L = 1.7e-3;
Nlayers = 1;

dispex('Starting ADS.\n');

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

dispex('Starting ADS_Real.\n');

for(erdes = 3:12)
%     dispex('Testing erdes = %f.\n', erdes);
    ads = DesignADS_Real(f0, p, L, erdes);
    [er, mu] = CST.GetEpsilonMu(ads, f0-1e9, Nlayers);
    if(isempty(er) || isempty(mu)); break; end
    [er2, mu2] = CST.GetEpsilonMu(ads, f0, Nlayers);
    if(isempty(er2) || isempty(mu2)); break; end
    [er3, mu3] = CST.GetEpsilonMu(ads, f0+1e9, Nlayers);
    if(isempty(er3) || isempty(mu3)); break; end
    dispex('erdes = %f, er.x = [%g, %g, %g] at [f0-1, f0 = %g, f0+1]\n', erdes, er.x, er2.x, er3.x, f0);
end

%% Output on 19/11/2019, in vacuum, where dx=p, at f0, with double first-layer spacing
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

%% Output on 20/11/2019, in vacuum, where dx = 2p, at f0, with double first-layer spacing
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

%% Output on 20/11/2019, with glue, with dx=2p at f0+-1, with double first-layer spacing
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

%% Output on 20/11/2019, in vacuum, with dx=2p at f0+-1, with Nlayers=2
% << erdes = 3.000000, er.x = [2.76795, 2.80709, 2.85876] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 4.000000, er.x = [4.38418, 5.25616, -0.0219256] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 5.000000, er.x = [3.15642, 3.64159, 3.48491] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 6.000000, er.x = [4.23027, 3.89266, 3.54495] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 7.000000, er.x = [4.01246, 3.63077, 3.2561] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 8.000000, er.x = [3.97158, 3.50154, 3.06543] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 9.000000, er.x = [3.79799, 3.28362, 2.79396] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 10.000000, er.x = [3.64577, 3.12537, 2.63755] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 11.000000, er.x = [3.42699, 2.85791, 2.33336] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 12.000000, er.x = [3.16889, 2.57804, 2.01308] at [f0-1, f0 = 2.2e+10, f0+1]

%% Output on 20/11/2019, in vacuum, with dx=2p at f0+-1, with Nlayers=1
% << Starting ADS.
% << erdes = 3.000000, er.x = [2.64178, 2.65519, 2.67127] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 4.000000, er.x = [3.79181, 3.82986, 3.8587] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 5.000000, er.x = [4.92946, 4.99782, 5.04813] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 6.000000, er.x = [6.08121, 6.16544, 6.27431] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 7.000000, er.x = [6.87045, 6.95207, 7.04165] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 8.000000, er.x = [8.07973, 8.18417, 8.33495] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 9.000000, er.x = [9.25939, 9.44644, 9.68317] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 10.000000, er.x = [9.67165, 9.87466, 10.098] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 11.000000, er.x = [10.8541, 11.156, 11.6078] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 12.000000, er.x = [12.0389, 12.516, 13.2701] at [f0-1, f0 = 2.2e+10, f0+1]
% << Starting ADS_Real.
% << erdes = 3.000000, er.x = [2.71914, 2.7316, 2.74575] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 4.000000, er.x = [3.89328, 3.93511, 3.95088] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 5.000000, er.x = [5.13972, 5.16856, 5.23723] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 6.000000, er.x = [6.42981, 6.51284, 6.58582] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 7.000000, er.x = [7.69582, 7.7966, 7.92706] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 8.000000, er.x = [8.79582, 8.91025, 8.86834] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 9.000000, er.x = [9.81502, 9.84872, 10.0067] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 10.000000, er.x = [10.8731, 11.0785, 11.0556] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 11.000000, er.x = [11.7516, 12.178, 12.2438] at [f0-1, f0 = 2.2e+10, f0+1]
% << erdes = 12.000000, er.x = [13.157, 13.3492, 13.6951] at [f0-1, f0 = 2.2e+10, f0+1]