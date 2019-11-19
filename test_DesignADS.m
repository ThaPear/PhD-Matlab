
f0 = 22e9;
p = 4.4e-3/2;
L = 1.7e-3;
erdes = 8;
for(erdes = 3:12)
    fprintf('Testing erdes = %f.\n', erdes);
    ads = DesignADS(f0, p, L, erdes);
    [er, mu] = CST.GetEpsilonMu(ads, f0);
    printvar er;
end

%% Output on 19/11/2019
% Testing erdes = 3.000000.
% Reading Touchstone file 'Temp/S0.s4p'.
% Reading Touchstone file 'Temp/S.s4p'.
% er.x = 2.5841
% er.y = 2.5833
% er.z = 0.94273
% Testing erdes = 4.000000.
% Reading Touchstone file 'Temp/S0.s4p'.
% Reading Touchstone file 'Temp/S.s4p'.
% er.x = 3.676
% er.y = 3.6743
% er.z = 0.95735
% Testing erdes = 5.000000.
% Reading Touchstone file 'Temp/S0.s4p'.
% Reading Touchstone file 'Temp/S.s4p'.
% er.x = 4.8014
% er.y = 4.8033
% er.z = 1.0177
% Testing erdes = 6.000000.
% Reading Touchstone file 'Temp/S0.s4p'.
% Reading Touchstone file 'Temp/S.s4p'.
% er.x = 5.4233
% er.y = 5.4214
% er.z = 1.0317
% Testing erdes = 7.000000.
% Reading Touchstone file 'Temp/S0.s4p'.
% Reading Touchstone file 'Temp/S.s4p'.
% er.x = 6.5953
% er.y = 6.5931
% er.z = 1.0618
% Testing erdes = 8.000000.
% Reading Touchstone file 'Temp/S0.s4p'.
% Reading Touchstone file 'Temp/S.s4p'.
% er.x = 7.7697
% er.y = 7.7663
% er.z = 1.086
% Testing erdes = 9.000000.
% Reading Touchstone file 'Temp/S0.s4p'.
% Reading Touchstone file 'Temp/S.s4p'.
% er.x = 8.9679
% er.y = 8.9695
% er.z = 1.1283
% Testing erdes = 10.000000.
% Reading Touchstone file 'Temp/S0.s4p'.
% Reading Touchstone file 'Temp/S.s4p'.
% er.x = 9.1341
% er.y = 9.1346
% er.z = 1.177
% Testing erdes = 11.000000.
% Reading Touchstone file 'Temp/S0.s4p'.
% Reading Touchstone file 'Temp/S.s4p'.
% er.x = 10.3823
% er.y = 10.3786
% er.z = 1.2035
% Testing erdes = 12.000000.
% Reading Touchstone file 'Temp/S0.s4p'.
% Reading Touchstone file 'Temp/S.s4p'.
% er.x = 11.574
% er.y = 11.5702
% er.z = 1.22