
f0 = 22e9;
p = 4.4e-3/2;
L = 1.7e-3;
erdes = 8;
ads = DesignADS(f0, p, L, erdes);

CST.GetEpsilonMu(ads, f0);