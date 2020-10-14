SetupPath;
clear;
close all;

C = inf;
zfeed = 80;

folder = sprintf('%s/ Reports/WP_2100&2400/%s/', resultdir_cst, mfilename);

filepath = [folder, 'final', '.cst'];
touchstonefile = [filepath(1:end-4), '.s4p'];

[parameters, S] = CST.LoadData(touchstonefile);
fs = parameters.frequencies;
fsnew = (9:33)*1e9;

[figS, axS] = a_figS([]);
[figVSWR, axVSWR] = a_figVSWR([]);

S11 = squeeze(S(1,1,:));
S22 = squeeze(S(2,2,:));
S33 = squeeze(S(3,3,:));
S44 = squeeze(S(4,4,:));

S11re = interp1(fs, real(S11.'), fsnew);
S11im = interp1(fs, imag(S11.'), fsnew);
S11 = S11re + 1j .* S11im;

S22re = interp1(fs, real(S22.'), fsnew);
S22im = interp1(fs, imag(S22.'), fsnew);
S22 = S22re + 1j .* S22im;

S33re = interp1(fs, real(S33.'), fsnew);
S33im = interp1(fs, imag(S33.'), fsnew);
S33 = S33re + 1j .* S33im;

S44re = interp1(fs, real(S44.'), fsnew);
S44im = interp1(fs, imag(S44.'), fsnew);
S44 = S44re + 1j .* S44im;

VSWR11 = S2VSWR(S11);
VSWR22 = S2VSWR(S22);
VSWR33 = S2VSWR(S33);
VSWR44 = S2VSWR(S44);

plot(axS, fsnew/1e9, 20*log10(abs(S11)));
plot(axS, fsnew/1e9, 20*log10(abs(S22)));
plot(axS, fsnew/1e9, 20*log10(abs(S33)), '--');
plot(axS, fsnew/1e9, 20*log10(abs(S44)), '--');

plot(axVSWR, fsnew/1e9, VSWR11);
plot(axVSWR, fsnew/1e9, VSWR22);
plot(axVSWR, fsnew/1e9, VSWR33, '--');
plot(axVSWR, fsnew/1e9, VSWR44, '--');

xlim(axVSWR, [12 32]);
ylim(axVSWR, [1 8]);

