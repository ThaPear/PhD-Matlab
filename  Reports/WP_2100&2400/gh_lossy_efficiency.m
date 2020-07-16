clear;
SetupPath;
close all;

folder = 'e:\ Reports\WP_2100&2400\gg_coax_lossy\';

% [parameters, tot2] = CST.LoadData([folder, 'tot2.txt']);
% tot2 = 10.^(tot2{1}(2,:)./10);
% 
% [parameters, tot4] = CST.LoadData([folder, 'tot4.txt']);
% tot4 = 10.^(tot4{1}(2,:)./10);

[parameters, rad00] = CST.LoadData([folder, 'rad00.txt']);
fs = rad00{1}(1,:);
rad00 = 10.^(rad00{1}(2,:)./10);

[parameters, rad600] = CST.LoadData([folder, 'rad600.txt']);
rad600 = 10.^(rad600{1}(2,:)./10);

[parameters, rad6090] = CST.LoadData([folder, 'rad6090.txt']);
rad6090 = 10.^(rad6090{1}(2,:)./10);

[parameters, S] = CST.LoadData([folder, 'S80_1.s2p']);
S00 = interp1(parameters.frequencies./1e9, squeeze(S(1,1,:)), fs);
[parameters, Sscan] = CST.LoadData([folder, 'S80_2.s2p']);
S600 = interp1(parameters.frequencies./1e9, squeeze(Sscan(1,1,:)), fs);
S6090 = interp1(parameters.frequencies./1e9, squeeze(Sscan(2,2,:)), fs);

mtot00 = rad00 .* (1 - abs(S00).^2);
mtot600 = rad600 .* (1 - abs(S600).^2);
mtot6090 = rad6090 .* (1 - abs(S6090).^2);

[figRad, axRad] = a_fig;
%     plot(axRad, fs, 10*log10(rad00));
%     plot(axRad, fs, 10*log10(rad600));
%     plot(axRad, fs, 10*log10(rad6090));
%     ylim(axRad, [-3 0]);
    plot(axRad, fs, 100*(rad00));
    plot(axRad, fs, 100*(rad600));
    plot(axRad, fs, 100*(rad6090));
    ylim(axRad, [0 100]);
    xlabel(axRad, 'Frequency [GHz]');
    ylabel(axRad, 'Efficiency [dB]');
    legend(axRad, axRad.Children(3:-1:1), {'Broadside', 'H-plane 60\circ', 'E-plane 60\circ'});
    drawnow; movelegend(axRad, 's');
    
[figTot, axTot] = a_fig;
    plot(axTot, fs, 10*log10(mtot00));
    plot(axTot, fs, 10*log10(mtot600));
    plot(axTot, fs, 10*log10(mtot6090));
    ylim(axTot, [-3 0]);
%     plot(axTot, fs, mtot00*100);
%     plot(axTot, fs, mtot600*100);
%     plot(axTot, fs, mtot6090*100);
%     ylim(axTot, [0 100]);
    xlabel(axTot, 'Frequency [GHz]');
    ylabel(axTot, 'Efficiency [dB]');
    legend(axTot, axTot.Children(3:-1:1), {'Broadside', 'H-plane 60\circ', 'E-plane 60\circ'});
    drawnow; movelegend(axTot, 's');

% [figS, axS] = a_figS;
%     plot(axS, fs, 20*log10(abs(S00)));
%     plot(axS, fs, 20*log10(abs(S600)));
%     plot(axS, fs, 20*log10(abs(S6090)));
    
[figVSWR, axVSWR] = a_figVSWR;
    plot(axVSWR, parameters.frequencies./1e9, S2VSWR(squeeze(S(1,1,:))));
    plot(axVSWR, parameters.frequencies./1e9, S2VSWR(squeeze(Sscan(1,1,:))));
    plot(axVSWR, parameters.frequencies./1e9, S2VSWR(squeeze(Sscan(2,2,:))));
%     plot(axVSWR, fs, S2VSWR(S600));
%     plot(axVSWR, fs, S2VSWR(S6090));
    legend(axVSWR, axVSWR.Children(3:-1:1), {'Broadside', 'H-plane 60\circ', 'E-plane 60\circ'});
    drawnow; movelegend(axVSWR, 'n');

