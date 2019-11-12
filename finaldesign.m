close all;
clear;
clear global;

folder = 'h:\My Documents\Thesis\CST\All Simulations\Report\final\';

f0 = 10e9;
lambda0 = 3e8/f0;

splot = 1;
zreplot     = splot+1;
zimplot     = splot+2;
gplot       = splot+3;
pat0plot    = splot+4;
pat90plot   = splot+5;
xplot       = splot+6;
xcanplot    = splot+7;
vswrplot    = splot+8;

styles = {'k', 'k--', 'k:'};

%%

[f, S00] = ReadCST_ASCIIRealImag([folder, 's00.txt']);
[f, S600] = ReadCST_ASCIIRealImag([folder, 's600.txt']);
[f, S6090] = ReadCST_ASCIIRealImag([folder, 's6090.txt']);

figureex(splot);
    plot(f, 20*log10(abs(S00)), styles{1});
    addlegendentry('Broadside');
    plot(f, 20*log10(abs(S600)), styles{2});
    addlegendentry('H-plane');
    plot(f, 20*log10(abs(S6090)), styles{3});
    addlegendentry('E-plane');
    xlim([1 12]);
    ylim([-30 0]);
    xlabel('Frequency [GHz]');
    ylabel('Active |\Gamma|');

[f, Z00] = ReadCST_ASCIIRealImag([folder, 'z00.txt']);
[f, Z600] = ReadCST_ASCIIRealImag([folder, 'z600.txt']);
[f, Z6090] = ReadCST_ASCIIRealImag([folder, 'z6090.txt']);

figureex(zreplot);
    plot(f, real(Z00), styles{1});
    addlegendentry('Broadside');
    plot(f, real(Z600), styles{2});
    addlegendentry('H-plane');
    plot(f, real(Z6090), styles{3});
    addlegendentry('E-plane');
    xlim([1 12]);
    ylim([0 200]);
    xlabel('Frequency [GHz]');
    ylabel('Active Input Resistance [\Omega]');

figureex(zimplot);
    plot(f, imag(Z00), styles{1});
    addlegendentry('Broadside');
    plot(f, imag(Z600), styles{2});
    addlegendentry('H-plane');
    plot(f, imag(Z6090), styles{3});
    addlegendentry('E-plane');
    xlim([1 12]);
    ylim([-150 100]);
    xlabel('Frequency [GHz]');
    ylabel('Active Input Reactance [\Omega]');
    
%%
[f, G] = ReadCST_ASCIIRealImag([folder, 'realizedgainvsfrequency.txt']);
lambda = 3e8./(f*1e9);
Nx = 32; Ny = Nx;
dx = 13.5e-3; dy = dx;
Gideal = 4*pi*(Nx*dx*Ny*dy)./lambda.^2;
figureex(gplot);
    plot(f, 10*log10(Gideal), 'k');
    plot(f, G, 'k:');
    
%%
f0 = 10e9;
lambda0 = 3e8./(f0);

% Phi = 0
[th, G00] = ReadCST_ASCIIRealImag([folder, 'pattern00.txt']);
G00 = imag(G00);
th = th(G00 > -200);
G00 = G00(G00 > -200);

[th, G600] = ReadCST_ASCIIRealImag([folder, 'pattern600.txt']);
G600 = imag(G600);
th = th(G600 > -200);
G600 = G600(G600 > -200);

[th, G300] = ReadCST_ASCIIRealImag([folder, 'pattern300.txt']);
G300 = imag(G300);
th = th(G300 > -200);
G300 = G300(G300 > -200);

Nx = 32; Ny = Nx;
dx = 13.5e-3; dy = dx;
Gideal = 4*pi*(Nx*dx*Ny*dy*cosd(th))./lambda0.^2;

figureex(pat0plot);
    plot(th, G00, 'k');
    addlegendentry('\theta_\text{scan}=\SI{0}{\degree}');
    plot(th, G300, 'k--');
    addlegendentry('\theta_\text{scan}=\SI{30}{\degree}');
    plot(th, G600, 'k:');
    addlegendentry('\theta_\text{scan}=\SI{60}{\degree}');
    plot(th, 10*log10(Gideal), 'k');
    addlegendentry('4\piA/\lambda^2');
    
    xlim([-90 90]);
    xlabel('\theta [\SI{\degree}]');
    ylim([5 35]);
    ylabel('Realized Gain [dB]');
    
% Phi = 90
[th, G090] = ReadCST_ASCIIRealImag([folder, 'pattern090.txt']);
G090 = imag(G090);
th = th(G090 > -200);
G090 = G090(G090 > -200);
    
[th, G3090] = ReadCST_ASCIIRealImag([folder, 'pattern3090.txt']);
G3090 = imag(G3090);
th = th(G3090 > -200);
G3090 = G3090(G3090 > -200);
    
[th, G6090] = ReadCST_ASCIIRealImag([folder, 'pattern6090.txt']);
G6090 = imag(G6090);
th = th(G6090 > -200);
G6090 = G6090(G6090 > -200);

figureex(pat90plot);
    plot(th, G090, 'k');
    addlegendentry('\theta_\text{scan}=\SI{0}{\degree}');
    plot(th, G3090, 'k--');
    addlegendentry('\theta_\text{scan}=\SI{30}{\degree}');
    plot(th, G6090, 'k:');
    addlegendentry('\theta_\text{scan}=\SI{60}{\degree}');
    plot(th, 10*log10(Gideal), 'k');
    addlegendentry('4\piA/\lambda^2');
    
    xlim([-90 90]);
    xlabel('\theta [\SI{\degree}]');
    ylim([5 35]);
    ylabel('Realized Gain [dB]');
    
%%
[f, Xx] = ReadCST_ASCIIRealImag([folder, 'crosspolx.txt']);
[f, Xy] = ReadCST_ASCIIRealImag([folder, 'crosspoly.txt']);

figureex(xplot);
    plot(f, Xx, 'k.-');
    plot(f, Xy, 'k.--');
    
    xlim([2 10]);
    xlabel('Frequency');
    ylabel('Crosspol [dB]');
    
[f, Xcan] = ReadCST_ASCIIRealImag([folder, 'crosspolcancelled.txt']);
[f, Xcan5] = ReadCST_ASCIIRealImag([folder, 'crosspolcancelled5.txt']);
[f, Xcan6] = ReadCST_ASCIIRealImag([folder, 'crosspolcancelled6.txt']);

figureex(xcanplot);
    plot(f, Xcan, 'k.-');
    plot(f, Xcan5, 'k.--');
    plot(f, Xcan6, 'k.:');
    
    xlim([2 10]);
    xlabel('Frequency');
    ylabel('Crosspol [dB]');
    
%%
[f, S00] = ReadCST_ASCIIRealImag([folder, 's00.txt']);
[f, S600] = ReadCST_ASCIIRealImag([folder, 's600.txt']);
[f, S6090] = ReadCST_ASCIIRealImag([folder, 's6090.txt']);

VSWR00 = S2VSWR(S00);
VSWR600 = S2VSWR(S600);
VSWR6090 = S2VSWR(S6090);

if(ishandle(vswrplot)); close(vswrplot); end
figureex(vswrplot);
    plot(f, VSWR00, styles{1});
    addlegendentry('Broadside');
    plot(f, VSWR600, styles{2});
    addlegendentry('H-plane');
    plot(f, VSWR6090, styles{3});
    addlegendentry('E-plane');
    xlim([1 12]);
    ylim([1 6]);
    xlabel('Frequency [GHz]');
    ylabel('VSWR');