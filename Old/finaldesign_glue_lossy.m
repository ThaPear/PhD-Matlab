close all;
clear;
clear global;

folder = 'h:\My Documents\Thesis\CST\All Simulations\Report\final_glue_lossy\';

f0 = 10e9;
lambda0 = 3e8/f0;

splot = 1;
zreplot     = splot+1;
zimplot     = splot+2;
gplot       = splot+3;
eplot       = splot+4;
eplot2      = splot+5;
vswrplot    = splot+6;

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
    alignplot(gcf, 3, 3, get(gcf, 'Number'), [], 2);

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
    alignplot(gcf, 3, 3, get(gcf, 'Number'), [], 2);

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
    alignplot(gcf, 3, 3, get(gcf, 'Number'), [], 2);
    
%%
if(ishandle(gplot)); close(gplot); end
[fG, G] = ReadCST_ASCIIRealImag([folder, 'realizedgainvsfrequency.txt']);
lambda = 3e8./(fG*1e9);
Nx = 32; Ny = Nx;
dx = 13.5e-3; dy = dx;
Gideal = 4*pi*(Nx*dx*Ny*dy)./lambda.^2;
figureex(gplot);
    plot(fG, 10*log10(Gideal), 'k');
    addlegendentry('Ideal Directivity');
    plot(fG, G, 'k:');
    addlegendentry('Realized Gain');
    xlim([2 10]);
    ylim([20 35]);
    xlabel('Frequency [GHz]');
    ylabel('Directivity [dB] / Realized gain [dB]  ');
    alignplot(gcf, 3, 3, get(gcf, 'Number'), [], 2);
    

%%
% if(ishandle(eplot)); close(eplot); end
% [f, e00] = ReadCST_ASCIIRealImag([folder, 'radefficiency00.txt']);
% [f, e600] = ReadCST_ASCIIRealImag([folder, 'radefficiency600.txt']);
% [f, e6090] = ReadCST_ASCIIRealImag([folder, 'radefficiency6090.txt']);
% 
% figureex(eplot);
%     plot(f, 10*log10(e00), styles{1});
%     addlegendentry('Broadside');
%     plot(f, 10*log10(e600), styles{2});
%     addlegendentry('H-plane');
%     plot(f, 10*log10(e6090), styles{3});
%     addlegendentry('E-plane');
%     alignplot(gcf, 3, 3, get(gcf, 'Number'), [], 2);
%     xlim([2 10]);
%     movelegend('sw');   drawnow;    movelegend('sw');

%%
if(ishandle(eplot2)); close(eplot2); end
[f, e00] = ReadCST_ASCIIRealImag([folder, 'totefficiency00.txt']);
[f, e600] = ReadCST_ASCIIRealImag([folder, 'totefficiency600.txt']);
[f, e6090] = ReadCST_ASCIIRealImag([folder, 'totefficiency6090.txt']);

figureex(eplot2);
    plot(f, 10*log10(e00), styles{1});
    addlegendentry('Broadside');
    plot(f, 10*log10(e600), styles{2});
    addlegendentry('H-plane');
    plot(f, 10*log10(e6090), styles{3});
    addlegendentry('E-plane');
%     plot(fG, G-10*log10(Gideal), 'r');
%     addlegendentry('G / D_{ideal}');
    alignplot(gcf, 3, 3, get(gcf, 'Number'), [], 2);
    xlim([2 10]);
    ylim([-1.5 0]);
    xlabel('Frequency [GHz]');
    ylabel('Efficiency [dB]');
    movelegend('sw');   drawnow;    movelegend('sw');
    
%%
[f, S00] = ReadCST_ASCIIRealImag([folder, 's00.txt']);
[f, S600] = ReadCST_ASCIIRealImag([folder, 's600.txt']);
[f, S6090] = ReadCST_ASCIIRealImag([folder, 's6090.txt']);

VSWR00 = S2VSWR(S00);
VSWR600 = S2VSWR(S600);
VSWR6090 = S2VSWR(S6090);

if(ishandle(vswrplot)); close(vswrplot); end
figureex(vswrplot);
    alignplot(figure(vswrplot), 3, 3, vswrplot, [], 2);
    plot(f, VSWR00, styles{1}, 'LineWidth', 1.5);
    addlegendentry('Broadside');
    plot(f, VSWR600, styles{2}, 'LineWidth', 1.5);
    addlegendentry('H-plane');
    plot(f, VSWR6090, styles{3}, 'LineWidth', 1.5);
    addlegendentry('E-plane');
    xlim([1 12]);
    ylim([1 6]);
    xlabel('Frequency [GHz]');
    ylabel('VSWR');

    
%%
for(i = 1:7)
    alignplot(figure(i), 3, 3, i, [], 2);
end