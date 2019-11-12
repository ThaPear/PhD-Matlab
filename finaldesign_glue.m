close all;
clear;
clear global;

folder = 'h:\My Documents\Thesis\CST\All Simulations\Report\final_glue\';

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
[f, Xx] = ReadCST_ASCIIRealImag([folder, 'xpolx.txt']);
[f, Xy] = ReadCST_ASCIIRealImag([folder, 'xpoly.txt']);

if(ishandle(xplot)); close(xplot); end
figureex(xplot);
    plot(f, Xx, 'k.-');
    addlegendentry('x-oriented slot');
    plot(f, Xy, 'k.--');
    addlegendentry('y-oriented slot');
    
    xlim([2 10]);
    xlabel('Frequency [GHz]');
    ylabel('Crosspol [dB]');
    
[f, Xcan6] = ReadCST_ASCIIRealImag([folder, 'xpol6.txt']);
[f, Xcan7] = ReadCST_ASCIIRealImag([folder, 'xpol7.txt']);
[f, Xcan8] = ReadCST_ASCIIRealImag([folder, 'xpol8.txt']);

if(ishandle(xcanplot)); close(xcanplot); end
figureex(xcanplot);
    plot(f, Xcan6, 'k.-');
    addlegendentry('Cancelled at 6 GHz');
    plot(f, Xcan7, 'k.--');
    addlegendentry('Cancelled at 7 GHz');
    plot(f, Xcan8, 'k.:');
    addlegendentry('Cancelled at 8 GHz');
    
    xlim([2 10]);
    xlabel('Frequency [GHz]');
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

%%
for(i = 1:100)
    if(~ishandle(i))
        break;
    end
    alignplot(figure(i), 3, 3, i, [], 2);
end