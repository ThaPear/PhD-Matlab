% close all;
closewaitbars;
clearvars -except array;
SetupPath;
clear global;

f0 = 31e9;
c0 = Constants.c0;
l0 = c0/f0;

dx = 0.45*l0;
dy = dx;
wslot = 0.05*l0;
dslot = 0.05*l0;
walled = 0;
dedge = 0.25*l0;
zfeed = 100;

% Number of unit cells.
Nx = 1;
Ny = 1;
fs = (26:2:40)*1e9;

[~, hostname] = system('hostname'); hostname = strsplit(hostname, '\n');
if(strcmp(hostname{1}, 'SRV539'))
    Nx = 1;
    Ny = 1;
    fs = (12:4:52)*1e9;
%     fs = [fs, (33:35)*1e9];
end

excitation = ones(Nx, Ny);

tlineup = FreeSpace();
tlinedown = FreeSpace();

slot = Slot(dx, dy, wslot, dslot, walled);


th = eps * pi/180;
ph = 0 * pi/180;

if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
    array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
end
array.InitializeZMatrix(fs);
% tic;
% array = FiniteArray_NoPrecompute(slot, Nx, Ny, dedge, zfeed);
% array.InitializeZMatrix(fs);
% toc

%%
Zas = array.GetInputImpedance(fs, excitation);
[hFig, hAx] = figureex;
    hAx.ColorOrder = lines(Nx*Ny);%reshape(repmat(lines(7), 1, 2).', [], 14).';
    plot(hAx, fs./1e9, real(reshape(Zas, [], length(fs))));
    addlegendentry(hAx, 'Real');
    plot(hAx, fs./1e9, imag(reshape(Zas, [], length(fs))), '--');
    addlegendentry(hAx, 'Imag');
    title(hAx,  'Zas');
    
%%
% Zmat = array.Zmat;
% 
% [hFig, hAx] = figureex;
%     hAx.ColorOrder = lines(((Nx+2)*Ny)^2);%reshape(repmat(lines(7), 1, 2).', [], 14).';
%     plot(hAx, fs./1e9, real(reshape(Zmat, [], length(fs))));
%     addlegendentry(hAx, 'Real');
%     plot(hAx, fs./1e9, imag(reshape(Zmat, [], length(fs))), '--');
%     addlegendentry(hAx, 'Imag');
%     title(hAx,  'Zmat');

    