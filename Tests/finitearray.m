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
Nx = 3;
Ny = 3;
fs = (9:4:35) * 1e9;
% fs2 = (10:4:35) * 1e9;
% fs3 = (11:4:35) * 1e9;
% fs4 = (12:4:35) * 1e9;
% fs = [fs fs2 fs3 fs4];

[~, hostname] = system('hostname'); hostname = strsplit(hostname, '\n');
if(strcmp(hostname{1}, 'SRV539'))
    Nx = 1;
    Ny = 1;
    fs = (12:4:52)*1e9;
%     fs = [fs, (33:35)*1e9];
end

excitation = ones(Nx, Ny);

% slab = StagedADS(dx/2, l0/4, 3, 1, f0);
% slab = slab.elements{1};
% slab = TLine({slab, slab, slab});
% tlineup = TerminatedTLine(slab, FreeSpace());
tlineup = FreeSpace();
tlinedown = FreeSpace();
% tlinedown = ShortedLine(1, 0.25*l0);

slot = Slot(dx, dy, wslot, dslot, walled);


th = eps * pi/180;
ph = 0 * pi/180;

if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
    array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
%     array = FiniteArray_FixedDs(slot, Nx, Ny, dedge, zfeed);
end
% Do the frequency points one-by-one.
% for(f = fs)
%     array.InitializeDs(f);
%     save(sprintf('%ix%i_Ds', Nx, Ny));
%     array.InitializeZMatrix(f);
%     save(sprintf('%ix%i_ZMatrix', Nx, Ny));
%     Zas = array.GetInputImpedance(f, excitation);
% end

array.InitializeDs(fs);
array.InitializeZMatrix(fs);
Zas = array.GetInputImpedance(fs, excitation);

%%
if(length(fs) > 1)
    [hFig, hAx] = figureex;
        hAx.ColorOrder = lines(Nx*Ny);%reshape(repmat(lines(7), 1, 2).', [], 14).';
        plot(hAx, fs./1e9, real(reshape(Zas, [], length(fs))));
        addlegendentry(hAx, 'Real');
        plot(hAx, fs./1e9, imag(reshape(Zas, [], length(fs))), '--');
        addlegendentry(hAx, 'Imag');
        title(hAx,  'Zas');
else
    [hFig, hAx] = figureex;
        hAx.ColorOrder = lines(Nx*Ny);%reshape(repmat(lines(7), 1, 2).', [], 14).';
        plot(hAx, fs./1e9, real(reshape(Zas, [], length(fs))), 'o');
        addlegendentry(hAx, 'Real');
        plot(hAx, fs./1e9, imag(reshape(Zas, [], length(fs))), 'x');
        addlegendentry(hAx, 'Imag');
        title(hAx,  'Zas');
end
    
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

    