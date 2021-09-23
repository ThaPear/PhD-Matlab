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

% slab = StagedADS(dx/2, l0/4, 3, 1, f0);
% slab = slab.elements{1};
% slab = TLine({slab, slab, slab});

% Number of unit cells.
Nx = 1;
Ny = 3;
% Nx = 3;
% Ny = 3;
excitation = ones(Nx, Ny);
% tlineup = TerminatedTLine(slab, FreeSpace());
tlineup = FreeSpace();
% tlinedown = FreeSpace();
tlinedown = ShortedLine(1, 0.25*l0);
% tlineup = TerminatedTLine(StagedADS(dx/2, l0/4, 3, 1, f0), FreeSpace());

slot = Slot(dx, dy, wslot, dslot, walled);

fs = (10:2.5:35) * 1e9;
th = eps * pi/180;
ph = 0 * pi/180;

tc = tic;

% if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
%     array = FiniteArray_Slow(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
% end
% array.InitializeDs(fs);
% save(sprintf('%ix%i_Ds', Nx, Ny));
% % array.InitializeKyInts(fs);
% % save(sprintf('%ix%i_Ds', Nx, Ny));
% array.InitializeZMatrix(fs);
% save(sprintf('%ix%i_Ds', Nx, Ny));
% Zas = array.GetInputImpedance(fs, excitation);




% if(~exist('array', 'var') || ~isa(array, 'FiniteArray') || array.Nx ~= Nx || array.Ny ~= Ny)
    array = FiniteArray(slot, tlineup, tlinedown, Nx, Ny, dedge, zfeed);
% end

array.InitializeDs(fs);
array.InitializeZMatrix(fs);
% array.InitializeZMatrix_fs(fs);
Zas = array.GetInputImpedance(fs, excitation);
dt = toc(tc);
dispex('Total computation took %.1fs.\n', dt);

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

hFig = figureex;
for(nx = 1:Nx)
    for(ny = 1:Ny)
        hAx = subplot(Nx, Ny, ny+(nx-1)*Ny);
        hold(hAx, 'on');
        grid(hAx, 'on');
        box(hAx, 'on');
        plot(hAx, fs/1e9, squeeze(real(Zas(nx, ny, :))), 'r');
        plot(hAx, fs/1e9, squeeze(imag(Zas(nx, ny, :))), 'r--');
    end
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

%{*
% global kxs;
% global vs;

% f = fs(1);
% 
% for(deformedpath = 0:1)
% %     close(100+deformedpath);
%     hFig = figureex(2+deformedpath);
%     for(dny = 0:Ny-1)
%         hAx = subplot(Ny, 1, dny+1);
%         kx = kxs{dny+1, deformedpath+1};
%         v = vs{dny+1, deformedpath+1};
%         if(isempty(kx))
%             continue;
%         end
%         [kx, I] = sort(kx, 'ComparisonMethod', 'real');
%         v = v(I);
% 
%         if(deformedpath)
%             % Go to -inf * 1j
%             k0 = 2*pi*f/c0;
%             delta = 0.01*k0;
%             lim1 = -5j.*k0-1.*delta;
%             lim2 = -5j.*k0+1.*delta+1.5*k0;
%             % Deform integration path around branch cuts.
%             integrationpath = [(-1-1j).*delta, (1+1j).*delta, 1.5*k0+1j.*delta, 1.5*k0+1.*delta];
% 
%             rightTail = logical((real(kx) > 0) & (imag(kx) <= imag(integrationpath(end))));
%             v(rightTail) = fliplr(v(rightTail));
%             kx(rightTail) = fliplr(kx(rightTail));
% 
%             kx = FiniteArray.UnfoldKVector(kx, integrationpath);
%         end
% 
%         k0 = 2*pi*f/c0;
%         plot(hAx, real(kx)./k0, real(v), real(kx)./k0, imag(v));
%         hold(hAx, 'on');
%         % figureex; plot(real(v)); plot(imag(v));
%         title(hAx, sprintf('%s, ny = %i, deformedpath = %i', class(array), dny, deformedpath));
%     end
% end

%}