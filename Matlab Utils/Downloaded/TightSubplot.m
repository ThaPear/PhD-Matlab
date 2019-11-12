% ---- TightSubPlot.m

function [ha, pos] = TightSubplot(Nx, Ny, gap, MarginX, MarginY)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nx, Ny, [GapX GapY], MarginX, MarginY)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(MarginY); MarginY = .05; end
if nargin<5; MarginX = .05; end

if numel(gap)==1; 
    gap = [gap gap];
end
if numel(MarginX)==1; 
    MarginX = [MarginX MarginX];
end
if numel(MarginY)==1; 
    MarginY = [MarginY MarginY];
end

gap = fliplr(gap);

axh = (1-sum(MarginY)-(Ny-1)*gap(1))/Ny; 
axw = (1-sum(MarginX)-(Nx-1)*gap(2))/Nx;


py = 1-MarginY(2)-axh; 

% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Ny
    px = MarginX(1);
    
    for ix = 1:Nx
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
