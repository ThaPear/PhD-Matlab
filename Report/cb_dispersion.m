SetupPath;
clear;
close all;

c0 = Constants.c0;
z0 = Constants.z0;
f0 = 31e9;          % 13.75 to 31 GHz.

z1 = 80;
z2 = z0;
zfeed = 80;

dx = 4.35e-3;
dy = dx;

p = dx / 2;
gamma = 0.2;
N = 2;
f0match = 21e9;
f0design = 29e9;

linewidth = 1;
axlinewidth = 1;

slab = ChebyshevADS(p, gamma, z1, z2, N, f0match, f0design, 1);

%% Perform simulation.
fs = 1e9 * (0.1:0.1:40);
slab.PlotEpsilons(fs);

hFig = gcf;
hAx = hFig.CurrentAxes;

alignplot(hFig, 8, 4, hFig.Number, [], 2);

ylim(hAx, [-20 10]);
xlim(hAx, [12 32]);

legend(hAx, {'Section 1', 'Section 2'});
ylabel('Error (%)');

hAx.LineWidth = 1;
for(i = 1:length(hAx.Children))
    hAx.Children(i).LineWidth = 1;
end

Zs = Chebyshev.GetImpedances(gamma, z1, z2, N);
ers = (z0 ./ Zs).^2;

plot(hAx, [min(fs) max(fs) nan min(fs) max(fs)]./1e9, [ers(1) ers(1) nan ers(2) ers(2)], 'k');






















