SetupPath;
clear;
close all;

z0 = Constants.z0;
z1 = 80;
z2 = 80;

f0 = 22e9;
lam0 = 3e8/f0;
fs = (10:0.01:35) .* 1e9;

Zs = Chebyshev.GetImpedances(0.05,80,80*2^3,5)
% Zs = Binomial.GetImpedances(80,80*2^3,5)
Zs = fliplr(Zs) ./ [8 4 4 2 1]
ers = (z0 ./ Zs).^2;

if(1)
    line1 = Line(lam0/4 / sqrt(ers(1)), ers(1));
    line2 = Line(lam0/4 / sqrt(ers(2)), ers(2));
    line3 = Line(lam0/4 / sqrt(ers(3)), ers(3));
    line4 = Line(lam0/4 / sqrt(ers(4)), ers(4));
    line5 = Line(lam0/4 / sqrt(ers(5)), ers(5));
    line5t = TerminatedTLine(line5, Impedance(80));

    lines = {line1, Shunt(TLine({line2, line3, Shunt(TLine({line4, Shunt(TLine({line5t})), ...
                                                                                line5t})), ...
                                                            line4, Shunt(TLine({line5t})), ...
                                                                                line5t})), ...
                                 line2, line3, Shunt(TLine({line4, Shunt(TLine({line5t})), ...
                                                                                line5t})), ...
                                                            line4, Shunt(TLine({line5t})), ...
                                                                                line5};
    line = TLine(lines);
else
    lines = {};
    for(i = 1:length(Zs))
        lines = [lines, {Line(lam0/4 / sqrt(ers(i)), ers(i))}];
    end
    line = TLine(lines);
end

th = eps * pi/180;
ph = 0 * pi/180;

% Calculate propagation constants.
[k0, kx0, ky0, kz0] = k(fs, 1, th, ph);
% [kd, ~, ~, kzd] = k(fs, (z0/z1)^2, th, ph);
kr = sqrt(kx0.^2 + ky0.^2);
% Calculate impedances.
[z0, z0te, z0tm] = z(1, k0, kz0);

% Calculate ABCD matrix - TE only since only H-plane is important.
isTE = 0;
ABCDte = line.GetABCD(isTE, fs, k0, kr);
Ste = ABCD2S(ABCDte, z1, z2);

[figS, axS] = figureex(1);
    plot(axS, fs/1e9, 20*log10(abs(Ste.s11)));
    ylim([-30 0]);

termline = TerminatedTLine(line, Impedance(80));
Zin = termline.GetInputImpedance(isTE, fs, k0, kr);
[figZ, axZ] = figureex(2);
    plot(axZ, fs/1e9, real(Zin), fs/1e9, imag(Zin));