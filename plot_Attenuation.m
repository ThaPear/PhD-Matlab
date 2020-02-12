distance = 4.35e-3 * (16 + 16);

tand = 0.0009;
er = 2.2;
z0 = 80;
b = (0.1:0.001:1)*1e-3;
W = Stripline.GetWoverb(z0, er) .* b;
Rs = 0.026;
t = 0.01e-3;
f = 22e9;
alpha = Stripline.Attenuation(er, W, b, t, Rs, f, tand);
figureex; plot(b, 20*log10(exp(alpha))*distance);

d = b;
w = Microstrip.GetWoverD(z0, er) .* d;
% w = w(1);
alpha = Microstrip.Attenuation(er, w, d, f, tand);
% yyaxis right;
plot(d, 20*log10(exp(alpha))*distance);

legend({'Stripline', 'Microstrip'});
xlabel('b or d [m]');
ylabel('\alpha [dB]');

ee = Microstrip.EpsilonEffective(er, w, d);
ff = (er .* (ee-1)) ./ (ee .* (er-1));