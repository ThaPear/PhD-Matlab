clear;
SetupPath;
close all;

%% Test struct2string
% strct = [];
% strct.a = 5;
% strct.b = 'hiya';
% strct.c = [];
% strct.c.hi = 10.123456789;
% strct.c.hi2 = 'foo';
% strct.aa = 20;
% str = struct2string(strct);


%% Cosine with varying frequency.
% x = 0.1:0.0001:10*pi;
% f = linspace(480,750,length(x))./100;
% % f = linspace(20,0.5,length(x));
% v(1) = 0;
% dx = x(2) - x(1);
% for(i = 2:length(x))
%     v(i) = v(i-1)+f(i)*dx;
% end
% 
% % x = x(1:find(cos(v) == 0, 1, 'last'));
% % v = v(1:find(cos(v) == 0, 1, 'last'));
% figure; plot(x, cos(v))

%% Reflection of 5-section Chebyshev through splitters.
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
    line1 = Line(ers(1), lam0/4 / sqrt(ers(1)));
    line2 = Line(ers(2), lam0/4 / sqrt(ers(2)));
    line3 = Line(ers(3), lam0/4 / sqrt(ers(3)));
    line4 = Line(ers(4), lam0/4 / sqrt(ers(4)));
    line5 = Line(ers(5), lam0/4 / sqrt(ers(5)));
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
        lines = [lines, {Line(ers(i), lam0/4 / sqrt(ers(i)))}];
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


%% Input impedance at arbitratry point in Chebyshev
% z0 = Constants.z0;
% z1 = 80;
% z2 = 80;
% 
% Zs = ChebyshevADS.GetImpedances(0.05,80,80*2^3,5);
% % Zs = fliplr(Zs);
% 
% Zin1 = Zs(5)^2 / z2
% Zin2 = Zs(4)^2 / (Zin1)
% Zin3 = Zs(3)^2 / (Zin2)
% Zin4 = Zs(2)^2 / (Zin3)
% Zin5 = Zs(1)^2 / (Zin4)
% 
% z0 = Constants.z0;
% 
% Zs = Zs ./ [1 2 4 4 8];
% ers = (z0 ./ Zs).^2;
% 
% f0 = 22e9;
% l0 = 3e8/f0;
% fs = (10:0.1:32)*1e9;
% ls = 3e8./fs;
% betas = 2*pi./ls;
% 
% l = l0/4;
% l = l ./ ones(1,5);%sqrt(ers);
% 
% z1 = 80;
% z2 = 80;
% 
% Zin1 = Zs(5)^2 / z2
% Zin2 = Zs(4)^2 / (Zin1/2)
% Zin3 = Zs(3)^2 / (Zin2/2)
% Zin4 = Zs(2)^2 / (Zin3)
% Zin5 = Zs(1)^2 / (Zin4/2)
% 
% Zin1 = Zs(5) .* (z2     + 1j .* Zs(5) .* tan(betas .* l(5))) ./ (Zs(5) + 1j .* z2     .* tan(betas .* l(5)));
% Zin2 = Zs(4) .* (Zin1/2 + 1j .* Zs(4) .* tan(betas .* l(4))) ./ (Zs(4) + 1j .* Zin1/2 .* tan(betas .* l(4)));
% Zin3 = Zs(3) .* (Zin2/2 + 1j .* Zs(3) .* tan(betas .* l(3))) ./ (Zs(3) + 1j .* Zin2/2 .* tan(betas .* l(3)));
% Zin4 = Zs(2) .* (Zin3   + 1j .* Zs(2) .* tan(betas .* l(2))) ./ (Zs(2) + 1j .* Zin3   .* tan(betas .* l(2)));
% Zin5 = Zs(1) .* (Zin4/2 + 1j .* Zs(1) .* tan(betas .* l(1))) ./ (Zs(1) + 1j .* Zin4/2 .* tan(betas .* l(1)));
% 
% Zin5(fs == f0)
% 
% gammas = (Zin5 - z1) ./ (Zin5 + z1);
% 
% figureex;
%     plot(fs/1e9, real(Zin5), fs/1e9, imag(Zin5));
% figureex;
%     plot(fs/1e9, 20*log10(abs(gammas)));
%     ylim([-30 0]);



% Zin1 = Z1 .* (Z0 + 1j .* Z1 .* tan(beta .* l1)) ./ (Z1 + 1j .* Z0 .* tan(beta .* l1));

%% Binomial test
% for(N = 1:5)
% z0 = Constants.z0;
% z1 = 50;
% z2 = 100;
% 
% f0 = 22e9;
% l0 = 3e8/f0;
% fs = (10:0.1:35) .* 1e9;
% 
% Zs = Binomial.GetImpedances(z1, z2, N);
% 
% lines = {};
% for(i = 1:N)
%     er = (z0 ./ Zs(i)).^2;
%     lines = [lines, {Line(er, l0/4 / sqrt(er))}];
% end
% line = TLine(lines);
% 
% th = eps * pi/180;
% ph = 0 * pi/180;
% 
% % Calculate propagation constants.
% [k0, kx0, ky0, kz0] = k(fs, 1, th, ph);
% [kd, ~, ~, kzd] = k(fs, (z0/z1)^2, th, ph);
% kr = sqrt(kx0.^2 + ky0.^2);
% % Calculate impedances.
% [z0, z0te, z0tm] = z(1, k0, kz0);
% 
% % Calculate ABCD matrix - TE only since only H-plane is important.
% isTE = 1;
% ABCDte = line.GetABCD(isTE, fs, k0, kr);
% Ste = ABCD2S(ABCDte, z1, z2);
% 
% [figS, axS] = figureex(1);
% %     plot(axS, fs/1e9, 20*log10(abs(Ste.s11)));
% %     ylim([-30 0]);
%     plot(axS, fs/1e9, (abs(Ste.s11)));
% end

%% Server hyperthreading benchmark
% delete(gcp('nocreate'));
% ppool = parpool('local', 32);
% tci = tic;
% for(iii = 1:50)
%     clearvars -except iii tci ppool;
% %     close all;
%     PreliminaryDesign;
% end
% dt = toc(tci);
% dispex('Using %i workers it took %f seconds.\n', ppool.NumWorkers, dt);
% delete(ppool);
% 
% ppool = parpool('local', 64);
% tci = tic;
% for(iii = 1:50)
%     clearvars -except iii tci ppool;
% %     close all;
%     PreliminaryDesign;
% end
% dt = toc(tci);
% dispex('Using %i workers it took %f seconds.\n', ppool.NumWorkers, dt);
% delete(ppool);

%% HFSS double 45-degree arc center piece length
% r = 1;
% l1 = 0.1;
% rounded = floor(10*(sqrt(2)-2)*r)/10
% ideal = (sqrt(2)-2)*r
% worst = (sqrt(2)-2)*r-1/2*sqrt(2)*2*r
% l = 2*r-(-worst+rounded)*sqrt(2);
% l = 2*r-(-((sqrt(2)-2)*r-1/2*sqrt(2)*2*r)+floor(10*(sqrt(2)-2)*r)/10)*sqrt(2);
% x = (sqrt(2)-2)*r-1/2*sqrt(2)*l;
% 
% 
% l3 = 2.5-l1-2*pi/4*r-l
% l3 = 2.3-(l1+sqrt(2)*r+sqrt(2)/2*l)
% y = -l1-sqrt(2)*r-sqrt(2)/2*l-l3;
% 
% fprintf('l = %.10f\n', l);
% fprintf('x = %.10f\n', x);
% fprintf('y = %.10f\n', y);

%% 5-section Chebyshev
% Zs = ChebyshevADS.GetImpedances(0.05,80,80*2^3,5)
% Zs = Zs./[1 2 4 4 8]
% ee = Microstrip.EpsilonEffective(2.2, Zs)
% qwavelength = 3e8./22e9./sqrt(ee)./4 * 1e3

%% Double Chebyshev with power combiner in between.
% fs = (12:0.01:32) * 1e9;
% lambdas = Constants.c0 ./ fs;
% beta = 2 .* pi ./ lambdas;
% 
% f0 = 21.8e9;
% f1 = f0;%14.125e9;
% f2 = f0;%29.5e9;
% lambda0 = Constants.c0 ./ f0;
% lambda1 = Constants.c0 ./ f1;
% lambda2 = Constants.c0 ./ f2;
% 
% l1 = lambda1 ./ 4;
% l2 = lambda2 ./ 4;
% 
% [Zs] = ChebyshevADS.GetImpedances(0.1,80,160,2);
% Z1 = Zs(1);
% Z2 = Zs(2)/2;
% 
% Zint = 113;
% Z0 = 80;
% % Z1 = (Z0 * Zint) ^ (1/2);
% % for(Z1 = 96:2:100)
% % Z2 = (Z0 * Zint/2) ^ (1/2);
% % for(Z2 = 60:65)
% 
% Zin1 = Z1 .* (Z0 + 1j .* Z1 .* tan(beta .* l1)) ./ (Z1 + 1j .* Z0 .* tan(beta .* l1));
% Zin2 = Zin1 / 2;
% Zin3 = Z2 .* (Zin2 + 1j .* Z2 .* tan(beta .* l2)) ./ (Z2 + 1j .* Zin2 .* tan(beta .* l2));
% 
% Gamma = (Zin3 - Z0) ./ (Zin3 + Z0);
% 
% [fig, ax] = figureex(1);
%     plot(ax, fs./1e9, 20*log10(abs(Gamma)));
%     addlegendentry(sprintf('Z_1: %.1f, Z_2: %.1f', Z1, Z2));
%     ylim(ax, [-30 0]);
% % end
% % end

%%

% Ah1 = 1.9699;
% Ah2 = 5.1679;
% phh1 = 239.87;
% phh2 = 144.44;
% 
% s1 = Ah1 * exp(phh1 * 1j * pi/180);
% s2 = Ah2 * exp(phh2 * 1j * pi/180);
% 
% ay = -s1 / s2
% abs(ay)
% angle(ay) * 180/pi


% tc = tic;
% for(i = 1:100000)
%     val = rand();
%     fprintf('val(%i) = %g\n', i, val);
% end
% dt = toc(tc);
% 
% tc = tic;
% for(i = 1:100000)
%     val = rand();
%     dispex('val(%i) = %g\n', i, val);
% end
% dt2 = toc(tc);
% clc;
% disp(['fprintf = ', num2str(dt), 's\n']);
% disp(['dispex = ', num2str(dt2), 's\n']);

% z0=50:200;
% Z11_00=79.427756+1j*5.5073314;
% 
% gamma = @(z0, Z) abs((Z-z0)./(Z+z0));
% 
% Z11_31=23.550842+1j*-18.277396;
% Z22_31=131.81712+1j*24.106141;
% G11_31 = gamma(z0, Z11_31);
% G22_31 = gamma(z0, Z22_31);
% 
% Z11_28=31.563439+1j*-43.602029;
% Z22_28=59.67009+1j*47.430883;
% G11_28 = gamma(z0, Z11_28);
% G22_28 = gamma(z0, Z22_28);
% 
% mx = max([20*log10(G11_31);20*log10(G22_31);20*log10(G11_28);20*log10(G22_28)],[],1);
% [m, idx] = min(mx);
% m
% z0(idx)
% G00 = 20*log10(abs((Z11_00-z0(idx))./(Z11_00+z0(idx))))


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for(i = 1:5)
%     tc = tic;
%     project = CST.Application.Active3D();
%     project.Rebuild();
%     fdsolver = project.FDSolver();
%     fdsolver.Start();
%     dt = toc(tc);
%     dispex('dt = %.5g\n', dt);
% end

%% On E-drive
% All 6 cores
% << dt = 110.09
% << dt = 109.03
% << dt = 106.48
% << dt = 106.65
% << dt = 107.71
% 5 cores
% << dt = 107.55
% << dt = 109.75
% << dt = 115.43
% << dt = 111.51
% << dt = 111.59
% 4 cores
% << dt = 114.46
% << dt = 116.65
% << dt = 113.69
% << dt = 115.18
% << dt = 116.75

%% On C-drive
% All 6 cores
% << dt = 108.93
% << dt = 109.71
% << dt = 109.45
% << dt = 108.35
% << dt = 108.8