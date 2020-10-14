SetupPath;
close all;
clear;

% Assumptions:
% Power on vertical slot:   A
% Power on horizontal slot: C
% Crosspol of vertical slot: a
% Crosspol of horizontal slot: b
% Horizontal: H = aA     + (1-d)C
% Vertical:   V = (1-a)A - dC

for(i = 1:2)
    maxP = [];
    if(i == 1)
        xpols = -20:5:-10; % dB
    else
        xpols = -50:0.1:0; % dB
    end
    for(xpol = xpols)
        a = 10^(xpol/20);
        d = 10^(xpol/20);

        % Linear polarization with angle th
        G = 1; % Amplitude.

        th = (0:359) * pi/180;
        g = cos(th) ./ sin(th);

        H = G .* cos(th);
        V = G .* sin(th);

        A = (H + V.*(1-d)./d) ./ (a + (1-d).*(1-a)./d);
        C = ((1-a).*A-V)./d;

        H = a.*A + (1-d).*C;
        V = (1-a).*A - d.*C;



        if(length(xpols) <= 5)
            fig = figureex(1);
                ax = fig.CurrentAxes;
                ax.ColorOrder = reshape(repmat(lines(7), 1, 3).', [], 7*3).';
                plot(th * 180/pi, A);
                addlegendentry([num2str(xpol)]);
                plot(th * 180/pi, C);
                plot(th * 180/pi, abs(A).^2+abs(C).^2);
                ax.XTick = 0:30:360;
            fig = figureex(2);
                ax = fig.CurrentAxes;
                plot(th * 180/pi, H);
                plot(th * 180/pi, V);
                plot(th * 180/pi, sqrt(H.^2+V.^2));
                legend({'H', 'V', 'G'});
                ax.XTick = 0:30:360;
        end
        maxP = [maxP, max(abs(A).^2+abs(C).^2)];
    end
end
figureex;
    plot(xpols, maxP);
    xlabel('Crosspol [dB]');
    ylabel('Power Ratio');