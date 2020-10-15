clear;
close all;

% [p, dat] = CST.LoadData('h:\Server\Exports\PairsBetweenSlots_60_-45.txt');
% 
% figureex;
% plot(p{1}.frequencies, 20*log10(abs(dat{1})));
% plot(p{2}.frequencies, 20*log10(abs(dat{2})));
% 
% xlim([12 32]);
% ylim([-30 0]);
% legend({'S_{11}', 'S_{22}'})
% xlabel('Frequency [GHz]'); ylabel('S_{ii}');




th0s = [0 60];
ph0s = [-45, 0, 45];%, 90];
for(th0 = th0s)
    for(ph0 = ph0s)
        if(th0 == 0 && ph0 ~= 0)
            continue;
        end

        filename = sprintf('h:\\Server\\Exports\\PairsBetweenSlots_%i_%i_Z.txt', th0, ph0);
%         filename = sprintf('h:\\Server\\Exports\\PairsAlongSlots_%i_%i_Z.txt', th0, ph0);
        if(~exist(filename, 'file'))
            continue;
        end
        [p, dat] = CST.LoadData(filename);
        f = p{1}.frequencies;
        z1 = dat{1};
        z2 = dat{2};



        zfeed = 100;
        S1 = (z1 - zfeed) ./ (z1 + zfeed);
        S2 = (z2 - zfeed) ./ (z2 + zfeed);
        zsum = z1 + z2;
        zfeedsum = 200;
        S = (zsum - zfeedsum) ./ (zsum + zfeedsum);

        [hFig, hAx] = figureex;
            if(length(hAx.Children) < 2)
                patch(hAx, [28 31 31 28], [-1000 -1000 1000 1000], [0 0 0], ...
                    'FaceAlpha', 0.1, 'EdgeColor', 'none');
                patch(hAx, [13.75 14.5 14.5 13.75], [-1000 -1000 1000 1000], [0 0 0], ...
                    'FaceAlpha', 0.1, 'EdgeColor', 'none');
            end
            repeatcolormap(2);
            
            plot(f, real(z1));
            addlegendentry('Z1');
            plot(f, imag(z1), '--');
            
            plot(f, real(z2));
            addlegendentry('Z2');
            plot(f, imag(z2), '--');
            
            plot(f, real(zsum));
            addlegendentry('Z1 + Z2');
            plot(f, imag(zsum), '--');
            
            plot(nan, nan, 'k');
            addlegendentry('Real');
            plot(nan, nan, 'k--');
            addlegendentry('Imag');
            
            xlabel('Frequency [GHz]');
            ylabel('Z_{act} [\Omega]');
            title(sprintf('\\theta = %i\\circ, \\phi = %i\\circ', th0, ph0));
            xlim([12 32]);
            ylim([-100 350]);
            legend(hAx, 'NumColumns', 5);
            legendlinelength(hAx, 13);
            alignplot(hFig, 5, 3, hFig.Number, [], 2);
            movelegend(hAx, 'se');
%         figureex;
%             repeatcolormap(2);
%             xlabel('Frequency [GHz]');
%             ylabel('Z_{act} [\Omega]');
%             title(sprintf('\\theta = %i\\circ, \\phi = %i\\circ', th0, ph0));
%             xlim([12 32]);
%             ylim([-30 0]);

        [hFig, hAx] = figureex(2);
            if(length(hAx.Children) < 2)
                patch(hAx, [28 31 31 28], [-1000 -1000 1000 1000], [0 0 0], ...
                    'FaceAlpha', 0.1, 'EdgeColor', 'none');
                patch(hAx, [13.75 14.5 14.5 13.75], [-1000 -1000 1000 1000], [0 0 0], ...
                    'FaceAlpha', 0.1, 'EdgeColor', 'none');
            end
            
%             plot(f, 20*log10(abs(S1)));
%             addlegendentry('S_{act,11}');
%             
%             plot(f, 20*log10(abs(S2)));
%             addlegendentry('S_{act,22}');
            
            plot(f, 20*log10(abs(S)));
            addlegendentry(sprintf('\\theta = %i\\circ, \\phi = %i\\circ', th0, ph0))
%             addlegendentry('S_{act,sum}');
            
            xlabel('Frequency [GHz]');
            ylabel('S_{act} [dB]');
%             title(sprintf('\\theta = %i\\circ, \\phi = %i\\circ', th0, ph0));
            xlim([12 32]);
            ylim([-30 0]);
            legend(hAx, 'NumColumns', 2);
            legendlinelength(hAx, 15);
            alignplot(hFig, 5, 3, hFig.Number, [], 2);
            drawnow;
%             movelegend(hAx, 'ne');
    end
end

%%
for(i = 1:100)
    if(ishandle(i))
        hFig = figure(i);
        hAx = hFig.CurrentAxes;
        alignplot(hFig, 5, 3, hFig.Number, [], 2);
        for(iLine = 1:length(hAx.Children))
            if(isa(hAx.Children(iLine), 'matlab.graphics.chart.primitive.Line'))
                hAx.Children(iLine).LineWidth = 1;
            end
        end
    end
end