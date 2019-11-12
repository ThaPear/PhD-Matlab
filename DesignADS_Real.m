% Suppress warnings:
% "For readability, consider using ones(x,y) instead of repmat(1,x,y)."
%#ok<*REPMAT>

% Design ADS with given unit cell size to have given desired epsilon.
function ads = DesignADS_Real(f0, p, L, erdes, prevlayer, nextlayer, hgap)
    if(nargin < 5)
        prevlayer = [];
    end
    if(nargin < 6)
        nextlayer = [];
    end
    wmin = Limits.wmin;
    wmax = Limits.wmaxoverp * p;
    Nw = 100;
    if(erdes < 2)
        wmin = 0.2*p;
    end
    wvals = linspace(wmin, wmax, Nw);
    
    Nmin = 1;
    Nmin = 2;
%     Nmin = 3;
    Nmax = 10;
    N = Nmin;
    
    if(erdes > 18 || erdes > 12 && erdes < 12.2)
        N = 2;
    end
    
    lambda0 = Constants.c0 / f0;
    
    % Maximum error allowed in epsilon.
    % Only used if the epsilon is too low to realize.
    maxerror = 0.25;
    
    erhost = 1;
    % Start with max shift.
    s = p/2;
    
    while(N < Nmax)
%         disp(['N = ', num2str(N)]);
        d = (L/N);
        ds = [d/2, repmat(d, 1, N-1), d/2];
        ss = repmat(s, 1, max(1, N));
        
        ers = zeros(size(wvals));
        parfor(wi = 1:Nw)
            w = wvals(wi);
            ws = repmat(w, 1, N);
            if(isempty(prevlayer))
                ws(1) = ws(1) * 2;
            end
            slab = ADS_Real(p, ds, ss, ws);
            % Introduce other layers.
%             slab = IntroduceOtherLayers(slab, prevlayer, nextlayer);
            slab = slab.SetNeighbours(prevlayer, nextlayer);
            
            [epsilon, ~] = slab.GetEpsilonMu(f0, 0, 0);
            ers(wi) = epsilon.x;
        end
%         figureex; plot(wvals/lambda0, ers, [wmin/lambda0, wmax/lambda0], [erdes, erdes]);
        % If all achieved er are too low.
        if(max(ers) < erdes)
            % Increase the number of layers.
            N = N + 1;
            continue;
        % If all achieved er are too high.
        elseif(min(ers) > erdes)
            if(N > Nmin)
                % Decrease the number of layers.
                N = N - 1;
                warning('Should not have to decrease number of layers.');
                breakpoint;
                continue;
            elseif(s > 0)
                % Remove shift.
                s = 0;
                continue;
            else
                % Minimum number of layers and no shift.
                % Check if the lowest possible er is low enough.
                if(abs(erdes - min(ers)) < maxerror)
                    warning('Accepting an epsilon error of %g - %g = %g (%g%%).',...
                        min(ers), erdes, abs(erdes - min(ers)), abs(erdes - min(ers))/erdes*100);
                    [~, wi] = min(ers);
                    w = wvals(wi);
                    ws = repmat(w, 1, N);
                    if(isempty(prevlayer))
                        ws(1) = ws(1) * 2;
                    end
                    ads = ADS_Real(p, ds, ss, ws);
                    ads = ads.SetNeighbours(prevlayer, nextlayer);
                    break;
                end
                error(['Desired er (%g) too low, lowest possible is %g.', ...
                       '\n\tIncrease dz (%.3gmm = %.3g%ceff)'], ...
                       erdes, min(ers), (L/N)*1e3, (L/N)/(lambda0/sqrt(erdes)), 955);
            end
        else
            % The desired er is in-range.
            
            
            
            
            P = InterX([wmin, wmax;     ...
                        erdes, erdes],  ...
                       [wvals;          ...
                        ers]);
%             if(size(P, 2) == 2)
%                 w = P(1, end);
%             else
                w = P(1, end);
%             end

%             figureex; plot(wvals/p, ers, [wmin/p, wmax/p], [erdes, erdes]);
%             xlabel('w / p');
%             ylabel('\varepsilon_r');
%             plot(w/p, erdes, 'ro');

%             disp(['final w = ', num2str(w/p), 'p = ', num2str(w*1e3), 'mm']);
            ws = repmat(w, 1, N);
            if(isempty(prevlayer))
                ws(1) = ws(1) * 2;
            end
            ads = ADS_Real(p, ds, ss, ws);
            ads = ads.SetNeighbours(prevlayer, nextlayer);
            break;
        end
    end
    if(~exist('ads', 'var') || isempty(ads))
        error(['Failed to design desired permittivity (%g) with up to %g layers (max is %g).', ...
              '\n\tDecrease dz (%.3gmm = %.3g%ceff)'], ...
              erdes, N, max(ers), (L/N)*1e3, (L/N)/(lambda0/sqrt(erdes)), 955);
    end
end