function PlotEpsilons(this, fs, hAx)
    if(nargin < 3 || isempty(hAx))
        hFig = figureex;
        hAx = hFig.CurrentAxes;
    end

    th = eps; ph = eps;
    ers = zeros(size(fs));
    parfor(fi = 1:length(fs))
        f = fs(fi);
        [epsilon, ~] = this.GetEpsilonMu(f, th, ph); %#ok<PFBNS>
        ers(fi) = epsilon.x;
    end
    f0 = 6e9;
    er0 = ers(fs == f0);
    plot(hAx, fs/1e9, (ers-er0)/er0);
%             plot(hAx, fs/1e9, ers);
end
