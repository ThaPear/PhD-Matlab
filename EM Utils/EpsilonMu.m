function [er, mu] = EpsilonMu(f, Ste, Stm, Ste0, Stm0, height, th0, th, ph)
%     printvar th*180/pi;
%     printvar th0*180/pi;
    d = height;
    
    k0 = 2*pi*f/3e8;
    k00 = 2*pi*f/3e8;
    
    etaTE  = sqrt(((1+Ste.s11 ).^2 - Ste.s21.^2 ) ./ ((1-Ste.s11 ).^2 - Ste.s21.^2 )) .* sec(th );
    etaTE0 = sqrt(((1+Ste0.s11).^2 - Ste0.s21.^2) ./ ((1-Ste0.s11).^2 - Ste0.s21.^2)) .* sec(th0);

    etaTM  = sqrt(((1+Stm.s11 ).^2 - Stm.s21.^2 ) ./ ((1-Stm.s11 ).^2 - Stm.s21.^2 )) .* cos(th );
    etaTM0 = sqrt(((1+Stm0.s11).^2 - Stm0.s21.^2) ./ ((1-Stm0.s11).^2 - Stm0.s21.^2)) .* cos(th0);

    if(real(etaTE)  < 0); etaTE  = -etaTE;  end
    if(real(etaTM)  < 0); etaTM  = -etaTM;  end
    if(real(etaTE0) < 0); etaTE0 = -etaTE0; end
    if(real(etaTM0) < 0); etaTM0 = -etaTM0; end

    zetaTE  = Ste.s21  ./ (1 - Ste.s11 .*((etaTE .*cos(th ) - 1) ./ (etaTE  .* cos(th ) + 1)));
    zetaTE0 = Ste0.s21 ./ (1 - Ste0.s11.*((etaTE0.*cos(th0) - 1) ./ (etaTE0 .* cos(th0) + 1)));

    zetaTM  = Stm.s21  ./ (1 - Stm.s11 .*((etaTM ./cos(th ) - 1) ./ (etaTM  ./ cos(th ) + 1)));
    zetaTM0 = Stm0.s21 ./ (1 - Stm0.s11.*((etaTM0./cos(th0) - 1) ./ (etaTM0 ./ cos(th0) + 1)));

    nTE  = @(m) sqrt( ( (log(abs(zetaTE )) + 1j.*(angle(zetaTE ) + 2*pi*m)) ./ (-1j.*k0 .*d)).^2 + sin(th ).^2);
    nTE0 = @(m) sqrt( ( (log(abs(zetaTE0)) + 1j.*(angle(zetaTE0) + 2*pi*m)) ./ (-1j.*k00.*d)).^2 + sin(th0).^2);

    nTM  = @(m) sqrt(((log(abs(zetaTM )) + 1j.*(angle(zetaTM ) + 2*pi*m)) ./ (-1j.*k0 .*d)).^2 + sin(th ).^2);
    nTM0 = @(m) sqrt(((log(abs(zetaTM0)) + 1j.*(angle(zetaTM0) + 2*pi*m)) ./ (-1j.*k00.*d)).^2 + sin(th0).^2);


%     eta  = @(S, th)     sqrt(((1+S.s11).^2 - S.s21.^2) ./ ((1-S.s11).^2 - S.s21.^2)) .* sec(th);
%     zeta = @(S, th)     S.s21 ./ (1 - S.s11.*((eta(S, th).*cos(th) - 1) ./ (eta(S, th) .* cos(th) + 1)));
%     n    = @(S, th, m)  sqrt(((log(abs(zeta(S, th))) + 1j.*angle(zeta(S, th)) + 2*pi*m) ./ (-1j.*k0.*d)).^2) + sin(th).^2;

    % Epsilon as function of m.
    erx  = @(m)         nTM0(m) ./ etaTM0;
    ery  = @(m)         nTE0(m) ./ etaTE0;
    erz  = @(m)         erx(m) .* sin(th).^2 ./ (sin(th).^2 - nTM(m).^2 + nTM0(m).^2);
    erz2 = @(m)         sin(th).^2 ./ (nTM0(m) .* etaTM0 - erx(m) .* etaTM.^2);

    mux  = @(m)         nTE0(m) .* etaTE0;
    muy  = @(m)         nTM0(m) .* etaTM0;
    muz  = @(m)         sin(th).^2 ./ ((nTE0(m) ./ etaTE0) - mux(m) ./ etaTE.^2);
    muz2 = @(m)         mux(m) .* sin(th).^2 ./ (sin(th).^2 - nTE(m).^2 + nTE0(m).^2);

%     ms = -50:1:50;
% 
%     errors = abs(erz(ms) - erz2(ms));
%     [~, idx] = min(errors); % Get index of lowest value
%     mbar = ms(idx);         % Get corresponding m.
% 
%     errors = abs(muz(ms) - muz2(ms));
%     [~, idx] = min(errors); % Get index of lowest value
%     m = ms(idx);            % Get corresponding m.

    % The error minimization doesn't work correctly.
    % Use m = mbar = 0 for all calculations.
    mbar = 0;
    m = 0;

    er.x = erx(mbar);
    er.y = ery(mbar);
    er.z = erz(mbar);
    mu.x = mux(m);
    mu.y = muy(m);
    mu.z = muz2(m);
    
    er.x = real(er.x);
    er.y = real(er.y);
    er.z = real(er.z);
    mu.x = real(mu.x);
    mu.y = real(mu.y);
    mu.z = real(mu.z);
end