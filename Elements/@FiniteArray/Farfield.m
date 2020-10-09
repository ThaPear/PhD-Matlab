function ff = Farfield(this, f, excitation, ths, phs, r)
    this.InitializeDs(f);
    this.InitializeZMatrix(f);

    dispex('Calculating Farfield for %i angles.\n', max(length(ths), length(phs)));
    tc = tic;


    N = length(ths);

    % If only a single theta or phi is given, repeat it to match the length of the other.
    if(length(ths) == 1)
        ths = repmat(ths, 1, length(phs));
    end
    if(length(phs) == 1)
        phs = repmat(phs, 1, length(ths));
    end

    if(any(ths < -pi/2) || any(ths > pi/2))
        error('theta outside [-pi/2, pi/2] not supported.\n');
    end

    z0 = Constants.z0;

    % Only the up-transmission line is relevant since we're calculating the pattern above
    % the array.
    tlineup_ = this.tlineup;
%     tlinedown_ = this.tlinedown;

%     ff = zeros(1, N);
%     for(ni = 1:N) % PARFOR
%         th = ths(ni);
%         ph = phs(ni);
%         
%         [k0, kx, ky, kz] = k(f, 1, th, ph);
%         
%         V = this.VoltageSpectrum(f, excitation, kx, ky);
%         
%         Ghm_xx = GreensFunction(f, k0, kx, ky, tlineup_, [], z0);
%         
%         ff(ni) = 1j .* Ghm_xx .* kz .* V .* exp(-1j .* k0 .* r) ./ (2.*pi.*r);
%     end

    [k0s, kxs, kys, kzs] = k(f, 1, ths, phs);

    V = this.VoltageSpectrum(f, excitation, kxs, kys);

    Ghm = GreensFunction3(f, k0s, kxs, kys, tlineup_, z0);
    ff_x = 1j .* Ghm.xx .* kzs .* V .* exp(-1j .* k0s .* r) ./ (2.*pi.*r);
    ff_y = 1j .* Ghm.yx .* kzs .* V .* exp(-1j .* k0s .* r) ./ (2.*pi.*r);
    ff_z = 1j .* Ghm.zx .* kzs .* V .* exp(-1j .* k0s .* r) ./ (2.*pi.*r);

    ff_theta = ff_x .* cos(ths) .* cos(phs) + ff_y .* cos(ths) .* sin(phs) - ff_z .* sin(ths);
    ff_phi   = -ff_x .* sin(phs) + ff_y .* cos(phs);

    ff = sqrt(abs(ff_theta) .^ 2 + abs(ff_phi) .^ 2);
    dt = toc(tc);
    dispex('Farfield took %.1fs for %i angles.\n', dt, max(length(ths), length(phs)));
end
function Ghm = GreensFunction3(f, k0, kx, ky, tlineup, z0)
    Vtm = 1;
    Vte = 1;
    
    kr = sqrt(kx.^2 + ky.^2);
    isTE = 1;    zteup = tlineup.GetInputImpedance(isTE, f, k0, kr);
    isTE = 0;    ztmup = tlineup.GetInputImpedance(isTE, f, k0, kr);

    iteup = 1 ./ zteup;
    itmup = 1 ./ ztmup;

    Ghm = SpectralGF.hm(z0, k0, kx, ky, Vtm, Vte, itmup, iteup);
end
