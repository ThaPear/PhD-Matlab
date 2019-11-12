classdef SpectralGF
    methods(Static)
        function [Gmat] = matrixform(G)
            % | 1,1  1,2 |    | xx yx |
            % | 2,1  2,2 |    | xy yy |
            % | 3,1  3,2 |    | xz yz |
            
            Gmat = zeros([3, 2, size(GF.xx)]);
            Gmat(1,1,:) = GF.xx(:);
            Gmat(1,2,:) = GF.yx(:);
            Gmat(2,1,:) = GF.xy(:);
            Gmat(2,2,:) = GF.yy(:);
            Gmat(3,1,:) = GF.xz(:);
            Gmat(3,2,:) = GF.yz(:);
            
        end
        function [Gej] = ej(zeta, k, kx, ky, vtm, vte, itm, ite)
            kr = sqrt(kx.^2 + ky.^2);
            
            Gej.xx = - (vtm .* kx.^2 + vte .* ky.^2) ./ kr.^2;
            Gej.yx =   ((vte - vtm) .* kx .* ky)     ./ kr.^2;
            Gej.xy =   ((vte - vtm) .* kx .* ky)     ./ kr.^2;
            Gej.yy = - (vte .* kx.^2 + vtm .* ky.^2) ./ kr.^2;
            Gej.xz =   zeta .* (kx ./ k) .* itm;
            Gej.yz =   zeta .* (ky ./ k) .* itm;
        end
        function [Ghj] = hj(zeta, k, kx, ky, vtm, vte, itm, ite)
            kr = sqrt(kx.^2 + ky.^2);
            
            Ghj.xx = - ((ite - itm) .* kx .* ky)     ./ kr.^2;
            Ghj.yx =   (ite .* kx.^2 + itm .* ky.^2) ./ kr.^2;
            Ghj.xy =   (itm .* kx.^2 + ite .* ky.^2) ./ kr.^2;
            Ghj.yy =   ((ite - itm) .* kx .* ky)     ./ kr.^2;
            Ghj.xz =   (ky ./ (k .* zeta)) .* vte;
            Ghj.yz =   (kx ./ (k .* zeta)) .* vte;
        end
        function [Ghm] = hm(zeta, k, kx, ky, vtm, vte, itm, ite)
            kr = sqrt(kx.^2 + ky.^2);
            
            kr((kx == 0) & (ky == 0)) = 1;
            
            Ghm.xx = - (ite .* kx.^2 + itm .* ky.^2) ./ kr.^2;
            Ghm.yx = - ((ite - itm) .* kx .* ky)     ./ kr.^2;
            Ghm.xy = - ((ite - itm) .* kx .* ky)     ./ kr.^2;
            Ghm.yy = - (itm .* kx.^2 + ite .* ky.^2) ./ kr.^2;
            Ghm.xz =   (kx ./ (k .* zeta)) .* vte;
            Ghm.yz =   (ky ./ (k .* zeta)) .* vte;
        end
        function [Gem] = em(zeta, k, kx, ky, vtm, vte, itm, ite)
            kr = sqrt(kx.^2 + ky.^2);
            
            Gem.xx = - ((vte - vtm) .* kx .* ky)     ./ kr.^2;
            Gem.yx = - (vte .* ky.^2 + vtm .* kx.^2) ./ kr.^2;
            Gem.xy =   (vte .* kx.^2 + vtm .* ky.^2) ./ kr.^2;
            Gem.yy = - ((vtm - vte) .* kx .* ky)     ./ kr.^2;
            Gem.xz = - zeta .* (kx ./ k) .* itm;
            Gem.yz =   zeta .* (kx ./ k) .* itm;
        end
    end
end