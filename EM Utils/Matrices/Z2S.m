function [S] = Z2S(Z, z01, z02)

z11n = Z.z11 ./ z01;
z12n = Z.z12 ./ sqrt(z01 .* z02);
z21n = Z.z21 ./ sqrt(z01 .* z02);
z22n = Z.z22 ./ z02;

const = 1 ./ ((z11n+1) .* (z22n+1) - z12n .* z21n);

S.s11 = const .* ((z11n-conj(z01)./z01) .* (z22n+1) - z12n.*z21n);
S.s12 = const .* (2 .* z12n .* sqrt(real(z01) .* real(z02) ./ (z01 .* z02)));
S.s21 = const .* (2 .* z21n .* sqrt(real(z01) .* real(z02) ./ (z01 .* z02)));
S.s22 = const .* ((z22n-conj(z02)./z02) .* (z11n+1) - z12n.*z21n);