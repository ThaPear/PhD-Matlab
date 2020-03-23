function [Z] = S2Z(S, z01, z02)

s11 = S.s11;
s12 = S.s12;
s21 = S.s21;
s22 = S.s22;

const = 1 ./ ((1-s11) .* (1-s22) - s12 .* s21);

Z.z11 = const .* ((conj(z01)./z01+s11) .* (1-s22) + s12 .* s21);
Z.z12 = const .* (2 .* s12 .* sqrt(real(z01) .* real(z02) ./ (z01 .* z02)));
Z.z21 = const .* (2 .* s21 .* sqrt(real(z01) .* real(z02) ./ (z01 .* z02)));
Z.z22 = const .* ((conj(z02)./z02+s22) .* (1-s11) + s12 .* s21);

