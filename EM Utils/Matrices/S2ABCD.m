function [ABCD] = S2ABCD(S, z01, z02)

s11 = S.s11;
s12 = S.s12;
s21 = S.s21;
s22 = S.s22;

const = 1 ./ (2 .* s21 .* sqrt(real(z01) .* real(z02)));

ABCD.A = const .* ( (conj(z01) + s11 .* z01) .* (1 - s22) + s12 .* s21 .* z01                       );
ABCD.B = const .* ( (conj(z01) + s11 .* z01) .* (conj(z02) + s22 .* z02) - s12 .* s21 .* z01 .* z02 );
ABCD.C = const .* ( (1 - s11) .* (1 - s22) - s12 .* s21                                             );
ABCD.D = const .* ( (1 - s11) .* (conj(z02) + s22 .* z02) + s12 .* s21 .* z02                       );

