function [ABCD] = Z2ABCD(Z, z01, z02)

S = Z2S(Z, z01, z02);
ABCD = S2ABCD(S, z01, z02);

