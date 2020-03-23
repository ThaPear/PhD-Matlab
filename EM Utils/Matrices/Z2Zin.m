function [Zin] = Z2Zin(Z, zL)

z11 = Z.z11;
z12 = Z.z12;
z21 = Z.z21;
z22 = Z.z22;

Zin = z11 - z21 .* z12 ./ (z22 + zL);