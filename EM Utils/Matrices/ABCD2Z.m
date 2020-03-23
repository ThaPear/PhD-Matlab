function [Z] = ABCD2Z(ABCD)

A = ABCD.A;
B = ABCD.B;
C = ABCD.C;
D = ABCD.D;

Z.z11 = A ./ C;
Z.z12 = (A .* D - B .* C) ./ C;
Z.z21 = 1 ./ C;
Z.z22 = D ./ C;
