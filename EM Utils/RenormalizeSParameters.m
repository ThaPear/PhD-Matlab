function Snew = RenormalizeSParameters(S, z1old, z2old, z1new, z2new);

ABCD = S2ABCD(S, z1old, z2old);
Snew = ABCD2S(ABCD, z1new, z2new);