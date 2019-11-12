function [S] = ABCD2S(ABCD, z01, z02)

if(nargin < 2)
    z01 = Constants.z0;
end
if(nargin < 3)
    z02 = Constants.z0;
end

A = ABCD.A;
B = ABCD.B;
C = ABCD.C;
D = ABCD.D;

% z0 = z01;
% 
% S.s11 = (A + B./z0 - C.*z0 - D)  ./ (A + B./z0 + C.*z0 + D);
% S.s12 = (2*(A.*D - B.*C))        ./ (A + B./z0 + C.*z0 + D);
% S.s21 = 2                        ./ (A + B./z0 + C.*z0 + D);
% S.s22 = (-A + B./z0 - C.*z0 + D) ./ (A + B./z0 + C.*z0 + D);

const = 1 ./ (A .* z02 + B + C .* z01 .* z02 + D .* z01);

S.s11 = const .* ( A .* z02 + B - C .* conj(z01) .* z02 - D .* conj(z01)  );
S.s12 = const .* ( 2 .* (A.*D - B.*C) .* sqrt(real(z01) .* real(z02))     );
S.s21 = const .* ( 2 .* sqrt(real(z01) .* real(z02))                      );
S.s22 = const .* ( -A .* conj(z02) + B - C .* z01 .* conj(z02) + D .* z01 );


