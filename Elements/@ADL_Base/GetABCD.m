function ABCD = GetABCD(this, isTE, f, k0, kr)
    Bn = this.GetLayerSusceptance(f);

    if(isTE)
        Y = 1j .* (Bn .* (1-kr.^2./k0.^2./2));
    else
        Y = 1j .* Bn;
    end

    ABCD = ABCDMatrix(1, 0, ...
                      Y, 1);

end