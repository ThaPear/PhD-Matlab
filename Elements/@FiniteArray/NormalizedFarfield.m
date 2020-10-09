function ff = NormalizedFarfield(this, f, excitation, ths, phs)
    r = 1;
    ff = Farfield(this, f, excitation, ths, phs, r);

    ff = ff ./ max(abs(ff));
end
