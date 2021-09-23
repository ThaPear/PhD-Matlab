function ff = NormalizedFarfield(this, f, excitation, ths, phs)
    r = 1;
    ff = this.Farfield(f, excitation, ths, phs, r);

    ff = ff ./ max(abs(ff));
end
