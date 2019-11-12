function h = GetEffectiveHeight(this, f)
    [epsilon, ~] = this.GetEpsilonMu(f, 0, 0);
    h = this.GetHeight() .* sqrt(epsilon.x);
end
