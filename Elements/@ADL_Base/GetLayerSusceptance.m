function Bn = GetLayerSusceptance(this, f)
    epdnext = Constants.ep0 .* this.ernext;
    epdprev = Constants.ep0 .* this.erprev;

    numM = 20;
    m = [-numM:-1, 1:numM];

    Sm = @(w)   abs(sinc(m.*w./this.p)).^2 ./ abs(m);
    fm = @(d)   -cot(-2 * 1j .* pi .* abs(m) .* d ./ this.p);
    gm = @(s,d) exp(1j*2*pi*m.*s./this.p) .* csc(-2 * 1j * pi .* abs(m) .* d ./ this.p);

    Bprevsum = Sm(this.w) .* fm(this.dprev)                       ...
             + Sm(this.wprev) .* gm(this.sprev, this.dprev);
    Bnextsum = Sm(this.w) .* fm(this.dnext)                       ...
             + Sm(this.wnext) .* gm(this.snext, this.dnext);

    Bn = 1j .* (f.*epdprev) .* this.p .* sum(Bprevsum)           ...
       + 1j .* (f.*epdnext) .* this.p .* sum(Bnextsum);
end