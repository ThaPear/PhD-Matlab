function [kd, kxd, kyd, kzd] = k(f, er, th, ph)
    c0 = Constants.c0;
    lambda = c0./f ./ sqrt(er);
    
    % In dielectric
    kd  = 2.*pi ./ lambda;
    kxd = kd.*sin(th).*cos(ph);
    kyd = kd.*sin(th).*sin(ph);
    kzd = -1j.*sqrt(-(kd.^2 - kxd.^2 - kyd.^2));
end