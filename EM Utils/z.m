function [zd, zdte, zdtm] = z(er, kd, kzd)

zd   = Constants.z0 ./ sqrt(er);
zdte = zd .* kd ./ kzd;
zdtm = zd .* kzd ./ kd;