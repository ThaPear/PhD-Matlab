classdef Binomial
    methods(Static)
        function [Zs] = GetImpedances(z1, z2, N)
            C = @(n, N) factorial(N) ./ (factorial(N-n) .* factorial(n));
            
            Zs = [];
            for(n = 1:N)
                if(n == 1)
                    Zs(n) = exp(log(z1) + 2 .^ -N .* C(n-1, N) .* log(z2./z1));
                else
                    Zs(n) = exp(log(Zs(n-1)) + 2 .^ -N .* C(n-1, N) .* log(z2./z1));
                end
            end
        end
    end
end