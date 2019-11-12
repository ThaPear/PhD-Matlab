classdef Constants
    properties(Constant)
        c0 = 3e8;%299792458;
        ep0 = 8.854e-12;%Constants.e^2 / (2*Constants.h*Constants.c0*Constants.alpha);%1/(Constants.mu0*Constants.c0^2); % 8.854187817620389e-12
                                 % Or, with c=3e8 8.841941282883074e-12
        % Vacuum permeability mu0 = 4*pi*1.00000000082*1e-7;
        mu0 = 4*pi*1e-7;%2*Constants.h*Constants.alpha/(Constants.c0*Constants.e^2);
        
        z0 = 120*pi;%2*Constants.h*Constants.alpha/Constants.e^2;%sqrt(Constants.mu0/Constants.ep0); % 120*pi
        h = 6.62607015e-34; % Planck
        e = 1.602176634e-19; % Elementary Charge
        kb = 1.380649e-23; % Boltzmann
        Na = 6.02214076e23; % Avogadro
        % Fine structure constant
        alpha = 0.007297352571289;%Constants.c0*Constants.mu0/(2*Constants.RK)%Constants.e^2 / (4*pi) * Constants.z0 / (Constants.h/2*pi); 
        RK = Constants.h / Constants.e^2; % von Klitzing constant
    end
end
