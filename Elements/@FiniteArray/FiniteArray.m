classdef FiniteArray < handle
    properties
        unitcell  % Unit cell in the finite array
        
        tlineup   % Stratification above the array
        tlinedown % Stratification below the array
        
        Nx % Number of unit cells in x
        Ny % Number of unit cells in y
        dedge % Distance to the slot termination
        zfeed % Impedance of the feeding network
        
        Zmat % Mutual impedance matrix
        Z_fs % Frequencies for impedance matrix
        
        Ds    % Precomputed Ds
        D_kxs % Kx values for precomputed Ds
        D_fs  % Frequencies for precomputed Ds
        D_interpolants % Interpolant to index at any given kx
        
        KyInts    % Precomputed ky integrals
        KyInt_kxs % Kx values for precomputed ky integrals
        KyInt_fs  % Frequencies for precomputed ky integrals
        KyInt_interpolants % Interpolant to index at any given kx
    end
    methods
        function this = FiniteArray(unitcell, tlineup, tlinedown, Nx, Ny, dedge, zfeed)
            this.unitcell = unitcell;
            this.tlineup = tlineup;
            this.tlinedown = tlinedown;
            this.Nx = Nx;
            this.Ny = Ny;
            this.dedge = dedge;
            this.zfeed = zfeed;
        end
        InitializeDs(this, fs);
        InitializeKyInts(this, fs);
        InitializeZMatrix(this, fs);
        Zas = GetInputImpedance(this, fs, excitation);
        Zred = ReducedZMatrix(this, fs);
        M = VoltageSpectrum(this, fs, excitation, kx, ky);
        [x, v] = Voltage(this, fs, excitation, ny, Npoints);
        ff = Farfield(this, f, excitation, ths, phs, r);
        ff = NormalizedFarfield(this, f, excitation, ths, phs);
        
        BuildCST(this, project, parentcomponent);
    end
    methods(Static)
        kvec = UnfoldKVector(kvec, integrationpath);
    end
end







%     pltx = kx;
%     plty = v;
%     [hFig, hAx] = figureex;
%         title(hAx, 'v');
%         hAx.ColorOrder = lines(min(size(plty)));
%         plot(real(pltx), real(plty));
%         plot(real(pltx), imag(plty), '--');
%         plot(real(pltx), abs(plty), ':');
%         alignplot(hFig, 5, 3, [], hFig.Number, 1);