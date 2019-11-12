% Builds a frequency-dispersive material to allow approximation of this ADS
% with a homogeneous dielectric slab.
function BuildCSTMaterial(this, project)
    material = project.Material();

    material.Name('dispersive1');
    material.Folder('Generated');

    material.Type('Anisotropic');
    material.DispersiveFittingFormatEps('Real_Imag');
    material.UseGeneralDispersionEps(1);
    material.UseGeneralDispersionMu(1);
    fs = (1:0.1:10)*1e9;
    for(f = fs)
        [eps, mu] = this.GetEpsilonMu(f, 0, 0);

        material.AddDispersionFittingValueXYZEps(f/1e9, real(eps.x), round(imag(eps.x), 4), ...
                                                        real(eps.y), round(imag(eps.y), 4), ...
                                                        real(eps.z), round(imag(eps.z), 4), 1);
        material.AddDispersionFittingValueXYZMu(f/1e9,  real(mu.x), round(imag(mu.x), 4), ...
                                                        real(mu.y), round(imag(mu.y), 4), ...
                                                        real(mu.z), round(imag(mu.z), 4), 1);
    end

    material.Create();
end
