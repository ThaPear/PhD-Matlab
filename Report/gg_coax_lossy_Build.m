ga_ms_Build;

project.StoreParameter('nsamplesperGHz', 1);

material = project.Material();
solid = project.Solid();

%% Define copper (Copied from history list)
material.Reset();
material.Name('Copper (pure)');
material.Folder('');
material.FrqType('all');
material.Type('Lossy metal');
material.MaterialUnit('Frequency', 'GHz');
material.MaterialUnit('Geometry', 'mm');
material.MaterialUnit('Time', 's');
material.MaterialUnit('Temperature', 'Kelvin');
material.Mu('1.0');
material.Sigma('5.96e+007');
material.Rho('8930.0');
material.ThermalType('Normal');
material.ThermalConductivity('401.0');
material.HeatCapacity('0.39');
material.MetabolicRate('0');
material.BloodFlow('0');
material.VoxelConvection('0');
material.MechanicsType('Isotropic');
material.YoungsModulus('120');
material.PoissonsRatio('0.33');
material.ThermalExpansionRate('17');
material.ReferenceCoordSystem('Global');
material.CoordSystemType('Cartesian');
material.NLAnisotropy('False');
material.NLAStackingFactor('1');
material.NLADirectionX('1');
material.NLADirectionY('0');
material.NLADirectionZ('0');
material.FrqType('static');
material.Type('Normal');
material.SetMaterialUnit('Hz', 'mm');
material.Epsilon('1');
material.Mu('1.0');
material.Kappa('5.96e+007');
material.TanD('0.0');
material.TanDFreq('0.0');
material.TanDGiven('False');
material.TanDModel('ConstTanD');
material.KappaM('0');
material.TanDM('0.0');
material.TanDMFreq('0.0');
material.TanDMGiven('False');
material.TanDMModel('ConstTanD');
material.DispModelEps('None');
material.DispModelMu('None');
material.DispersiveFittingSchemeEps('Nth Order');
material.DispersiveFittingSchemeMu('Nth Order');
material.UseGeneralDispersionEps('False');
material.UseGeneralDispersionMu('False');
material.Colour('1', '1', '0');
material.Wireframe('False');
material.Reflection('False');
material.Allowoutline('True');
material.Transparentoutline('False');
material.Transparency('0');
material.Create();

%% Add losses to 2.2
material.Reset();
material.Name('2.2');
material.Folder('Generated');
clr = [50 50 50]./256;
material.Colour(clr(1), clr(2), clr(3));
material.Epsilon(2.2);
material.Transparency(0.5);
material.TanD('0.0009');
material.TanDGiven('True');
material.Create();

%% Change all PEC to Copper
objects = {'ADS1:Metal', 'ADS2:Metal', 'BackingReflector:Metal', 'Slot:Metal', 'Vias:Via_2', 'Vias:Via_2_1', 'Vias:Via_3', 'Vias:Via_3_3'};
feedobjects = {'Core', 'Core2ms', 'Microstrip', 'Microstrip_Buried', 'Patch', 'Shield', 'Shield_1', 'Shield_2', 'SlotFeed'};
objects = [objects, strcat('Feed/Xfeed:', feedobjects), strcat('Feed/Yfeed:', feedobjects)];

for(i = 1:length(objects))
    solid.ChangeMaterial(objects{i}, 'Copper (pure)');
end



