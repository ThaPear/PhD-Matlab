classdef Materials
    properties(Constant)
        Rogers4350 = struct('permittivity', 3.66, ...
                            'color', [200 200 200])
        ArlonCuClad6250 = struct('thickness', 38e-6, ...
                                 'permittivity', 2.32, ...
                                 'color', [225 220 50])
        DupontPyraluxAP = struct('thickness', 25.4e-6, ...
                                 'permittivity', 3.4, ...
                                 'color', [225 50 50])
        Rohacell31HF = struct('permittivity', 1.045, ...
                              'color', [200 200 200])
        Rohacell51HF = struct('permittivity', 1.065, ...
                              'color', [200 200 200])
        Rohacell71HF = struct('permittivity', 1.093, ...
                              'color', [200 200 200])
        Rogers6002 = struct('permittivity', 2.94, ...
                            'color', [50 50 50]);
        Rogers5008 = struct('permittivity', 2.2, ...
                            'color', [50 50 50]);
        QuartzEpoxy = struct('thickness', 0.8e-3, ...
                             'permittivity', 3.3, ...
                             'color', [255, 100, 0]);
                        
        Superstrate = Materials.QuartzEpoxy
        Substrate = Materials.Rogers5008
        ADLSubstrate = Materials.DupontPyraluxAP
        Glue = Materials.ArlonCuClad6250
        Foam = Materials.Rohacell31HF
    end
    methods(Static)
        function BuildCST(project)
            material = project.Material();
            material.Folder('Generated');
            materialnames = properties(Materials);
            for(i = 1:length(materialnames))
                mat = Materials.(materialnames{i});
                name = num2str(mat.permittivity, 5);
                material.Name(name);
                material.Folder('Generated');
                material.Colour(mat.color(1), mat.color(2), mat.color(3));
                material.Epsilon(mat.permittivity);
                material.Transparency(0.5);
                material.Create();
            end
        end
    end
end