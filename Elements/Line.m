classdef Line < Element
    properties
        er
        L
        CSTLengthParameter
    end
    methods
        function this = Line(er, L)
            this.er = er;
            this.L = L;
        end
        function Z = GetInputImpedance(this, isTE, f, k0, kr)
            error('%s::GetInputImpedance:\n\tInput impedance is not valid on a Line element.\n\tUse a shunt element or an open/shorted/terminated line.', mfilename);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
%             [kd, ~, ~, kzd] = k(f, obj.er, 0, 0);
            kd = k0 .* sqrt(this.er);
            
            kzd = -1j .* sqrt(-(kd.^2 - kr.^2));
            
            % Since kzd can be 0, the sin(kzd*L)/(kzd*L) are combined into
            % sinc such that they result in 1 instead of NaN.
            % Replaced code was kept for clarity.
            % [~, zcte, zctm] = z(this.er, kd, kzd);
            % 
            % if(isTE)
            %     zc = zcte;
            % else
            %     zc = zctm;
            % end
            % 
            % ABCD = ABCDMatrix(                                          ...
            %     cos(kzd.*this.L),                                        ...
            %     1j.*zc.*sin(kzd.*this.L),                                ...
            %     1j./zc.*sin(kzd.*this.L),                                ...
            %     cos(kzd.*this.L)                                         ...
            % );
            % Replacement code below.
            zd = Constants.z0 ./ sqrt(this.er);
            if(isTE)
                ABCD = ABCDMatrix(                                      ...
                    cos(kzd.*this.L),                                    ...
                    1j .* zd .* kd .* this.L .* sinc(kzd.*this.L./pi),    ... % Divide by pi for matlab
                    1j ./ (zd .* kd ./ kzd) .* sin(kzd.*this.L),         ...
                    cos(kzd.*this.L)                                     ...
                );
            else
                ABCD = ABCDMatrix(                                      ...
                    cos(kzd.*this.L),                                    ...
                    1j .* (zd .* kzd ./ kd) .* sin(kzd.*this.L),         ...
                    1j .* kd ./ zd .* this.L .* sinc(kzd.*this.L./pi),    ... % Divide by pi for matlab
                    cos(kzd.*this.L)                                     ...
                );
            end
        end
        function h = GetHeight(this)
            h = this.L;
        end
        function h = GetEffectiveHeight(this, f)
            h = this.L .* sqrt(this.er);
        end
        function BuildCST(this, project)
            if(~isempty(this.CSTLengthParameter))
                project.StoreParameter(this.CSTLengthParameter, this.L*1e3);
            end
            
            if(this.er ~= 1)
                % Create necessary material
                material = project.Material();
                material.Reset();
                materialname = num2str(this.er, 5);
                material.Name(materialname);
                material.Folder('Generated');
                material.Colour(0, min(1, this.er/20), 1);
                material.Epsilon(this.er);
                material.Transparency(0.5);
                material.Create();


                % Create the dielectric.
                brick = project.Brick();
                brick.Reset();
                brick.Component('Lines');
                brickname = ['Line ', num2str(floor(rand()*1e4))];
                brick.Name(brickname);
                brick.Xrange('-dx/2', 'dx/2');
                brick.Yrange('-dy/2', 'dy/2');
                if(~isempty(this.CSTLengthParameter))
                    brick.Zrange(0, this.CSTLengthParameter);
                else
                    brick.Zrange(0, this.L*1e3);
                end
                brick.Material(['Generated/', materialname]);
                brick.Create();
                
                solid = project.Solid();
                solid.MergeMaterialsOfComponent(['Lines:', brickname]);
            end
            
            wcs = project.WCS();
            if(~isempty(this.CSTLengthParameter))
                wcs.MoveWCS('local', 0, 0, this.CSTLengthParameter);
            else
                wcs.MoveWCS('local', 0, 0, this.L*1e3);
            end
        end
    end
end