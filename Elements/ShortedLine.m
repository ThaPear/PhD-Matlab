classdef ShortedLine < Element
    properties
        er
        L
        erCST
    end
    methods
        function this = ShortedLine(ermatlab, L, erCST)
            if(nargin < 3)
                erCST = ermatlab;
            end
            this.er = ermatlab;
            this.L = L;
            this.erCST = erCST;
        end
        function zin = GetInputImpedance(this, isTE, f, k0, kr)
            kd = k0 .* sqrt(this.er);
            
            kzd = -1j .* sqrt(-(kd.^2 - kr.^2));
            
            [~, zcte, zctm] = z(this.er, kd, kzd);
            
            if(isTE)
                zc = zcte;
            else
                zc = zctm;
            end
            
            zin = 1j .* zc .* tan(kzd .* this.L);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
            error('%s::GetABCD:\n\tShould not be called on ShortedLine.\n\tABCD matrix is not valid.', mfilename);
        end
        function h = GetHeight(this)
            h = this.L;
        end
        function BuildCST(this, project, parentcomponent)
            if(nargin < 3 || isempty(parentcomponent))
                parentcomponent = '';
            else
                if(~strcmp(parentcomponent(end), '/'))
                    parentcomponent = [parentcomponent, '/'];
                end
            end
            componentname = [parentcomponent, 'BackingReflector'];
            project.MakeSureParameterExists('hback', this.L*1e3);
            
            % Create the metal plate.
            brick = project.Brick();
            brick.Reset();
            brick.Component(componentname);
            brick.Name('Metal');
            brick.Xrange('-dx/2', 'dx/2');
            brick.Yrange('-dy/2', 'dy/2');
            brick.Zrange('hback', 'hback');
            brick.Material('PEC');
            brick.Create();
            
            if(1)%obj.er ~= 1)
                % Create the material of the dielectric.
                materialname = num2str(this.erCST, 5);
                material = project.Material();
                material.Name(materialname);
                material.Folder('Generated');
                clr = 0.8-this.erCST/10;
                material.Colour(clr, clr, clr);
                material.Epsilon(this.erCST);
                material.Transparency(0.5);
                material.CreateConditional(['Not Material.Exists("Generated/', num2str(this.erCST), '")']);
%                 material.Create();
                
%                 % Create necessary material
%                 material.Reset();
%                 material.Name(['" & ', adlierhost, ' & "']);
%                 material.Folder('Generated');
%                 clr = ['" & 0.8-', adlierhost, '/10 & "'];
%                 material.Colour(clr, clr, clr);
%                 material.Epsilon(['" & ', adlierhost, ' & "']);
%                 material.Transparency(50);
%                 material.CreateConditional(['Not Material.Exists("Generated/" & ', adlierhost, ')']);
                
                % Create the dielectric.
                brick.Reset();
                brick.Component(componentname);
                brick.Name('Dielectric');
                brick.Xrange('-dx/2', 'dx/2');
                brick.Yrange('-dy/2', 'dy/2');
                brick.Zrange(0, 'hback');
                brick.Material(['Generated/', materialname]);
                brick.Create();
            end
        end
    end
end