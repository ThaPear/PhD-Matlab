
%      o       o
%      |       |
% -   ---     ---
% |   | |     | |
% | L | |     | | er
% |   | |     | |
% -   ---     ---
%      |       |
%      o       o

classdef TXLine < Element
    properties
        L
        z0
        % cst.LengthParameter % Store the length of the line as a parameter in CST with this name.
    end
    methods
        function this = TXLine(L, z0)
            this.L = L;
            this.z0 = z0;
        end
        function Z = GetInputImpedance(this, isTE, f, k0, kr)
            error('%s::GetInputImpedance:\n\tInput impedance is not valid on a Line element.\n\tUse a shunt element or an open/shorted/terminated line.', mfilename);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
            kline = k0;
            ABCD = ABCDMatrix(                                          ...
                cos(kline.*this.L),                                     ...
                1j.*this.z0.*sin(kline.*this.L),                        ...
                1j./this.z0.*sin(kline.*this.L),                        ...
                cos(kline.*this.L)                                      ...
            );
        end
        function h = GetHeight(this)
            h = this.L;
        end
        function h = GetEffectiveHeight(this, f)
            h = this.L .* sqrt(this.er);
        end
        function BuildCST(this, project, parentcomponent)
            if(nargin < 3 || isempty(parentcomponent))
                parentcomponent = '';
            else
                if(~strcmp(parentcomponent(end), '/'))
                    parentcomponent = [parentcomponent, '/'];
                end
            end
            componentname = [parentcomponent, 'Lines'];
            if(isfield(this.cst, 'LengthParameter'))
                project.StoreParameter(this.cst.LengthParameter, this.L*1e3);
            end
            
            if(length(this.er) > 1 || this.er ~= 1 || length(this.mu) > 1 || this.mu ~= 1)
                anisotropic = (length(unique(this.er)) == 1) && (length(unique(this.mu)) == 1);
                % Create necessary material
                material = project.Material();
                material.Reset();
                if(anisotropic)
                    er_ = repmat(this.er, 1, 4-length(this.er)); % Ensure er is 1 by 3.
                    mu_ = repmat(this.mu, 1, 4-length(this.mu)); % Ensure mu is 1 by 3.
                    materialname = [num2str(er_(1), 5), '_', num2str(er_(2), 5), '_', num2str(er_(3), 5)];
                    material.Type('Anisotropic');
                    material.EpsilonX(er_(1));
                    material.EpsilonY(er_(2));
                    material.EpsilonZ(er_(3));
                    material.MuX(mu_(1));
                    material.MuY(mu_(2));
                    material.MuZ(mu_(3));
                else
                    material.Epsilon(this.er(1));
                    material.Mu(this.mu(1));
                    materialname = num2str(this.er(1));
                end
                material.Name(materialname);
                material.Folder('Generated');
                material.Colour(0, min(1, this.er(1)/20), 1);
                material.Transparency(0.5);
                material.Create();


                % Create the dielectric.
                brick = project.Brick();
                brick.Reset();
                brick.Component(componentname);
                brickname = ['Line ', num2str(floor(rand()*1e4))];
                brick.Name(brickname);
                brick.Xrange('-dx/2', 'dx/2');
                brick.Yrange('-dy/2', 'dy/2');
                if(isfield(this.cst, 'LengthParameter'))
                    brick.Zrange(0, this.cst.LengthParameter);
                else
                    brick.Zrange(0, this.L*1e3);
                end
                brick.Material(['Generated/', materialname]);
                brick.Create();
                
                solid = project.Solid();
                solid.MergeMaterialsOfComponent([componentname, ':', brickname]);
            end
            
            wcs = project.WCS();
            if(isfield(this.cst, 'LengthParameter'))
                wcs.MoveWCS('local', 0, 0, this.cst.LengthParameter);
            else
                wcs.MoveWCS('local', 0, 0, this.L*1e3);
            end
        end
    end
end