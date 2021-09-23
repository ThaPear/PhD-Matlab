
%      o       o
%      |       |
% -   ---     ---
% |   | |     | |
% | L | |     | | er
% |   | |     | |
% -   ---     ---
%      |       |
%      o       o

classdef Line < Element
    properties
        L
        er
        mu
        cst
        % cst.LengthParameter % Store the length of the line as a parameter in CST with this name.
    end
    methods
        function this = Line(L, er, mu)
            this.L = L;
            
            if(length(er) == 1)
                % Scalar value
                this.er = er;%ones(1,3)*er;
                if(nargin > 2)
                    this.er(1:length(mu)) = er;
                end
            elseif(length(er(:)) == 3)
                % Vector value
                this.er(1) = er(1);
                this.er(2) = er(2);
                this.er(3) = er(3);
            else
                % Matrix value [erx, ..., ...
                %               ..., ery, ...
                %               ..., ..., erz]
                this.er(1) = er(1,1);
                this.er(2) = er(2,2);
                this.er(3) = er(3,3);
            end
            
            if(nargin < 3)
                this.mu = ones(size(this.er));
            elseif(length(mu) == 1)
                % Scalar value
                this.mu(1:length(this.er)) = mu;
            elseif(length(mu(:)) == 3)
                % Vector value
                this.mu(1) = mu(1);
                this.mu(2) = mu(2);
                this.mu(3) = mu(3);
            else
                % Matrix value [mux, ..., ...
                %               ..., muy, ...
                %               ..., ..., muz]
                this.mu(1) = mu(1,1);
                this.mu(2) = mu(2,2);
                this.mu(3) = mu(3,3);
            end
        end
        function Z = GetInputImpedance(this, isTE, f, k0, kr)
            error('%s::GetInputImpedance:\n\tInput impedance is not valid on a Line element.\n\tUse a shunt element or an open/shorted/terminated line.', mfilename);
        end
        function ABCD = GetABCD(this, isTE, f, k0, kr)
            % ABCD = GetABCD(this, isTE, f, k0, kr)
%             [kd, ~, ~, kzd] = k(f, obj.er, 0, 0);
            
            if(length(this.er) == 1)
                % Scalar value [er]
                kd = k0 .* sqrt(this.er);
                kzd = -1j .* sqrt(-(kd.^2 - kr.^2));
                zd = Constants.z0 ./ sqrt(this.er);
            else
                % Vector value [er.x, er.y, er.z]
                % Get incidence angle of fundamental mode.
                thi = asin(kr(round(end/2), round(end/2))./k0);
                if(isTE)
                    n = sqrt(this.er(2) .* this.mu(1) + (1-this.mu(1)./this.mu(3)).*sin(thi));
                else
                    n = sqrt(this.er(1) .* this.mu(2) + (1 - this.er(1)./this.er(3)) .* sin(thi));
                end
                kd = k0 .* n;
                kzd = -1j .* sqrt(-(kd.^2-kr.^2));
                zd = Constants.z0 ./ n;
            end
            
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