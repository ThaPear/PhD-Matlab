classdef Vias
    properties(SetAccess = protected)
        via_distance
        via_radius
    end
    methods
        function this = Vias(via_distance, via_radius)
            this.via_distance = via_distance;
            this.via_radius = via_radius;
        end
        
        function BuildCST(this, project, parentcomponent)
            if(nargin < 3 || isempty(parentcomponent))
                parentcomponent = '';
            else
                if(~strcmp(parentcomponent(end), '/'))
                    parentcomponent = [parentcomponent, '/'];
                end
            end
            componentname = [parentcomponent, 'Vias'];
            
            brick = project.Brick();
            cylinder = project.Cylinder();
            transform = project.Transform();
            solid = project.Solid();
            
            if(isnumeric(this.via_distance))
                project.StoreParameter('via_distance', this.via_distance*1e3);
            else
                project.StoreParameter('via_distance', this.via_distance);
            end
            if(isnumeric(this.via_radius))
                project.StoreParameter('via_radius', this.via_radius*1e3);
            else
                project.StoreParameter('via_radius', this.via_radius);
            end
            
            cylinder.Reset();
            cylinder.Component(componentname);
            cylinder.Name('Via');
            cylinder.Axis('z');
            cylinder.Material('PEC');
            cylinder.Xcenter(0);
            cylinder.Ycenter(['-via_distance']);
            cylinder.Zrange('-hback', 0);
            cylinder.Innerradius(0);
            cylinder.Outerradius('via_radius');
            cylinder.Create();
            
            transform.Reset();
            transform.Name([componentname, ':Via']);
            transform.Angle(0, 0, 90);
            transform.Origin('Free');
            transform.Center('0', '-dy/2', 0);
            transform.Repetitions(3);
            transform.MultipleObjects(1);
            transform.Transform('Shape', 'Rotate');
            
            transform.Reset();
            transform.Name([componentname, ':Via']);
            transform.AddName([componentname, ':Via_1']);
            transform.AddName([componentname, ':Via_2']);
            transform.AddName([componentname, ':Via_3']);
            transform.Vector('dx', 0, 0);
            transform.Repetitions(1);
            transform.MultipleObjects(1);
            transform.Transform('Shape', 'Translate');
            
            transform.Reset();
            transform.Name([componentname, ':Via']);
            transform.AddName([componentname, ':Via_1']);
            transform.AddName([componentname, ':Via_2']);
            transform.AddName([componentname, ':Via_3']);
            transform.Vector('dx', 'dy', 0);
            transform.Repetitions(1);
            transform.MultipleObjects(1);
            transform.Transform('Shape', 'Translate');
            
            transform.Reset();
            transform.Name([componentname, ':Via']);
            transform.AddName([componentname, ':Via_1']);
            transform.AddName([componentname, ':Via_2']);
            transform.AddName([componentname, ':Via_3']);
            transform.Vector(0, 'dy', 0);
            transform.Repetitions(1);
            transform.MultipleObjects(1);
            transform.Transform('Shape', 'Translate');
            
            brick.Reset();
            brick.Component(componentname);
            brick.Name('Substrate');
            brick.Xrange('-dx/2', 'dx/2');
            brick.Yrange('-dy/2', 'dy/2');
            brick.Zrange('-hback', 0);
            brick.Create();
            
            brick.Reset();
            brick.Component(componentname);
            brick.Name('InverseSubstrate');
            brick.Xrange('-dx*2', 'dx*2');
            brick.Yrange('-dy*2', 'dy*2');
            brick.Zrange('-hback', 0);
            brick.Create();
            
            solid.Subtract([componentname, ':InverseSubstrate'], [componentname, ':Substrate']);
            
            for(i = 1:6)
                if(i < 4)
                    for(j = 1:3)
                        solid.Insert([componentname, ':Via_', num2str(i), '_', num2str(j)], [componentname, ':InverseSubstrate']);
                    end
                end
                solid.Insert([componentname, ':Via_', num2str(i)], [componentname, ':InverseSubstrate']);
            end
            solid.Subtract([componentname, ':Via'], [componentname, ':InverseSubstrate']);
        end
    end
end