classdef Cavity
    properties(SetAccess = protected)
        cavity_width
        cavity_width_diagonal
    end
    methods
        function this = Cavity(cavity_width, cavity_width_diagonal)
            this.cavity_width = cavity_width;
            this.cavity_width_diagonal = cavity_width_diagonal;
        end
        
        function BuildCST(this, project, parentcomponent)
            if(nargin < 3 || isempty(parentcomponent))
                parentcomponent = '';
            else
                if(~strcmp(parentcomponent(end), '/'))
                    parentcomponent = [parentcomponent, '/'];
                end
            end
            componentname = [parentcomponent, 'Cavity'];
            
            brick = project.Brick();
            transform = project.Transform();
            solid = project.Solid();
            
            if(isnumeric(this.cavity_width))
                project.StoreParameter('cavity_width', this.cavity_width*1e3);
            else
                project.StoreParameter('cavity_width', this.cavity_width);
            end
            if(isnumeric(this.cavity_width_diagonal))
                project.StoreParameter('cavity_width_diagonal', this.cavity_width_diagonal*1e3);
            else
                project.StoreParameter('cavity_width_diagonal', this.cavity_width_diagonal);
            end
            project.MakeSureParameterExists('dms', 0.127);
            
            brick.Reset();
            brick.Component(componentname);
            brick.Name('Cavity');
            brick.Xrange('-dx/2', 'dx/2');
            brick.Yrange('-cavity_width/2', 'cavity_width/2');
            brick.Zrange('-hback', '-dms');
            brick.Material('PEC');
            brick.Create();
            
            brick.Reset();
            brick.Component(componentname);
            brick.Name('CavityDiag');
            brick.Xrange('dx/2-cavity_width_diagonal/2', 'dx/2+cavity_width_diagonal/2');
            brick.Yrange('-cavity_width_diagonal/2', 'cavity_width_diagonal/2');
            brick.Zrange('-hback', '-dms');
            brick.Material('PEC');
            brick.Create();
        
            transform.Reset();
            transform.Name([componentname, ':CavityDiag']);
            transform.Angle(0, 0, 45);
            transform.Origin('ShapeCenter');
            transform.Repetitions(1);
            transform.Transform('Shape', 'Rotate');
        
            transform.Reset();
            transform.Name([componentname, ':CavityDiag']);
            transform.Vector('-dx', 0, 0);
            transform.Repetitions(1);
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Transform('Shape', 'Translate');
        
            transform.Reset();
            transform.Name([componentname, ':CavityDiag']);
            transform.Vector(0, '-dy', 0);
            transform.Repetitions(1);
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Transform('Shape', 'Translate');
        end
    end
end