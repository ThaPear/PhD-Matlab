classdef Cavity_Dualpol
    properties(SetAccess = protected)
        cavity_width
        cavity_width_diagonal
    end
    methods
        function this = Cavity_Dualpol(cavity_width, cavity_width_diagonal)
            this.cavity_width = cavity_width;
            this.cavity_width_diagonal = cavity_width_diagonal;
        end
        
        function BuildCST(this, project)
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
            brick.Component('Cavity');
            brick.Name('Cavity');
            brick.Xrange('-dx/2', 'dx/2');
            brick.Yrange('slot_s0*dy-cavity_width/2', 'slot_s0*dy+cavity_width/2');
            brick.Zrange('-hback', '-dms');
            brick.Material('PEC');
            brick.Create();
            
            brick.Reset();
            brick.Component('Cavity');
            brick.Name('Cavity2');
            brick.Xrange('slot_s0*dx-cavity_width/2', 'slot_s0*dx+cavity_width/2');
            brick.Yrange('-dy/2', 'dy/2');
            brick.Zrange('-hback', '-dms');
            brick.Material('PEC');
            brick.Create();
            
            solid.Add('Cavity:Cavity', 'Cavity:Cavity2');
            
            brick.Reset();
            brick.Component('Cavity');
            brick.Name('CavityDiag');
            brick.Xrange('slot_s0*dx-cavity_width_diagonal/2', 'slot_s0*dx+cavity_width_diagonal/2');
            brick.Yrange('slot_s0*dy-cavity_width_diagonal/2', 'slot_s0*dy+cavity_width_diagonal/2');
            brick.Zrange('-hback', '-dms');
            brick.Material('PEC');
            brick.Create();
        
            transform.Reset();
            transform.Name('Cavity:CavityDiag');
            transform.Angle(0, 0, 45);
            transform.Origin('ShapeCenter');
            transform.Repetitions(1);
            transform.Transform('Shape', 'Rotate');
        
            transform.Reset();
            transform.Name('Cavity:CavityDiag');
            transform.Vector('-dx', 0, 0);
            transform.Repetitions(1);
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Transform('Shape', 'Translate');
        
            transform.Reset();
            transform.Name('Cavity:CavityDiag');
            transform.Vector(0, '-dy', 0);
            transform.Repetitions(1);
            transform.MultipleObjects(1);
            transform.GroupObjects(1);
            transform.Transform('Shape', 'Translate');
        end
    end
end