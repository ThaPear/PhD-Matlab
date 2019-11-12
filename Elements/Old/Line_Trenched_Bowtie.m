classdef Line_Trenched_Bowtie < Line
    properties
        % Parameters for the trench in CST.
        erreal
    end
    methods
        function this = Line_Trenched_Bowtie(er, L, erreal)
            this@Line(er, L);
            
            this.erreal = erreal;
        end
        function BuildCST(this, project)
            material = project.Material();
            brick = project.Brick();
            solid = project.Solid();
%             cylinder = project.Cylinder();
            transform = project.Transform();
            extrudecurve = project.ExtrudeCurve();
            rectangle = project.Rectangle();
            blendcurve = project.BlendCurve();
            polygon3d = project.Polygon3D();
            project.MakeSureParameterExists('drillradius', Limits.drillradius*1e3);
            project.MakeSureParameterExists('trench2metal', Limits.trench2metal*1e3);
            project.MakeSureParameterExists('hback', this.L*1e3);
            project.MakeSureParameterExists('htrench', 'hback');
            
            if(this.erreal ~= 1)
                % Create necessary material
                material.Reset();
                materialname = num2str(this.erreal, 5);
                material.Name(materialname);
                material.Folder('Generated');
                material.Colour(0, min(1, this.erreal/20), 1);
                material.Epsilon(this.erreal);
                material.Transparency(0.5);
                material.Create();


                % Create the dielectric.
                brick.Reset();
                brick.Component('TrenchedLines');
                brickname = ['Line ', num2str(floor(rand()*1e4))];
                brick.Name(brickname);
                brick.Xrange('-dx/2', 'dx/2');
                brick.Yrange('-dy/2', 'dy/2');
                brick.Zrange(0, this.L*1e3);
                brick.Material(['Generated/', materialname]);
                brick.Create();
                
                % Draw the trench shape.
                polygon3d.Reset();
                polygon3d.Curve('TrenchCut');
                polygon3d.Name('TrenchCut');
                polygon3d.Point('dx/2'   , '-wslot/2+trench2metal', '0');
                polygon3d.Point('lbtout/2+trench2metal', '-wslot/2+trench2metal', '0');
                polygon3d.Point('lbt/2+trench2metal', '-wfeed/2+trench2metal', '0');
                polygon3d.Point('lbt/2+trench2metal', 'wfeed/2-trench2metal', '0');
                polygon3d.Point('lbtout/2+trench2metal', 'wslot/2-trench2metal', '0');
                polygon3d.Point('dx/2'   , 'wslot/2-trench2metal', '0');
                polygon3d.Point('dx/2'   , '-wslot/2+trench2metal', '0');
                polygon3d.Create();
                
                % Extrude the trench shape.
                extrudecurve.Reset();
                extrudecurve.Curve('TrenchCut');
                extrudecurve.Component('TrenchedLines');
                extrudecurve.Name('TrenchCut');
                extrudecurve.Material('PEC');
                extrudecurve.Thickness('-htrench');
                extrudecurve.DeleteProfile(1);
                extrudecurve.Create();
                
                % Mirror in x=0.
                transform.Reset();
                transform.Name('TrenchedLines:TrenchCut');
                transform.PlaneNormal(1, 0, 0);
                transform.Center(0, 0, 0);
                transform.Origin('Free');
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Transform('Shape', 'Mirror');
                
                % Subtract from dielectric.
                solid.Subtract(['TrenchedLines:', brickname], 'TrenchedLines:TrenchCut');
            
%                 solid.MergeMaterialsOfComponent(['Lines:', brickname]);
            end
            
            wcs = project.WCS();
            wcs.MoveWCS('local', 0, 0, 'htrench');
        end
    end
end