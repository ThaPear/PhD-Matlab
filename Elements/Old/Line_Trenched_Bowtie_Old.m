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
            project.MakeSureParameterExists('bowtieratio', 'wfeed/wslot');
            
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
                
                %% Square
                % Cut the trench.
                rectangle.Reset();
                rectangle.Name('Rectangle');
                rectangle.Curve('TrenchCut');
                rectangle.Xrange('-dx/2', '-lfeed/2-trench2metal');
                rectangle.Yrange('-wslot/2+trench2metal', 'wslot/2-trench2metal');
                rectcode = ['With Rectangle', newline, rectangle.history, '     .Create', newline, 'End With'];
                
                blendcurve.Reset();
                blendcurve.Name('blend1');
                blendcurve.Radius('drillradius');
                blendcurve.Curve('TrenchCut');
                blendcurve.CurveItem1('Rectangle');
                blendcurve.CurveItem2('Rectangle');
                blendcurve.EdgeId1(2);
                blendcurve.EdgeId2(3);
                blendcurve.VertexId1(3);
                blendcurve.VertexId2(3);
                rectcode = [rectcode, newline, 'With BlendCurve', newline, blendcurve.history, '     .Create', newline, 'End With'];
                blendcurve.history = strrep(blendcurve.history, 'blend1', 'blend2');
                rectcode = [rectcode, newline, 'With BlendCurve', newline, blendcurve.history, '     .Create', newline, 'End With'];
                
                rectcode = ['If (bowtieratio <= 1) And (wslot > 2*trench2metal+2*drillradius+0.1) Then', newline, rectcode, newline, 'End If'];
                project.AddToHistory('Conditional Rectangle', rectcode);
                
                %% House
%                 polygoncode = ['alpha = wslot/(2*trench2metal + 2*drillradius)'];
                polygon3d.Reset();
%                 polygon3d.Version(10);
                polygon3d.Name('Polygon');
                polygon3d.Curve('TrenchCut');
                polygon3d.Point('-dx/2', '-wslot/2+trench2metal', '0');
                polygon3d.Point('-lfeed/2-trench2metal', '-wslot/2+trench2metal', '0');
                polygon3d.Point('-lfeed/2/Min(bowtieratio, wslot/(2*trench2metal + 2*drillradius))-trench2metal', ...
                                '-wslot/2/Min(bowtieratio, wslot/(2*trench2metal + 2*drillradius))+trench2metal', ...
                                '0');
                polygon3d.Point('-lfeed/2/Min(bowtieratio, wslot/(2*trench2metal + 2*drillradius))-trench2metal',  ...
                                'wslot/2/Min(bowtieratio, wslot/(2*trench2metal + 2*drillradius))-trench2metal', ...
                                '0');
                polygon3d.Point('-lfeed/2-trench2metal', 'wslot/2-trench2metal', '0');
                polygon3d.Point('-dx/2', 'wslot/2-trench2metal', '0');
                polygon3d.Point('-dx/2', '-wslot/2+trench2metal', '0');
                polygoncode = ['With Polygon3D', newline, polygon3d.history, '     .Create', newline, 'End With'];
                
                blendcurve.Reset();
                blendcurve.Name('blend1');
                blendcurve.Radius('drillradius');
                blendcurve.Curve('TrenchCut');
                blendcurve.CurveItem1('Polygon');
                blendcurve.CurveItem2('Polygon');
                blendcurve.EdgeId1(3);
                blendcurve.EdgeId2(4);
                blendcurve.VertexId1(4);
                blendcurve.VertexId2(4);
                
                polygoncode = [polygoncode, newline, 'With BlendCurve', newline, blendcurve.history, '     .Create', newline, 'End With'];
                blendcurve.history = strrep(blendcurve.history, 'blend1', 'blend2');
                polygoncode = [polygoncode, newline, 'With BlendCurve', newline, blendcurve.history, '     .Create', newline, 'End With'];
                % Cut off the edge & vertex IDs.
                blendcurve.history = blendcurve.history(1:end-76);
                blendcurve.EdgeId1(5);
                blendcurve.EdgeId2(6);
                blendcurve.VertexId1(7);
                blendcurve.VertexId2(7);
                blendcurve.history = strrep(blendcurve.history, 'blend2', 'blend3');
                polygoncode = [polygoncode, newline, 'With BlendCurve', newline, blendcurve.history, '     .Create', newline, 'End With'];
                % Cut off the edge & vertex IDs.
                blendcurve.history = blendcurve.history(1:end-76);
                blendcurve.EdgeId1(5);
                blendcurve.EdgeId2(6);
                blendcurve.VertexId1(8);
                blendcurve.VertexId2(8);
                blendcurve.history = strrep(blendcurve.history, 'blend3', 'blend4');
                polygoncode = [polygoncode, newline, 'With BlendCurve', newline, blendcurve.history, '     .Create', newline, 'End With'];
                
                polygoncode = ['If (bowtieratio > 1) And (wslot > 2*trench2metal+2*drillradius+0.1) Then', newline, polygoncode, newline, 'End If'];
                project.AddToHistory('Conditional Polygon', polygoncode);
                
                %% Make into brick.
                extrudecurve.Reset();
                extrudecurve.Curve('TrenchCut');
                extrudecurve.Component('TrenchedLines');
                extrudecurve.Name('TrenchCut');
                extrudecurve.Material('PEC');
                extrudecurve.Thickness('htrench');
                extrudecurve.DeleteProfile(1);
                extrudecode = ['With ExtrudeCurve', newline, extrudecurve.history, '     .Create', newline, 'End With'];
                
                %% Mirror in X-axis.
                transform.Reset();
                transform.Name('TrenchedLines:TrenchCut');
                transform.PlaneNormal(1, 0, 0);
                transform.Center(0, 0, 0);
                transform.Origin('Free');
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                extrudecode = [extrudecode, newline, 'With Transform', newline, transform.history, '     .Transform "Shape", "Mirror"', newline, 'End With'];
                extrudecode = ['If (wslot > 2*trench2metal+2*drillradius+0.1) Then', newline, extrudecode, newline, 'End If'];
                project.AddToHistory('Conditional Extrude', extrudecode);
                
% Old code for cutting the trench
%{
                brick.Name('TrenchCut');
                brick.Xrange('-dx/2', '-lfeed/2-trench2metal-drillradius');
                brick.Yrange('-wslot/2+trench2metal', 'wslot/2-trench2metal');
                brick.Zrange(0, obj.L*1e3);
                brick.Create();
                
                brick.Name('TrenchPart2');
                brick.Xrange('-dx/2', '-lfeed/2-trench2metal');
                brick.Yrange('-wslot/2+trench2metal+drillradius', 'wslot/2-trench2metal-drillradius');
                brick.Zrange(0, obj.L*1e3);
                brick.Create();
                
                cylinder.Reset();
                cylinder.Name('TrenchPart3');
                cylinder.Component('TrenchedLines');
                cylinder.Axis('z');
                cylinder.Xcenter('-lfeed/2-trench2metal-drillradius');
                cylinder.Ycenter('wslot/2-trench2metal-drillradius');
                cylinder.Zrange(0, obj.L*1e3);
                cylinder.Outerradius('drillradius');
                cylinder.Create();
                
                cylinder.Name('TrenchPart4');
                cylinder.Ycenter('-wslot/2+trench2metal+drillradius');
                cylinder.Create();
                
                solid.Add('TrenchedLines:TrenchCut', 'TrenchedLines:TrenchPart2')
                solid.Add('TrenchedLines:TrenchCut', 'TrenchedLines:TrenchPart3')
                solid.Add('TrenchedLines:TrenchCut', 'TrenchedLines:TrenchPart4')
                
                transform.Reset();
                transform.Name('TrenchedLines:TrenchPart1');
                transform.PlaneNormal(1, 0, 0);
                transform.Center(0, 0, 0);
                transform.Origin('Free');
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Transform('Shape', 'Mirror');
%}
                subtractcode = ['Solid.Subtract "', ['TrenchedLines:', brickname], '", "', 'TrenchedLines:TrenchCut', '"'];
                subtractcode = ['If (wslot > 2*trench2metal+2*drillradius+0.1) Then', newline, subtractcode, newline, 'End If'];
                project.AddToHistory('Conditional Subtract', subtractcode);
                
                solid.MergeMaterialsOfComponent(['Lines:', brickname]);
            end
            
            wcs = project.WCS();
            wcs.MoveWCS('local', 0, 0, 'htrench');
        end
    end
end