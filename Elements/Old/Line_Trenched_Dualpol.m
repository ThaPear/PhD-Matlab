classdef Line_Trenched_Dualpol < Line
    properties
        % Parameters for the trench in CST.
        erreal
    end
    methods
        function this = Line_Trenched_Dualpol(er, L, erreal)
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
%             rectangle = project.Rectangle();
            blendcurve = project.BlendCurve();
            polygon3d = project.Polygon3D();
            wcs = project.WCS();
            project.MakeSureParameterExists('drillradius', Limits.drillradius*1e3);
            project.MakeSureParameterExists('trench2metal', Limits.trench2metal*1e3);
            project.MakeSureParameterExists('htrench', this.L*1e3);
            project.MakeSureParameterExists('wms', 'lfeed');
            
            if(this.erreal ~= 1)
                
                if(Globals.exists('slot_s0'))
                    s0 = Globals.slot_s0;
                else
                    s0 = 0;
                end
                
                project.MakeSureParameterExists('slot_s0', s0);
                
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
                brick.Zrange(0, 'htrench');
                brick.Material(['Generated/', materialname]);
                brick.Create();
                
                wcs.Enable();
                wcs.Store('Pre-Trench');
                wcs.MoveWCS('local', 'slot_s0 * dx', 'slot_s0 * -dy', 0);
                
                %% Polygon
                condition = 'wfeed > (2*trench2metal+2*drillradius) And wms < lfeed';
                NOTcondition = 'wfeed <= (2*trench2metal+2*drillradius) Or wms >= lfeed';
%                 polygoncode = ['alpha = wslot/(2*trench2metal + 2*drillradius)'];
                polygon3d.Reset();
%                 polygon3d.Version(10);
                polygon3d.Name('Polygon');
                polygon3d.Curve('TrenchCut');
                polygon3d.Point('-wslot/2+trench2metal', '-wslot/2+trench2metal', 0);
                    polygon3d.Point('-wfeed/2+trench2metal', '-dy/2+lfeed/2+trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      If ', condition, ' Then', newline];
                        polygon3d.Point('-wfeed/2+trench2metal', '-dy/2+wms/2+trench2metal', 0);
                        polygon3d.Point(' wfeed/2-trench2metal', '-dy/2+wms/2+trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      End If', newline];
                    polygon3d.Point(' wfeed/2-trench2metal', '-dy/2+lfeed/2+trench2metal', 0);
                polygon3d.Point(' wslot/2-trench2metal', '-wslot/2+trench2metal', 0);
                    polygon3d.Point('dx/2-lfeed/2-trench2metal', '-wfeed/2+trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      If ', condition, ' Then', newline];
                        polygon3d.Point('dx/2-wms/2-trench2metal', '-wfeed/2+trench2metal', 0);
                        polygon3d.Point('dx/2-wms/2-trench2metal', ' wfeed/2-trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      End If', newline];
                    polygon3d.Point('dx/2-lfeed/2-trench2metal', 'wfeed/2-trench2metal', 0);
                polygon3d.Point(' wslot/2-trench2metal', 'wslot/2-trench2metal', 0);
                    polygon3d.Point( 'wfeed/2-trench2metal', 'dy/2-lfeed/2-trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      If ', condition, ' Then', newline];
                        polygon3d.Point(' wfeed/2-trench2metal', 'dy/2-wms/2-trench2metal', 0);
                        polygon3d.Point('-wfeed/2+trench2metal', 'dy/2-wms/2-trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      End If', newline];
                    polygon3d.Point('-wfeed/2+trench2metal', 'dy/2-lfeed/2-trench2metal', 0);
                polygon3d.Point('-wslot/2+trench2metal', 'wslot/2-trench2metal', 0);
                    polygon3d.Point('-dx/2+lfeed/2+trench2metal', 'wfeed/2-trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      If ', condition, ' Then', newline];
                        polygon3d.Point('-dx/2+wms/2+trench2metal', ' wfeed/2-trench2metal', 0);
                        polygon3d.Point('-dx/2+wms/2+trench2metal', '-wfeed/2+trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      End If', newline];
                    polygon3d.Point('-dx/2+lfeed/2+trench2metal', '-wfeed/2+trench2metal', 0);
                polygon3d.Point('-wslot/2+trench2metal', '-wslot/2+trench2metal', 0);
                polygon3d.Create();
                
                edges =    [[ 2, 3]; [ 1,20]; [ 4, 5]; [ 5,19]; [ 5, 6]; [ 9,18]; [ 6, 7]; [13,17]];
                vertices = [[ 3, 3]; [ 1, 1]; [ 6, 6]; [ 6, 6]; [ 8, 8]; [11,11]; [10,10]; [16,16]];
                
                for(i = 1:length(edges))
                    blendcurve.Reset();
                    blendcurve.Name(['blend', num2str(i)]);
                    blendcurve.Radius('drillradius');
                    blendcurve.Curve('TrenchCut');
                    blendcurve.CurveItem1('Polygon');
                    blendcurve.CurveItem2('Polygon');
                    blendcurve.EdgeId1(edges(i, 1));
                    blendcurve.EdgeId2(edges(i, 2));
                    blendcurve.VertexId1(vertices(i, 1));
                    blendcurve.VertexId2(vertices(i, 2));
                    blendcurve.CreateConditional(condition);
                end
                
                edges =    [[ 1, 2]; [ 1,12]; [ 2, 3]; [ 3,11]; [ 3, 4]; [ 5,10]; [ 4, 5]; [ 7, 9]];
                vertices = [[ 2, 2]; [ 1, 1]; [ 4, 4]; [ 4, 4]; [ 6, 6]; [ 7, 7]; [ 8, 8]; [10,10]];
                
                for(i = 1:length(edges))
                    blendcurve.Reset();
                    blendcurve.Name(['blend', num2str(i)]);
                    blendcurve.Radius('drillradius');
                    blendcurve.Curve('TrenchCut');
                    blendcurve.CurveItem1('Polygon');
                    blendcurve.CurveItem2('Polygon');
                    blendcurve.EdgeId1(edges(i, 1));
                    blendcurve.EdgeId2(edges(i, 2));
                    blendcurve.VertexId1(vertices(i, 1));
                    blendcurve.VertexId2(vertices(i, 2));
                    blendcurve.CreateConditional(NOTcondition);
                end
                
                extrudecurve.Reset();
                extrudecurve.Curve('TrenchCut');
                extrudecurve.Component('TrenchedLines');
                extrudecurve.Name('TrenchCut');
                extrudecurve.Material('PEC');
                extrudecurve.Thickness('htrench');
                extrudecurve.DeleteProfile(1);
                extrudecurve.Create();
                
                % Mirror the trench in x = 0.
                transform.Reset();
                transform.Name('TrenchedLines:TrenchCut');
                transform.Vector('0', 'dy', '0');
                transform.MultipleObjects(1);
                transform.GroupObjects(1);
                transform.Transform('Shape', 'Translate');

                transform.Vector('-dx', '0', '0');
                transform.Transform('Shape', 'Translate');
                
                solid.Subtract(['TrenchedLines:', brickname], 'TrenchedLines:TrenchCut');
                
                wcs.Restore('Pre-Trench');
            end
            
            wcs = project.WCS();
            wcs.MoveWCS('local', 0, 0, 'htrench');
        end
    end
end