classdef Line_Trenched_Dualpol_Bowtie < Line
    properties
        % Parameters for the trench in CST.
        erreal
    end
    methods
        function this = Line_Trenched_Dualpol_Bowtie(er, L, erreal)
            this@Line(L, er);
            
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
            project.MakeSureParameterExists('wms', 'lfeed/2');
            
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
                condition = 'wfeed > (2*trench2metal+2*drillradius) And wms < lfeed And wms < lfeed';
                NOTcondition = 'wfeed <= (2*trench2metal+2*drillradius) Or wms >= lfeed Or wms >= lfeed';
%                 polygoncode = ['alpha = wslot/(2*trench2metal + 2*drillradius)'];
                polygon3d.Reset();
%                 polygon3d.Version(10);
                polygon3d.Name('Polygon');
                polygon3d.Curve('TrenchCut');
                % Down
                polygon3d.Point('-wslot/2+trench2metal', '-wslot/2+trench2metal', 0);
                    polygon3d.Point('-wslot/2+trench2metal', '-dy/2+lbtout/2+trench2metal', 0);
                    polygon3d.Point('-wfeed/2+trench2metal', '-dy/2+lbt/2+trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      If ', condition, ' Then', newline];
                        polygon3d.Point('-wfeed/2+trench2metal', '-dy/2+wms/2+trench2metal', 0);
                        polygon3d.Point('wfeed/2-trench2metal', '-dy/2+wms/2+trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      End If', newline];
                    polygon3d.Point('wfeed/2-trench2metal', '-dy/2+lbt/2+trench2metal', 0);
                    polygon3d.Point('wslot/2-trench2metal', '-dy/2+lbtout/2+trench2metal', 0);
                % Right
                polygon3d.Point('wslot/2-trench2metal', '-wslot/2+trench2metal', 0);
                    polygon3d.Point('dx/2-lbtout/2-trench2metal', '-wslot/2+trench2metal', 0);
                    polygon3d.Point('dx/2-lbt/2-trench2metal', '-wfeed/2+trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      If ', condition, ' Then', newline];
                        polygon3d.Point('dx/2-wms/2-trench2metal', '-wfeed/2+trench2metal', 0);
                        polygon3d.Point('dx/2-wms/2-trench2metal', 'wfeed/2-trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      End If', newline];
                    polygon3d.Point('dx/2-lbt/2-trench2metal', 'wfeed/2-trench2metal', 0);
                    polygon3d.Point('dx/2-lbtout/2-trench2metal', 'wslot/2-trench2metal', 0);
                % Up
                polygon3d.Point('wslot/2-trench2metal', 'wslot/2-trench2metal', 0);
                    polygon3d.Point('wslot/2-trench2metal', 'dy/2-lbtout/2-trench2metal', 0);
                    polygon3d.Point('wfeed/2-trench2metal', 'dy/2-lbt/2-trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      If ', condition, ' Then', newline];
                        polygon3d.Point('wfeed/2-trench2metal', 'dy/2-wms/2-trench2metal', 0);
                        polygon3d.Point('-wfeed/2+trench2metal', 'dy/2-wms/2-trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      End If', newline];
                    polygon3d.Point('-wfeed/2+trench2metal', 'dy/2-lbt/2-trench2metal', 0);
                    polygon3d.Point('-wslot/2+trench2metal', 'dy/2-lbtout/2-trench2metal', 0);
                % Left
                polygon3d.Point('-wslot/2+trench2metal', 'wslot/2-trench2metal', 0);
                    polygon3d.Point('-dx/2+lbtout/2+trench2metal', 'wslot/2-trench2metal', 0);
                    polygon3d.Point('-dx/2+lbt/2+trench2metal', 'wfeed/2-trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      If ', condition, ' Then', newline];
                        polygon3d.Point('-dx/2+wms/2+trench2metal', 'wfeed/2-trench2metal', 0);
                        polygon3d.Point('-dx/2+wms/2+trench2metal', '-wfeed/2+trench2metal', 0);
                    polygon3d.history = [polygon3d.history, '      End If', newline];
                    polygon3d.Point('-dx/2+lbt/2+trench2metal', '-wfeed/2+trench2metal', 0);
                    polygon3d.Point('-dx/2+lbtout/2+trench2metal', '-wslot/2+trench2metal', 0);
                polygon3d.Point('-wslot/2+trench2metal', '-wslot/2+trench2metal', 0);
                polygon3d.Create();
                    
                
                
                edges =    [[ 1,  2]; [ 1,  2]; [ 3, 28]; [ 2,  3]; [ 5,  6]; ...
                            [ 4,  5]; [ 7, 25]; [ 5,  6]; [ 9, 10]; [ 7,  8]; ...
                            [11, 22]; [ 8,  9]; [13, 14]; [10, 11]; [15, 19]; ...
                            [11, 12]];
                vertices = [[ 2,  2]; [ 2,  2]; [ 4,  4]; [ 4,  4]; [ 8,  8]; ...
                            [ 7,  7]; [10, 10]; [ 9,  9]; [14, 14]; [12, 12]; ...
                            [16, 16]; [14, 14]; [20, 20]; [17, 17]; [22, 22]; ...
                            [19, 19]];
                
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
                
%                 edges =    [[]; []; []; []; []; []; []; []; []; []; []; []; []; []; 
%                 vertices = [[]; []; []; []; []; []; []; []; []; []; []; []; []; []; 
%                 
%                 for(i = 1:length(edges))
%                     blendcurve.Reset();
%                     blendcurve.Name(['blend', num2str(i)]);
%                     blendcurve.Radius('drillradius');
%                     blendcurve.Curve('TrenchCut');
%                     blendcurve.CurveItem1('Polygon');
%                     blendcurve.CurveItem2('Polygon');
%                     blendcurve.EdgeId1(edges(i, 1));
%                     blendcurve.EdgeId2(edges(i, 2));
%                     blendcurve.VertexId1(vertices(i, 1));
%                     blendcurve.VertexId2(vertices(i, 2));
%                     blendcurve.CreateConditional(NOTcondition);
%                 end
                
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