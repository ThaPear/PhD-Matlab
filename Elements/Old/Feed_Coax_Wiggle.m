classdef Feed_Coax
    properties
        
    end
    methods
        function this = Feed_Coax()
            
        end
        function BuildCST(this, project)
            cylinder = project.Cylinder();
            component = project.Component();
            brick = project.Brick();
            solid = project.Solid();
            port = project.Port();
            pick = project.Pick();
            discretefaceport = project.DiscreteFacePort();
            polygon3d = project.Polygon3D();
            blendcurve = project.BlendCurve();
            tracefromcurve = project.TraceFromCurve();
            
            lambda0 = Constants.c0/11e9;
            
            shield_ncyl = 5;
            shield_angles = linspace(0, 180, shield_ncyl);
            ratio = Coaxial.GetRadiusRatio(50, 3.66);
            
%             lms = 3.58e-3 - core_radius - core_ms_radius;

            project.StoreParameter('dms', '0.25');
            project.StoreParameter('lms', '2');
            project.StoreParameter('wms', 'dms * 28/25'); % 80 ohm
            project.StoreParameter('zms', 50);
            project.StoreParameter('res_w', 'dms * 39/25'); % 63 Ohm
            
            project.StoreParameter('core_height', 'hback');
            project.StoreParameter('core_radius', ['max(hback/12, 0.55/2)']);
            project.StoreParameter('core_ms_radius', 0.1);%'res_w/2 / sin(atn(res_w/2 / core_radius)) - core_radius + 0.001');
%             project.StoreParameter('core_ms_radius', 'max(0.001, wms/2 - core_radius)');
            project.StoreParameter('core_tophole_inner', 'core_radius+0.1');
            project.StoreParameter('core_tophole_outer', 'core_radius+0.3');

            
            project.StoreParameter('shield_radius', ['max(core_height/12, 0.25/2)']);
            project.StoreParameter('shield_height', 'hback');
            project.StoreParameter('shield_distance', [num2str(ratio), ' * core_radius + shield_radius + 0.1']);
            
            project.StoreParameter('shield_Nvias', '5');
            project.StoreParameter('shield_thetastart', 'pi/4');
            project.StoreParameter('shield_thetaend', 'pi + shield_thetastart');
            
            
            project.StoreParameter('res_ltotal', '0.25 * lambda0 / sqr(3.66) * 9.5/5.5');
            project.StoreParameter('res_l_one', '2');
            project.StoreParameter('res_r_one', '0.5');
            project.StoreParameter('res_l_two', '0.1');
            project.StoreParameter('res_r_two', '0.5');
            project.StoreParameter('res_l_three', 'res_ltotal - res_l_one - pi*res_r_one - res_l_two - pi*res_r_two');
            
            project.StoreParameter('ground_szx', 7);
            project.StoreParameter('ground_szy', 7);
            project.StoreParameter('lground', 'res_l_one - res_l_two + res_l_three + core_radius + lms');
            
            component.New('Feed');
            %% Ground planes & Substrate
            % Substrate
            brick.Reset();
            brick.Component('Feed');
            brick.Material('Generated/3.66');
            brick.Name('Dielectric');
            brick.Xrange('-ground_szx/2', 'ground_szx/2');
            brick.Yrange('-lground', 'ground_szy/2');
            brick.Zrange('-hback', '0');
            brick.Create();
            
            % Top ground plane
            brick.Reset();
            brick.Component('Feed');
            brick.Material('PEC');
            brick.Name('TopGround');
            brick.Xrange('-ground_szx/2', 'ground_szx/2');
            brick.Yrange('-lground', 'ground_szy/2');
            brick.Zrange('0', '0');
            brick.Create();
            
            % Bottom ground plane
            brick.Name('BottomGround');
            brick.Zrange('-hback', '-hback');
            brick.Create();
            
            cylinder.Reset();
            cylinder.Component('Feed');
            cylinder.Material('PEC');
            cylinder.Axis('z');
            project.AddToHistory('Set cylinder settings', cylinder.history)
            cylinder.history = '';
            
            %% Coaxial
            % Hole in bottom ground plane
            cylinder.Name('Hole');
            cylinder.Innerradius('core_radius');
            cylinder.Outerradius('shield_distance - shield_radius - 0.1');
            cylinder.Zrange('-hback', '-hback');
            cylinder.Xcenter(0);
            cylinder.Ycenter(0);
            cylinder.Zcenter(0);
            cylinder.Create();
            
            solid.Subtract('Feed:BottomGround', 'Feed:Hole');
            
            % Hole in top ground plane
            cylinder.Innerradius('core_tophole_inner');
            cylinder.Outerradius('core_tophole_outer');
            cylinder.Zrange('0', '0');
            cylinder.Create();
            
            solid.Subtract('Feed:TopGround', 'Feed:Hole');
            
            % Outer cylinders
            cylinder.Name(['Shield " & i & "']);
            cylinder.Innerradius(0);
            cylinder.Outerradius('shield_radius');
            cylinder.Zrange('-shield_height', '0');
            cylinder.Xcenter(['" & cos(shield_thetastart + (shield_thetaend - shield_thetastart) * i / (shield_Nvias-1)) * shield_distance & "']);
            cylinder.Ycenter(['" & sin(shield_thetastart + (shield_thetaend - shield_thetastart) * i / (shield_Nvias-1)) * shield_distance & "']);
            cylinder.Zcenter(0);
            cylinder.CreateForLoop('i', 1, 'shield_Nvias');
            
            % Core cylinder
            cylinder.Name('Core');
            cylinder.Innerradius(0);
            cylinder.Outerradius('core_radius');
            cylinder.Zrange('-hback', '-(hback-core_height)');
            cylinder.Xcenter(0);
            cylinder.Ycenter(0);
            cylinder.Zcenter(0);
            cylinder.Create();
            
            % Transfer from cylinder to microstrip
            cylinder.Name('Core_to_ms');
            cylinder.Innerradius('core_radius');
            cylinder.Outerradius('core_radius + core_ms_radius');
            cylinder.Zrange('-dms', '-dms');
            cylinder.Create();
            
            % Coax port
            port.Reset();
            port.NumberOfModes(1);
            port.Xrange('-shield_distance', 'shield_distance');
            port.Yrange('-shield_distance', 'shield_distance');
            port.Zrange('-hback', '-hback');
            port.PortOnBound(0);
            port.Orientation('zmin');
            port.PortNumber(2);
            port.Create();
            
            
            
            %% Microstrips
            % Quarter-wave resonator
            polygon3d.Reset();
            polygon3d.Name('FeedResonator');
            polygon3d.Curve('FeedResonator');
%             polygon3d.Point(x,y,z);
            polygon3d.Point(0                               , '-core_radius'                                        , '-dms');
            polygon3d.Point(0                               , '-core_radius - res_l_one - res_r_one'                , '-dms');
            polygon3d.Point('2 * res_r_one'                 , '-core_radius - res_l_one - res_r_one'                , '-dms');
            polygon3d.Point('2 * res_r_one'                 , '-core_radius - res_l_one + res_l_two + res_r_two'    , '-dms');
            polygon3d.Point('2 * res_r_one + 2 * res_r_two' , '-core_radius - res_l_one + res_l_two + res_r_two'    , '-dms');
            polygon3d.Point('2 * res_r_one + 2 * res_r_two' , '-core_radius - res_l_one + res_l_two - res_l_three'  , '-dms');
            polygon3d.Create();
            
            % Blend curves of microstrip
            blendcurve.Reset();
            blendcurve.Name('FeedResonator_Corner1');
            blendcurve.Radius('res_r_one');
            blendcurve.Curve('FeedResonator');
            blendcurve.CurveItem1('FeedResonator');
            blendcurve.CurveItem2('FeedResonator');
            blendcurve.EdgeId1(1);
            blendcurve.EdgeId2(2);
            blendcurve.VertexId1(2);
            blendcurve.VertexId2(2);
            blendcurve.CreateConditional('res_r_one > 0');
            
            blendcurve.Name('FeedResonator_Corner2');
            blendcurve.EdgeId1(2);
            blendcurve.EdgeId2(5);
            blendcurve.VertexId1(3);
            blendcurve.VertexId2(3);
            blendcurve.CreateConditional('res_r_one > 0');
            
            blendcurve.Name('FeedResonator_Corner3');
            blendcurve.Radius('res_r_two');
            blendcurve.EdgeId1(1);
            blendcurve.EdgeId2(3);
            blendcurve.VertexId1(1);
            blendcurve.VertexId2(1);
            blendcurve.CreateConditional('res_r_one > 0');
            
            blendcurve.Name('FeedResonator_Corner4');
            blendcurve.EdgeId1(3);
            blendcurve.EdgeId2(4);
            blendcurve.VertexId1(6);
            blendcurve.VertexId2(6);
            blendcurve.CreateConditional('res_r_one > 0');
            
            % Change curve to microstrip
            tracefromcurve.Reset();
            tracefromcurve.Name('Resonator')
            tracefromcurve.Component('Feed');
            tracefromcurve.Material('PEC');
            tracefromcurve.Curve('FeedResonator');
%             tracefromcurve.Thickness(0);
            tracefromcurve.Width('res_w');
%             tracefromcurve.RoundStart(0);
%             tracefromcurve.RoundEnd(0);
%             tracefromcurve.DeleteCurve(1);
%             tracefromcurve.GapType(2);
            tracefromcurve.Create();
            
            solid.Add('Feed:Resonator', 'Feed:Core_to_ms');
            
            brick.Reset();
            brick.Component('Feed');
            brick.Material('PEC');
            brick.Name('Microstrip');
            brick.Xrange('2 * res_r_one + 2 * res_r_two - wms/2', '2 * res_r_one + 2 * res_r_two + wms/2');
            brick.Yrange('-core_radius - res_l_one + res_l_two - res_l_three', '-core_radius - res_l_one + res_l_two - res_l_three - lms');
            brick.Zrange('-dms', '-dms');
            brick.Create();
            
            % Microstrip port
            port.Reset();
            port.NumberOfModes(1);
            port.Xrange('2 * res_r_one + 2 * res_r_two - wms/2 * 4', '2 * res_r_one + 2 * res_r_two + wms/2 * 4');
            port.Yrange('-core_radius - res_l_one + res_l_two - res_l_three - lms', '-core_radius - res_l_one + res_l_two - res_l_three - lms');
            port.Zrange('-dms * 2', '0');
            port.PortOnBound(0);
            port.Orientation('ymin');
            port.PortNumber(1);
            port.Create();
            
%             % Microstrip port
%             brick.Name('FeedPort');
%             brick.Material('Vacuum');
%             brick.Xrange('-wms/2', 'wms/2');
%             brick.Yrange('- lms', '- lms');
%             brick.Zrange('-dms', '0');
%             brick.Create();
% 
%             % Define the port.
%             pick.PickEdgeFromId('Feed:FeedPort', 2, 2);
%             pick.PickEdgeFromId('Feed:FeedPort', 4, 4);
%             discretefaceport.Reset();
%             discretefaceport.PortNumber(1);
%             discretefaceport.Label('SlotFeedPort');
%             discretefaceport.Impedance('zms');
%             discretefaceport.SetP1(1);
%             discretefaceport.SetP2(1);
%             discretefaceport.CenterEdge(1);
%             discretefaceport.Monitor(1);
%             discretefaceport.Create();
        end
    end
end
    