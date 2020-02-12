
project = CST.Application.Active3D();
plot = project.Plot();
plot.RestoreView('Bottom');
plot.RotationAngle(360-125); plot.Rotate('Right');
plot.RotationAngle(20); plot.Rotate('Down');
plot.Update();
plot.ResetZoom();
% 6 ticks zoom in