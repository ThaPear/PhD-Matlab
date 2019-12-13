clear;
SetupPath

desktop = HFSS.Application.Desktop();
hDesktop = desktop.hDesktop;
hProject = hDesktop.invoke('GetActiveProject');
if(isempty(hProject))
    hProject = hDesktop.invoke('NewProject');
end
% hProject.invoke('InsertDesign', "HFSS 3D Layout Design", "test", "","");

path = GetFullPath('testscript.vbs');

val = hDesktop.invoke('RunScript', path);
val
% hDefinitionManager = hProject.invoke('GetDefinitionManager');
% hDefinitionManager.invoke('AddMaterial', {"NAME:Material2", "dielectric_loss_tangent:=", "44"})