% Check if CST is desired.
ok = questdlg('Build CST project?','','Yes','No','Yes');
if(~strcmp(ok, 'Yes'))
    fprintf('Cancelled building of CST project.\n');
    return;
end

% Is the slot dualpol?
dualpol = ~strcmp(class(slot), 'Slot'); %#ok<STISA>

%% CST setup
project = CST.Util.InitializeProject();
solid = project.Solid();

dvia = 1.3e-3;
rvia = 0.35e-3;
wcavity = 1.65e-3;
wcavitydiag = 3.65e-3;

if(~dualpol)
    vias = Vias(dvia, rvia);
    cavity = Cavity(wcavity, wcavitydiag);
else
    vias = Vias_Dualpol(dvia, rvia);
    cavity = Cavity_Dualpol(wcavity, wcavitydiag);
end

project.StoreParameter('f0', f0/1e9);
project.StoreParameter('lambda0', 'c0/f0/1e9');
project.StoreParameter('slot_impedance', zfeed);
project.StoreParameter('dms', 0.127*2);

project.StoreParameter('slot_width', wslot*1e3);
project.StoreParameter('slot_feedwidth', ['slot_width']);
project.StoreParameter('slot_feedlength', dslot*1e3);
project.StoreParameter('slot_bowtie_outer', ['slot_feedlength']);
project.StoreParameter('slot_s0', 0.25);

slot.BuildCST(project);
vias.BuildCST(project);
cavity.BuildCST(project);
Materials.BuildCST(project);

% If it's single-pol, delete the walls and set shift to zero.
if(~dualpol)
    solid.Delete('Slot:Walls');
	project.StoreParameter('slot_s0', 0);
    project.Rebuild();
end

solid.Subtract('BackingReflector:Dielectric', 'Cavity:Cavity');
solid.Subtract('BackingReflector:Dielectric', 'Cavity:CavityDiag');