% Check if CST is desired.
% ok = questdlg('Build CST project?','','Yes','No','Yes');
% if(~strcmp(ok, 'Yes'))
%     dispex('Cancelled building of CST project.\n');
%     return;
% end

% Is the slot dualpol?
dualpol = contains(lower(class(slot)), 'dualpol');

%% CST setup
project = CST.InitializePeriodicProject();
solid = project.Solid();

dvia = 1.3e-3;
rvia = 0.35e-3;
wcavity = 1.65e-3;
wcavitydiag = 3.65e-3;
wcavity = 0e-3;
wcavitydiag = 3e-3;

adl_s0 = 0;

if(dualpol)
    viasobj = Vias_Dualpol(dvia, rvia);
    cavityobj = Cavity_Dualpol(wcavity, wcavitydiag);
else
    viasobj = Vias(dvia, rvia);
    cavityobj = Cavity(wcavity, wcavitydiag);
end

project.StoreParameter('f0', f0/1e9);
project.StoreParameter('lambda0', 'c0/f0/1e9');
project.StoreParameter('slot_impedance', params.zfeed);
project.StoreParameter('dms', 0.127);

slot.BuildCST(project);
if(vias)
    viasobj.BuildCST(project);
end
if(cavity)
    cavityobj.BuildCST(project);
    solid.Subtract('BackingReflector:Dielectric', 'Cavity:Cavity');
    solid.Subtract('BackingReflector:Dielectric', 'Cavity:CavityDiag');
end
Materials.BuildCST(project);

project.StoreParameter('slot_width', params.wslot*1e3);
project.StoreParameter('slot_feedwidth', ['slot_width']);
project.StoreParameter('slot_feedlength', params.dslot*1e3);
project.StoreParameter('slot_bowtie_outer', ['slot_feedlength']);
project.StoreParameter('slot_s0', 0.25);
project.StoreParameter('adl_s0', adl_s0);

% If it's single-pol, delete the walls and set shift to zero.
if(~dualpol)
	project.StoreParameter('slot_s0', 0);
    project.Rebuild();
end

dsproject = CST.Application.ActiveDS();
CST.BuildSchematic(dsproject, params.zfeed, C);
CST.AddSchematicTask(dsproject);