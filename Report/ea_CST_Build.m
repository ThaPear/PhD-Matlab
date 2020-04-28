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

project.StartBulkMode();

solid = project.Solid();
cstplot = project.Plot();

% cstplot.DrawWorkplane(0);


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
project.StoreParameter('slot_impedance', zfeed);
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

project.StoreParameter('slot_width', wslot*1e3);
project.StoreParameter('slot_feedwidth', ['slot_width']);
project.StoreParameter('slot_feedlength', dslot*1e3);
project.StoreParameter('slot_bowtie_outer', ['dx-slot_width']);
project.StoreParameter('slot_s0', 0.25);
project.StoreParameter('adl_s0', adl_s0);

if(dualpol)
    project.StoreParameter('slot_feedwidth', wfeed*1e3);
	project.StoreParameter('adl_s0', -1);
end

% If it's single-pol, set shift to zero.
if(~dualpol)
	project.StoreParameter('slot_s0', 0);
end

dsproject = CST.Application.ActiveDS();
CST.BuildSchematic(dsproject, zfeed, C);
CST.AddSchematicTask(dsproject);

if(exist('hgap', 'var'))
    project.StoreParameter('hgap', hgap*1e3);
    for(i = 1:N)
        transform = project.Transform();
        transform.Reset();
        transform.Name(['ADS', num2str(i)]);
        transform.Vector(0, 0, 'hgap');
        transform.MultipleObjects(0);
        transform.GroupObjects(0);
        transform.Repetitions(1);
        transform.Transform('Shape', 'Translate');
    end
end

project.EndBulkMode();
project.Rebuild();
