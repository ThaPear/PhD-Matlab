function BuildCST(this, project)
    brick = project.Brick();
    transform = project.Transform();
    solid = project.Solid();

    project.StoreParameter('wcavity', this.wcavity*1e3);
    project.StoreParameter('wcavitydiag', this.wcavitydiag*1e3);
    project.MakeSureParameterExists('dms', 0.254);

    % x-oriented cavity
    brick.Reset();
    brick.Component('Cavity');
    brick.Name('Cavity');
    brick.Xrange('-dx', 'dx');
    brick.Yrange('slot_s0*dy-wcavity/2', 'slot_s0*dy+wcavity/2');
    brick.Zrange('-hback', '-dms');
    brick.Material('PEC');
    brick.Create();

    % Rotate by 90 degrees 4 times to make the entire cut.
    transform.Reset();
    transform.Name('Cavity:Cavity');
    transform.Angle(0, 0, 90);
    transform.Origin('Free');
    transform.Center('(slot_s0-0.5)*dx', '(slot_s0-0.5)*dy', 0);
    transform.Repetitions(3);
    transform.MultipleObjects(1);
    transform.GroupObjects(1);
    transform.Transform('Shape', 'Rotate');

    % Diagonal cavity
    brick.Reset();
    brick.Component('Cavity');
    brick.Name('CavityDiag');
    brick.Xrange('slot_s0*dx-wcavitydiag/2', 'slot_s0*dx+wcavitydiag/2');
    brick.Yrange('slot_s0*dy-wcavitydiag/2', 'slot_s0*dy+wcavitydiag/2');
    brick.Zrange('-hback', '-dms');
    brick.Material('PEC');
    brick.Create();

    % Rotate it 45 degrees
    transform.Reset();
    transform.Name('Cavity:CavityDiag');
    transform.Angle(0, 0, 45);
    transform.Origin('ShapeCenter');
    transform.Repetitions(1);
    transform.Transform('Shape', 'Rotate');

    % Rotate by 90 degrees 4 times to make the entire cut.
    transform.Reset();
    transform.Name('Cavity:CavityDiag');
    transform.Angle(0, 0, 90);
    transform.Origin('Free');
    transform.Center('(slot_s0-0.5)*dx', '(slot_s0-0.5)*dy', 0);
    transform.Repetitions(3);
    transform.MultipleObjects(1);
    transform.GroupObjects(1);
    transform.Transform('Shape', 'Rotate');

    % Add the diagonal cavity to the straight one.
    solid.Add('Cavity:Cavity', 'Cavity:CavityDiag');
end
