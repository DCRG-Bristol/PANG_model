function [M_glob,F_glob] = buildAero(runObj,M_glob,F_glob, L, D, M, ang, q_all, p0, p)

%% runObj...

q = runObj.transF*q_all(runObj.structDisp);
qt = runObj.transF*q_all(runObj.structVel);

%%

vy = -cos(ang + runObj.geom.twist(runObj.basis.xColloc(:)));
vz = -sin(ang + runObj.geom.twist(runObj.basis.xColloc(:)));

vyBlk = diag(vy);
vzBlk = diag(vz);

%% get displacement functions...

proj_lift = analysis.aero.disps.liftDisp(q,qt,p0,p,vyBlk,vzBlk);
proj_moment = analysis.aero.disps.momentDisp(q,qt,p0,p,vyBlk,vzBlk);
proj_drag = analysis.aero.disps.dragDisp(q,qt,p0,p,vyBlk,vzBlk);

%% apply structural aerodynamic forces to structure...

F_glob(runObj.structVel) = F_glob(runObj.structVel) + runObj.transF'*(...
    proj_lift*L + proj_moment*M + proj_drag*D);

end

















