function [M_glob,F_glob] = buildStruct(runObj,M_glob,F_glob,q_all,p0,p)

% assign variables...

q = runObj.transF*q_all(runObj.structDisp);
qt = runObj.transF*q_all(runObj.structVel);

if nargin<5
    p = [];
end

Mstr = analysis.inertia.massMatr(q,qt, p0, p);
Fstr = analysis.elas.elasForces(q,qt, p0, p) + analysis.grav.gravForces(q,qt, p0, p) +...
    -runObj.dampMatr*qt;
% 
% Mstr = runObj.massMatr(q,qt, p0, p);
% Fstr = runObj.elasForces(q,qt, p0, p) + runObj.gravForces(q,qt, p0, p) +...
%     -runObj.dampMatr*qt;

M_glob(runObj.structDisp,runObj.structDisp) = eye(runObj.Nstr);
M_glob(runObj.structVel,runObj.structVel) = runObj.transF'*Mstr*runObj.transF;

F_glob(runObj.structDisp) = q_all(runObj.Nstr+1:2*runObj.Nstr);
F_glob(runObj.structVel) = F_glob(runObj.structVel) + runObj.transF'*Fstr;

end

















