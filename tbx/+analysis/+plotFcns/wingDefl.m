function [u, v, w] = wingDefl(obj, q_all, xi)

x_cont = xi;
b = obj.geom.b;
a = obj.geom.a;

for j=1:length(x_cont)
    uv{j} = project.basis.fv(x_cont(j)); uw{j} = project.basis.fw(x_cont(j));
end

q = obj.transF*q_all(obj.structDisp);

% 
% u = zeros(2, length(xi));
% w = zeros(2, length(xi));
% v = zeros(2, length(xi));

rLE = zeros(3, length(xi));
rTE = zeros(3, length(xi));

for j=1:length(x_cont)
    x = x_cont(j);

    bLE = b(x)*(1+a(x));
    bTE = -b(x)*(1-a(x));

    dr = [1-0.5*(q'*(project.basis.W1(x).*project.basis.W1(x)')*q + q'*(project.basis.V1(x).*project.basis.V1(x)')*q);...
        q'*project.basis.V1(x);...
        q'*project.basis.W1(x)];

    ey = [-project.basis.V1(x)'*q*(1-0.5*(project.basis.T(x)'*q)^2 + 0.5*(project.basis.W1(x)'*q)^2) +...
        -(project.basis.W1(x)'*q)*(project.basis.T(x)'*q);...
        1-0.5*(project.basis.T(x)'*q)^2+...
        -0.5*(project.basis.V1(x)'*q)^2 - (project.basis.V1(x)'*q)*(project.basis.T(x)'*q)*(project.basis.W1(x)'*q);...
        (project.basis.T(x)'*q)*(1-0.5*(project.basis.W1(x)'*q)^2) +...
        -(1/6)*(project.basis.T(x)'*q)^3];

    rEA = [x - 0.5*(q'*uw{j}*q + q'*uv{j}*q ); project.basis.V(x)'*q; project.basis.W(x)'*q];
    
    ez = cross(dr,ey);

    yLE = bLE*cos(obj.geom.twist(x));
    zLE = bLE*sin(obj.geom.twist(x));
    yTE = bTE*cos(obj.geom.twist(x));
    zTE = bTE*sin(obj.geom.twist(x));

    rLE(:,j) = rEA + yLE*ey + zLE*ez;
    rTE(:,j) = rEA + yTE*ey + zTE*ez;

    %%
%     w(1,j) = project.basis.W(x)'*q + bLE*((project.basis.T(x)'*q)*(1-0.5*(project.basis.W1(x)'*q)^2) +...
%         -(1/6)*(project.basis.T(x)'*q)^3 );
% 
%     v(1,j) = project.basis.V(x)'*q + bLE*(1-0.5*(project.basis.T(x)'*q)^2+...
%         -0.5*(project.basis.V1(x)'*q)^2 - (project.basis.V1(x)'*q)*(project.basis.T(x)'*q)*(project.basis.W1(x)'*q));
% 
%     u(1,j) = x - 0.5*(q'*uw{j}*q + q'*uv{j}*q ) +...
%         bLE*(-project.basis.V1(x)'*q*(1-0.5*(project.basis.T(x)'*q)^2 + 0.5*(project.basis.W1(x)'*q)^2) +...
%         -(project.basis.W1(x)'*q)*(project.basis.T(x)'*q));
% 
%     w(2,j) = project.basis.W(x)'*q - bTW*((project.basis.T(x)'*q)*(1-0.5*(project.basis.W1(x)'*q)^2) +...
%         -(1/6)*(project.basis.T(x)'*q)^3);
% 
%     v(2,j) = project.basis.V(x)'*q - bTW*(1-0.5*(project.basis.T(x)'*q)^2+...
%         -0.5*(project.basis.V1(x)'*q)^2 - (project.basis.V1(x)'*q)*(project.basis.T(x)'*q)*(project.basis.W1(x)'*q));
% 
%     u(2,j) = x - 0.5*(q'*uw{j}*q + q'*uv{j}*q ) +...
%         -bTW*(-project.basis.V1(x)'*q*(1-0.5*(project.basis.T(x)'*q)^2 + 0.5*(project.basis.W1(x)'*q)^2) +...
%         -(project.basis.W1(x)'*q)*(project.basis.T(x)'*q));

end

u = [rLE(1,:); rTE(1,:)];
v = [rLE(2,:); rTE(2,:)];
w = [rLE(3,:); rTE(3,:)];

end