
function [obj] = genBasis(obj)

fPath = ['+project\+basis\'];

%geometric,
L = obj.geom.L;
b = obj.geom.b;
a = obj.geom.a;

%basis sizes...
N = obj.basis;
Nw = N.Nw; Nv = N.Nv; Nthet = N.Nthet; 

%% aerodynamic grids..

NA = obj.basis.Ngam;
Ngam = NA;

%lifting line descritisations
xi = obj.basis.xi; %linspace(0.005,0.655,Ngam+1);
    dx = (xi(2:end)-xi(1:end-1)); %element widths
    xGrid = 0.5*(xi(2:end)+xi(1:end-1)); %control points
    xTrail = xi(2:end); %trailing vortex points

% % retrieving matrix for circs, GAM = X*A;
% X = eye(NA,NA);%cos((pi*xGrid'/L).*[0:1:NA-1]); 
% 
% aer.X = X; 

aer.xGrid = xGrid; aer.dx = dx; aer.xTrail = xTrail;
aer.bAer = b(xGrid); aer.NA = NA;
obj.aer=aer;

%% Shape functions

syms x 
assume(x, 'real')

%Scaled/Shifted Chebyshev polynomials (I) - faster
syms Y(x) %Shifting
s=-2*(x/L)+1;
gen1=1; gen2=s; 
Y(x)=[gen1; gen2];
if max([Nw,Nv,Nthet])>2
    for j=3:max([Nw,Nv,Nthet])
        Y(x)=[Y(x); 2*s*gen2-gen1];
        g1=gen2; g2=2*s*gen2-gen1;
        gen1=g1; gen2=g2;
    end
end
y=Y(x); yL=Y(L);

%flap bending functions and derivatives
symW=[y(1:Nw).*x^2; zeros(Nv+Nthet,1)]; %symW = symW.*(wdip./subs(symW, [x,0], [L,1]));
symW1 = diff(symW,x); symW2 = diff(symW1,x); symW3 = diff(symW2,x);

%chord bending functions and derivatives
symV = [zeros(Nw,1); y(1:Nv).*x^2; zeros(Nthet,1)]; %symV = symV.*(vdip./subs(symV, [x,0], [L,1]));
symV1 = diff(symV,x); symV2 = diff(symV1,x); symV3 = diff(symV2,x);

%planar pitching functions and derivatives
symThet = [zeros(Nw+Nv,1); y(1:Nthet)*x]; 
symThet1 = diff(symThet,x); symThet2 = diff(symThet1,x);

%convert to matlab functions...
W = matlabFunction(symW, 'File', [fPath,'W'], 'vars', {x}, 'Optimize', false); 
W1 = matlabFunction(symW1, 'File', [fPath,'W1'], 'vars', {x}, 'Optimize', false);
W2 = matlabFunction(symW2, 'File', [fPath,'W2'], 'vars', {x}, 'Optimize', false); 
W3 = matlabFunction(symW3, 'File', [fPath,'W3'], 'vars', {x}, 'Optimize', false);

V = matlabFunction(symV, 'File', [fPath,'V'], 'vars', {x}, 'Optimize', false);
V1 = matlabFunction(symV1, 'File', [fPath,'V1'], 'vars', {x}, 'Optimize', false);
V2 = matlabFunction(symV2, 'File', [fPath,'V2'], 'vars', {x}, 'Optimize', false); 
V3 = matlabFunction(symV3, 'File', [fPath,'V3'], 'vars', {x}, 'Optimize', false);

T = matlabFunction(symThet, 'File', [fPath,'T'], 'vars', {x}, 'Optimize', false);
T1 = matlabFunction(symThet1, 'File', [fPath,'T1'], 'vars', {x}, 'Optimize', false);

fw = int(symW1.*symW1', [0,x]);
fv = int(symV1.*symV1', [0,x]);

%subfunctions for horizontal displacements
Fw = matlabFunction(fw, 'File', [fPath,'fw'], 'vars', {x}, 'Optimize', false);
Fv = matlabFunction(fv, 'File', [fPath,'fv'], 'vars', {x}, 'Optimize', false);

%% Basis and save

% obj.basis.w = {W,W1,W2}; obj.basis.v = {V,V1,V2}; obj.basis.t = {T,T1};
% obj.basis.u = {Fw,Fv};

%% functions for plotting...
% 
% x_cont = obj.basis.xi;
% for j=1:length(x_cont)
%     uv{j} = project.basis.fv(x_cont(j)); uw{j} = project.basis.fw(x_cont(j));
% end
% idx_tot = 1:1:obj.basis.Ntot;
% %displaements..
% for j=1:length(x_cont)
%     x = x_cont(j);
% 
%     if j==1
%         w = @(q)([project.basis.W(x)'*q(idx_tot,1) + b(x)*(1+a(x))*((project.basis.T(x)'*q(idx_tot,1))*(1-0.5*(project.basis.W1(x)'*q(idx_tot,1))^2) +...
%             -(1/6)*(project.basis.T(x)'*q(idx_tot,1))^3 )]);
% 
%         v = @(q)([project.basis.V(x)'*q(idx_tot,1) + b(x)*(1+a(x))*(1-0.5*(project.basis.T(x)'*q(idx_tot,1))^2+...
%             -0.5*(project.basis.V1(x)'*q(idx_tot,1))^2 - (project.basis.V1(x)'*q(idx_tot,1))*(project.basis.T(x)'*q(idx_tot,1))*(project.basis.W1(x)'*q(idx_tot,1)))]);
% 
%         u = @(q)([x - 0.5*(q(idx_tot,1)'*uw{j}*q(idx_tot,1) + q(idx_tot,1)'*uv{j}*q(idx_tot,1) ) +...
%             b(x)*(1+a(x))*(-project.basis.V1(x)'*q(idx_tot,1)*(1-0.5*(project.basis.T(x)'*q(idx_tot,1))^2 + 0.5*(project.basis.W1(x)'*q(idx_tot,1))^2) +...
%             -(project.basis.W1(x)'*q(idx_tot,1))*(project.basis.T(x)'*q(idx_tot,1))) ]);
%     else
% 
%         w = @(q)([w(q);...
%             project.basis.W(x)'*q(idx_tot,1) + b(x)*(1+a(x))*((project.basis.T(x)'*q(idx_tot,1))*(1-0.5*(project.basis.W1(x)'*q(idx_tot,1))^2) +...
%             -(1/6)*(project.basis.T(x)'*q(idx_tot,1))^3 )]);
% 
%         v = @(q)([v(q);...
%             project.basis.V(x)'*q(idx_tot,1) + b(x)*(1+a(x))*(1-0.5*(project.basis.T(x)'*q(idx_tot,1))^2+...
%             -0.5*(project.basis.V1(x)'*q(idx_tot,1))^2 - (project.basis.V1(x)'*q(idx_tot,1))*(project.basis.T(x)'*q(idx_tot,1))*(project.basis.W1(x)'*q(idx_tot,1)))]);
% 
%         u = @(q)([u(q);...
%             x - 0.5*(q(idx_tot,1)'*uw{j}*q(idx_tot,1) + q(idx_tot,1)'*uv{j}*q(idx_tot,1) ) +...
%             b(x)*(1+a(x))*(-project.basis.V1(x)'*q(idx_tot,1)*(1-0.5*(project.basis.T(x)'*q(idx_tot,1))^2 + 0.5*(project.basis.W1(x)'*q(idx_tot,1))^2) +...
%             -(project.basis.W1(x)'*q(idx_tot,1))*(project.basis.T(x)'*q(idx_tot,1))) ]);
%     end
% end
% 
% for j=length(x_cont):-1:1
%     x = x_cont(j);
% 
%     w = @(q)([w(q);...
%         project.basis.W(x)'*q(idx_tot,1) - b(x)*(a(x)-1)*((project.basis.T(x)'*q(idx_tot,1))*(1-0.5*(project.basis.W1(x)'*q(idx_tot,1))^2) +...
%     -(1/6)*(project.basis.T(x)'*q(idx_tot,1))^3 )]);
% 
%     v = @(q)([v(q);...
%         project.basis.V(x)'*q(idx_tot,1) - b(x)*(a(x)-1)*(1-0.5*(project.basis.T(x)'*q(idx_tot,1))^2+...
%     -0.5*(project.basis.V1(x)'*q(idx_tot,1))^2 - (project.basis.V1(x)'*q(idx_tot,1))*(project.basis.T(x)'*q(idx_tot,1))*(project.basis.W1(x)'*q(idx_tot,1)))]);
% 
%     u = @(q)([u(q);...
%         x - 0.5*(q(idx_tot,1)'*uw{j}*q(idx_tot,1) + q(idx_tot,1)'*uv{j}*q(idx_tot,1) ) +...
%         -b(x)*(a(x)-1)*(-project.basis.V1(x)'*q(idx_tot,1)*(1-0.5*(project.basis.T(x)'*q(idx_tot,1))^2 + 0.5*(project.basis.W1(x)'*q(idx_tot,1))^2) +...
%         -(project.basis.W1(x)'*q(idx_tot,1))*(project.basis.T(x)'*q(idx_tot,1))) ]);
% end
% 
% Q = sym('q',[obj.basis.Ntot,1]); assume(Q, 'real');
% matlabFunction(w(Q), 'File', [fPath,'zMesh'], 'vars', {Q},...
%         'Optimize', true);
% matlabFunction(v(Q), 'File', [fPath,'yMesh'], 'vars', {Q},...
%         'Optimize', true);
% matlabFunction(u(Q), 'File', [fPath,'xMesh'], 'vars', {Q},...
%         'Optimize', true);

% plotF.W = w; plotF.V = v; plotF.U = u;
% plotF.zag = zag;
% 
% plotF.arm.comp = {xx, yy, zz}; plotF.arm.atch = atch;
% 
% save('model.mat', 'plotF', '-append');

end