function [beta_y] = G1_1_beta_y(obj, ang, varargin)
 
% arguments
%     obj analysisBase
%     ang=0;
%     varargin
% end

%% fixed parameters...

s_I = 0.025; %measurment location for WRBM
s_II = 0.050;
dy1 = 0.075; 
dy2 = -0.075;
opts = optimoptions('fsolve', 'Display', 'off');

%% alpha=90 case for displ offset....

% alp=90*pi/180;
% q=obj.q0_struct;
% q = fsolve(@(q_all)(...
%     obj.structDeriv(q_all,'alpha0', alp, varargin{:})), q, opts);
% 
% [x,y,z] = obj.getDisplField(q, obj.geom.L,'beamModel');
% 
% delta_LE0 = z(1)*cos(alp) + y(1)*sin(alp) +...
%     -dy1*sin(alp);

% delta_TE0 = z(2)*cos(alp) + y(2)*sin(alp) +...
%     -dy2*sin(alp);

%% solve for alpha=0 case...

alp=0;
q=obj.q0_struct;
q = fsolve(@(q_all)(...
    obj.structDeriv(q_all,'alpha0', alp, varargin{:})), q, opts);

%record wing root bending moments....
beta_y0 = project.basis.W2(s_I)'*obj.transF*q(obj.structDisp);
% beta_x0 = project.basis.T1(s_II)'*obj.transF*q(obj.structDisp);

%% solve through angle vector...



for a_idx=1:length(ang)

    alp = ang(a_idx);

    if isempty(varargin)
        q = fsolve(@(q_all)(...
            obj.structDeriv(q_all,'alpha0', alp)),...
            q, opts);
    else
        q = fsolve(@(q_all)(...
            obj.structDeriv(q_all,'alpha0', alp,...
            varargin{:})),...
            q, opts);
    end

    % [x,y,z] = obj.getDisplField(q, obj.geom.L,'beamModel');

    % delta_LE(1,a_idx) = z(1)*cos(alp) + y(1)*sin(alp) +...
    %      -dy1*sin(alp) - delta_LE0;

    % delta_TE(1,a_idx) = z(2)*cos(alp) + y(2)*sin(alp) +...
    %     -dy2*sin(alp) - delta_TE0;

    %record wing root bending moments....
    beta_y(1,a_idx) = project.basis.W2(s_I)'*obj.transF*q(obj.structDisp)...
        - beta_y0;
    % beta_x(1,a_idx) = project.basis.T1(s_II)'*obj.transF*q(obj.structDisp)...
    %     - beta_x0;


end



end