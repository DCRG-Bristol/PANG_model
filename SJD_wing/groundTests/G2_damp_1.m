function [damp] = G2_damp_1(obj, ang, varargin)
 
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

%% solve for alpha=0 case...

q=obj.q0_struct;
q = fsolve(@(q_all)(...
    obj.structDeriv(q_all,'alpha0', 0, varargin{:})), q, opts);

%record wing root bending moments....
beta_y0 = project.basis.W2(s_I)'*obj.transF*q(obj.structDisp);
beta_x0 = project.basis.T1(s_II)'*obj.transF*q(obj.structDisp);

%% solve through angle vector...



for a_idx=1:length(ang)

    alp = ang(a_idx);

    if isempty(varargin)
        q = fsolve(@(q_all)(...
            obj.structDeriv(q_all,'alpha0', alp)),...
            q, opts);
        [shp, evals, Kmat, Cmat, Mmat] = obj.getStructModes(q);
    else
        q = fsolve(@(q_all)(...
            obj.structDeriv(q_all,'alpha0', alp,...
            varargin{:})),...
            q, opts);
        [shp, evals, Kmat, Cmat, Mmat] = obj.getStructModes(q,...
            varargin{:});
    end
    evals=evals(:);

    frqs(1,a_idx) = abs(evals(1))/(2*pi);
    damp(1,a_idx) = -real(evals(1))./abs(evals(1));
end



end