function [U_bd, beta_ybd, beta_xbd] = W_stabEnvel(obj, coco_log)

%% function to compute the entire stability envelope.....

%% fixed parameters...

s_I = 0.025; %measurment location for WRBM
s_II = 0.050;
dy1 = 0.075; 
dy2 = -0.075;
opts = optimoptions('fsolve', 'Display', 'off');

%% solve for alpha=0 case...

q=obj.q0_struct;
q = fsolve(@(q_all)(...
    obj.structDeriv(q_all,'alpha0', 0, coco_log.inputPars{:})), q, opts);

%record wing root bending moments....
beta_y0 = project.basis.W2(s_I)'*obj.transF*q(obj.structDisp);
beta_x0 = project.basis.T1(s_II)'*obj.transF*q(obj.structDisp);

%%

prob_HB = coco_prob();
prob_HB = coco_set(prob_HB,'cont','PtMX', [-500,500],...
    'ResTOL',1e-3, 'h_min', 0.1,'h_max', 1000);
prob_HB = coco_set(prob_HB,'corr','PtMX', [0,2000], 'ItMX', 250,...
    'ResTOL',1e-3, 'h_max', 1000);

equib = coco_bd_read(coco_log.name);
HB_labs = coco_bd_labs(equib, 'HB');
prob_HB = ode_HB2HB(prob_HB, '', coco_log.name, HB_labs(1));

hb2name = [coco_log.name, '_hb2hb'];

bd_HB = coco(prob_HB, hb2name, [], 1, {'U', 'alp'},...
    {[20, 30],[0, 5]*pi/180});

%% post proc..

HBEquib = coco_bd_read(hb2name);
Ur = coco_bd_col(HBEquib, 'U');
q = coco_bd_col(HBEquib, 'x');

for i=1:length(Ur)
    beta_ybd(1,u_idx) = rbm(q(:,i));
    beta_xbd(1,u_idx) = rtm(q(:,i));
end
U_bd = Ur;

%% utility functions...
    function beta_y = rbm(qin)
        beta_y =...
            project.basis.W2(s_I)'*obj.transF*qin(obj.structDisp) - beta_y0;
    end

    function beta_x = rtm(qin)
        beta_x =...
            project.basis.T1(s_II)'*obj.transF*qin(obj.structDisp) - beta_x0;
    end

end