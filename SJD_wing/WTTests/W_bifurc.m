function [statResp, lcoResp, coco_log] = W_bifurc(obj, ang, lamL, alpL, varargin)

%function to call COCO (continuation core) to compute the full bifuration
%diagram (equib and LCO) at a fixed angle of attack 'ang'

%% fixed parameters...
s_I = 0.025; %measurment location for WRBM
s_II = 0.050;
dy1 = 0.075;
dy2 = -0.075;
opts = optimoptions('fsolve', 'Display', 'off');

%record and keep original valyes of ONERA parameters...
lam0 = obj.lam;
ML0 = obj.ML;

obj.lam = lamL;
obj.ML = alpL;

%outputs...
lcoResp = [];
statResp = [];
coco_log = [];

%% solve for alpha=0 case...
q = obj.q0_struct;
q = fsolve(@(q_all)(...
    obj.structDeriv(q_all,'alpha0', ang+0.8*pi/180, varargin{:})), q, opts);

%record wing root bending moments....
beta_y0 = project.basis.W2(s_I)'*obj.transF*q(obj.structDisp);
beta_x0 = project.basis.T1(s_II)'*obj.transF*q(obj.structDisp);

%% equilibrium continuation....

obj.dispFlag = false;

%ep tool box problem structure...
prob_equib = coco_prob();
prob_equib = coco_set(prob_equib,'cont','PtMX', [200,200],...
    'ResTOL',1e-3, 'TOL', 1e-1, 'h_min', 0.01, 'h_max', 5);
prob_equib = coco_set(prob_equib,'corr',...
    'ResTOL',1e-2, 'TOL', 1e-1);
prob_equib = coco_set(prob_equib,'ode','vectorized', false,...
    'ResTOL',1e-3, 'TOL', 1e-1);

U0 = 18; %initial velocity
q = fsolve(@(q_all)(...
    obj.aero_structDeriv(q_all,...
    'alpha0', ang+0.8*pi/180, 'alpha', ang, 'U', U0,...
    varargin{:})), obj.q0_aeroStruct, opts);

prob_equib = ode_isol2ep(prob_equib, '',...
    @(q,p)(...
    obj.aero_structDeriv(q,...
    'U', p(1), 'alpha', ang, 'alpha0', ang+0.8*pi/180, varargin{:})),...
    q, {'U'}, [U0]); %continuation variables and initial values...

coco_log.ep.name = ['ep_',datestr(now, 'HH_MM_SS_dd_mm_yy')];
coco_log.ang=ang;
coco_log.inputPars=varargin;

%call coco.....
coco(prob_equib, coco_log.ep.name, [],...
    1, {'U'}, [15, 30]); %run continuuation between 10<U<21

equib = coco_bd_read(coco_log.ep.name);
qStat = coco_bd_col(equib, 'x');
Ur = coco_bd_col(equib, 'U');

statResp.beta_y(1,:) = rbm(qStat);
statResp.beta_x(1,:) = rtm(qStat);
statResp.U=Ur;
hopf_idx = coco_bd_idxs(equib, 'HB');

if ~isempty(hopf_idx)
    statResp.beta_yf(1,:) = rbm(qStat(:,hopf_idx));
    statResp.beta_xf(1,:) = rtm(qStat(:,hopf_idx));
    statResp.Uf=Ur(hopf_idx);
else
    statResp.beta_yf(1,:) = [NaN, NaN];
    statResp.beta_xf(1,:) = [NaN, NaN];
    statResp.Uf=[NaN, NaN];
end

%% PO continuation

if ~isempty(hopf_idx) %..LCO if hopf exists...

    coco_log.po.name = ['po_',datestr(now, 'HH_MM_SS_dd_mm_yy')];
    LCO_prob = coco_prob();
    LCO_prob = coco_set(LCO_prob,'cont','PtMX', [0,75], 'MXCL', false,...
        'ResTOL',1e-3, 'TOL', 1e-3, 'h_min', 0.0001, 'h0', 0.001);
    LCO_prob = coco_set(LCO_prob,'corr', 'ItMX', 500, 'ResTOL',1e-3, 'TOL', 1e-3);
    LCO_prob = coco_set(LCO_prob, 'coll', 'NTST', 30,'NAdapt', 2,...
        'ResTOL',1e-3, 'TOL', 1e-3);
    LCO_prob = coco_set(LCO_prob,'ode','vectorized', false,...
        'ResTOL',1e-3, 'TOL', 1e-1);

    HB_labs = coco_bd_labs(equib, 'HB');
    LCO_prob = ode_HB2po(LCO_prob, '', coco_log.ep.name, HB_labs(1));

    [data uidx] = coco_get_func_data(LCO_prob, 'po.orb.coll', 'data', 'uidx');
    maps = data.coll_seg.maps;
    uidx = uidx([maps.xbp_idx]);
    LCO_prob = coco_add_func(LCO_prob,'mx_rbm', @mx_rbm_lco, data,...
        'regular','mx_rbm', 'uidx', uidx);
    LCO_prob = coco_add_func(LCO_prob,'mn_rbm', @mn_rbm_lco, data,...
        'regular','mn_rbm', 'uidx', uidx);

    coco(LCO_prob, coco_log.po.name, [], 1, {'U', 'mx_rbm, mn_rbm'}, [17.5,35]);

    po_bd = coco_bd_read(coco_log.po.name);

    %min/max bending response
    lcoResp.mx_beta_y = coco_bd_col(po_bd, 'mx_rbm');
    lcoResp.mn_beta_y = coco_bd_col(po_bd, 'mn_rbm');

    %lco airpseed vector..
    lcoResp.U = coco_bd_col(po_bd, 'U');
    lco_stab = coco_bd_col(po_bd,'po.test.USTAB'); %# unstable floquet multipliers
    lcoResp.isStab_LCO = (lco_stab==0); %..LO stability indicator...

    %lco prd..d..
    prd = coco_bd_col(po_bd, 'po.period'); %LCO period
    lcoResp.frq=1./prd;%..frequency

end

%% re set orignal ONERA parameters...

obj.lam = lam0;
obj.ML = ML0;

%% utility functions...

    function [data, y] = mx_rbm_lco(prob, data, u)
        shp = data.coll_seg.maps.xbp_shp;
        u_grp = reshape(u,shp);
        bd = rbm(u_grp);
        y = max(bd);
    end

    function [data, y] = mn_rbm_lco(prob, data, u)
        shp = data.coll_seg.maps.xbp_shp;
        u_grp = reshape(u,shp);
        bd = rbm(u_grp);
        y = min(bd);
    end

    function beta_y = rbm(qin)
        beta_y =...
            project.basis.W2(s_I)'*obj.transF*qin(obj.structDisp,:) - beta_y0;
    end

    function beta_x = rtm(qin)
        beta_x =...
            project.basis.T1(s_II)'*obj.transF*qin(obj.structDisp,:) - beta_x0;
    end


end