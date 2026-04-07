function [statResp, Uf, beta_yf, beta_xf, coco_log] = W_statStab(obj, isCont, ang, lamL, alpL, varargin)

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

coco_log = [];

%% solve for alpha=0 case...

q=obj.q0_struct;
q = fsolve(@(q_all)(...
    obj.structDeriv(q_all,'alpha0', 0, varargin{:})), q, opts);

%record wing root bending moments....
beta_y0 = project.basis.W2(s_I)'*obj.transF*q(obj.structDisp);
beta_x0 = project.basis.T1(s_II)'*obj.transF*q(obj.structDisp);


q = obj.q0_aeroStruct;
beta_yf = [NaN, NaN];
beta_xf = [NaN, NaN];
Uf = [NaN, NaN];

if isCont
    %% solution with numerical continuation.... (COC0)
    obj.dispFlag = false;

    %ep tool box problem structure...
    prob_equib = coco_prob();
    prob_equib = coco_set(prob_equib,'cont','PtMX', [200,200],...
        'ResTol',1e-3, 'Tol', 1e-1, 'h_min', 0.01, 'h_max', 5);
    prob_equib = coco_set(prob_equib,'corr',...
        'ResTol',1e-2, 'Tol', 1e-1);
    prob_equib = coco_set(prob_equib,'ode','vectorized', false,...
        'ResTol',1e-3, 'Tol', 1e-1);

    U0 = 11; %initial velocity

    prob_equib = ode_isol2ep(prob_equib, '',...
        @(q,p)(...
        obj.aero_structDeriv(q,...
        'U', p(1), 'alpha', p(2), 'alpha0', p(2)+0.8*pi/180, varargin{:})),... 
        q,... %initial estimate
        {'U', 'alp'}, [U0, ang]); %continuation variables and initial values...

    coco_log.name = ['eq_',datestr(now, 'HH_MM_SS_dd_mm_yy')];

    %call coco.....
    bd_TrvEquib = coco(prob_equib, coco_log.name, [],...
        1, {'U'}, [10, 29]); %run continuuation between 10<U<21

    equib = coco_bd_read(coco_log.name);
    q = coco_bd_col(equib, 'x');
    eigs_all = coco_bd_col(equib, 'eigs');
    Ur = coco_bd_col(equib, 'U');
    hopf_idx = coco_bd_idxs(equib, 'HB');

    statResp.U = Ur;

    for u_idx=1:length(Ur)
        %record static responses>>>>>>>>>>>>>>>>>>>>
        statResp.beta_y(1,u_idx) = rbm(q(:,u_idx));
        statResp.beta_x(1,u_idx) = rtm(q(:,u_idx));

        %select eigenvalues>>>>>>>>>>>>>>>>>>>>>>>>>>
        D = eigs_all(:,u_idx);
        %find oscillatory modes..
        osc_idx = find(abs(imag(D))>0.01*abs(real(D)));
        D = D(osc_idx);

        %sort by frequency..
        [~,odr] = sort(abs(D));
        eval=D(odr(2:2:end));

        statResp.frqs(:,u_idx) = abs(eval(1:4))/(2*pi);
        statResp.damp(:,u_idx) = -real(eval(1:4))./abs(eval(1:4));
    end

    for ii=1:min([2, length(hopf_idx)])
        Uf(1,ii) = Ur(hopf_idx(ii));
        beta_yf(1,ii) = statResp.beta_y(1,hopf_idx(ii));
        beta_xf(1,ii) = statResp.beta_x(1,hopf_idx(ii));
    end

    coco_log.inputPars = varargin;
    coco_log.ang = ang;

else

    %% solution based on manual computation...

    Ur = [linspace(10, 19, 10), linspace(20, 29, 30)];
    isCurrFlut =  false; %is flutter happening currently
    statResp.U = Ur;

    for u_idx = 1:length(Ur)

        q = fsolve(@(q_all)(...
            obj.aero_structDeriv(q_all,...
            'alpha0', ang+0.8*pi/180, 'alpha', ang, 'U', Ur(u_idx),...
            varargin{:})), q, opts);

        [eval, evec, isFlut] = obj.AE_modes(q,...
            'alpha0', ang+0.8*pi/180, 'alpha', ang, 'U', Ur(u_idx),...
            varargin{:});

        eval=eval(:);

        %record static responses...
        statResp.beta_y(1,u_idx) = rbm(q);
        statResp.beta_x(1,u_idx) = rtm(q);
        statResp.isStab(1,u_idx) = ~isFlut;
        statResp.frqs(:,u_idx) = abs(eval(1:4))/(2*pi);
        statResp.damp(:,u_idx) = -real(eval(1:4))./abs(eval(1:4));

        %going into flutter.....
        if ~isCurrFlut
            if isFlut
                %find mode in flutter...
                mode = find(statResp.damp(:,u_idx)<0);

                cr = statResp.damp(mode,u_idx);
                cr0 = statResp.damp(mode,u_idx-1);

                Uf(1,1) = Ur(u_idx-1) +...
                    (0-cr0)*(Ur(u_idx)-Ur(u_idx-1))/(cr-cr0);

                beta_yf(1,1) = statResp.beta_y(1,u_idx-1) +...
                    (0-cr0)*(statResp.beta_y(u_idx)-statResp.beta_y(u_idx-1))/(cr-cr0);

                beta_xf(1,1) = statResp.beta_x(1,u_idx-1) +...
                    (0-cr0)*(statResp.beta_x(u_idx)-statResp.beta_x(u_idx-1))/(cr-cr0);

                isCurrFlut = true;
            end
        end

        %coming out of flutter...
        if isCurrFlut
            if ~isFlut
                %find mode that was in flutter...
                mode = find(statResp.damp(:,u_idx-1)<0);

                cr = statResp.damp(mode,u_idx);
                cr0 = statResp.damp(mode,u_idx-1);

                Uf(1,2) = Ur(u_idx-1) +...
                    (0-cr0)*(Ur(u_idx)-Ur(u_idx-1))/(cr-cr0);

                beta_yf(1,2) = statResp.beta_y(1,u_idx-1) +...
                    (0-cr0)*(statResp.beta_y(u_idx)-statResp.beta_y(u_idx-1))/(cr-cr0);

                beta_xf(1,2) = statResp.beta_x(1,u_idx-1) +...
                    (0-cr0)*(statResp.beta_x(u_idx)-statResp.beta_x(u_idx-1))/(cr-cr0);

                isCurrFlut = false;
            end
        end
    end
end


%% re set orignal ONERA parameters...

obj.lam = lam0;
obj.ML = ML0;

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