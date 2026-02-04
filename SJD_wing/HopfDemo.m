clear all; clc;


EI1 = 2.268; %EI stored in this variable is used to compute My
s_I = 0.025; %measurment location for WRBM

load('run_ONERA.mat'); run = run_ONERA;

%set initial parameters - note mach =0 for SJD wing
run = run.setPars('EI_1', 1.1, 'EI_2', 1, 'G', 0.9, 'alpha0', 0*pi/180, 'alpha', 0, 'g', 9.81, 'mach', 0);

run.dispFlag = false; %this suppresses unnecessary message displaying

qstr0 = fsolve(@(q_all)(run.structDeriv(q_all,'alpha', 0, 'alpha0', 0)),...
    run.q0_struct);

%wing root bnding momment at self-weight deformation...
My0 = EI1*project.basis.W2(s_I)'*run.transF*qstr0(run.structDisp);

%get linear modes and matrices....
[shp, evals, Kmat, Cmat, Mmat] = run.getStructModes(qstr0);

%add damping (values from paper)
Dmat = 0.071*Mmat + (10e-5)*Kmat;
run.dampMatr = Dmat;

%set to work with modal form...
run = run.setTransform('modal', 14); %...using the first 6 modes

%% equilirbium tracing.. using coco...

%this block performs continution alng an equilibrium path with increasing
%airspeed and identifies the Hopf bifurcations. This uses the 'ep' toolbox
% in the continuation core (COCO) package - see
% (https://sourceforge.net/projects/cocotools/)

%starting points..
alp0 = ([0.3, 1.0, 1.2])*pi/180; %angles of attack to run...
U0 =  16; %initial speed...

%figure for plotting equilibrium paths....
figure;
tipDispAx = subplot(1,2,1);
xlabel('Airspeed, [U]'); ylabel('w_L/L, [-]'); hold on;

dMyAx = subplot(1,2,2);
xlabel('Airspeed, [U]'); ylabel('\DeltaM_y, [Nm]'); hold on;

for a_idx=1%:length(alp0)

    %collect and initial solution for U=11 using fsolve...
    q0 = fsolve(@(q_all)(run.aero_structDeriv(q_all, 'U', 11,...
        'alpha', alp0(a_idx))),...
        run.q0_aeroStruct);

    %ep tool box problem structure...
    %see coco po example...
    prob_equib = coco_prob();
    prob_equib = coco_set(prob_equib,'cont','PtMX', [200,200]);
    prob_equib = coco_set(prob_equib,'ode','vectorized', false);

    %see coco po example....
    prob_equib = ode_isol2ep(prob_equib, '',...
        @(q,p)(...
        run.aero_structDeriv(q, 'U', p(1), 'alpha', p(2), 'alpha0', p(2))),... %state space system
        q0,... %initial estimae
        {'U', 'alp'}, [U0, alp0(a_idx)]); %continuation variables and initial values...

    %string for run id....
    equib_run_id{a_idx} = ['equilibrium_branch_AoA=',...
        num2str(alp0(a_idx)*180/pi,3)];

    %call coco.....
    bd_TrvEquib = coco(prob_equib, equib_run_id{a_idx}, [],...
        1, {'U'}, [15, 31]); %run continuuation between 10<U<21

    %.read coco solutions....see coco examples...
    equib = coco_bd_read(equib_run_id{a_idx});
    q_stat{a_idx} = coco_bd_col(equib, 'x');
    U_stat{a_idx} = coco_bd_col(equib, 'U');

    %plot equilibrium path....
    for i=1:length(q_stat{a_idx}(1,:))
        [x,y,z] = run.getDisplField(q_stat{a_idx}(:,i), run.geom.L,...
            'beamModel');
        tipDisp_equib{a_idx}(i) = z(1,1)/run.geom.L;  
    end

    %record wing root bending moments....
    dMy_equib{a_idx} =...
        EI1*project.basis.W2(s_I)'*...
        run.transF*q_stat{a_idx}(run.structDisp,:) - My0;

    %record flutter speeds...
    hopf_idx = coco_bd_idxs(equib, 'HB');

    if ~isempty(hopf_idx)
        tipDisp_1dim_HB{a_idx} = tipDisp_equib{a_idx}(hopf_idx);
        U_1dim_HB{a_idx} = U_stat{a_idx}(hopf_idx);
        dMy_1dim_HB{a_idx} = dMy_equib{a_idx}(hopf_idx);
    else
        tipDisp_1dim_HB{a_idx} = [];
        U_1dim_HB{a_idx} = [];
        dMy_1dim_HB{a_idx} = [];
    end

    plot(tipDispAx, U_stat{a_idx}, tipDisp_equib{a_idx}, 'k-');
    hold(tipDispAx, 'on');
    plot(tipDispAx, U_1dim_HB{a_idx}, tipDisp_1dim_HB{a_idx},...
        'k^', 'markerFaceColor', 'r');
    hold(tipDispAx, 'on');

    plot(dMyAx, U_stat{a_idx}, dMy_equib{a_idx}, 'k-');
    hold(dMyAx, 'on');
    plot(dMyAx, U_1dim_HB{a_idx}, dMy_1dim_HB{a_idx},...
        'k^', 'markerFaceColor', 'r');
    hold(dMyAx, 'on');
    drawnow;
end

%the follolowing could be read as teh output at this point....
%alp0 - angles of attach used to run...index a_idx=1,...
%U_stat{a_idx} - airspeed vector
%tipDisp_equib{a_idx} - tip deflection vector
%tipDisp_1dim_HB{a_idx} - tip deflections at hopf points...
%dMy_1dim_HB{a_idx} - bending moments at hopf points...
%U_1dim_HB{a_idx} - hopf airspeeds



%% co-dim2 HB tracing...

%old coco (2020)
prob_HB = coco_prob();
prob_HB = coco_set(prob_HB,'cont','PtMX', [-500,500],...
    'ResTOL',1e-2, 'h_min', 10,'h_max', 1000);
prob_HB = coco_set(prob_HB,'corr','PtMX', [0,2000], 'ItMX', 250,...
    'ResTOL',1e-2, 'h_max', 1000);

%settings for new COCO (2025)
% prob_HB = coco_prob();
% prob_HB = coco_set(prob_HB,'cont','PtMX', [-500,500],...
%     'ResTOL',1e-4, 'TOL',1e-4, 'h_min', 50,'h_max', 500);
% prob_HB = coco_set(prob_HB,'corr','PtMX', [0,2000], 'ItMX', 250,...
%     'ResTOL',1e-4, 'TOL',1e-4, 'h_max', 100);

HB_labs = coco_bd_labs(bd_TrvEquib, 'HB');
prob_HB = ode_HB2HB(prob_HB, '', equib_run_id{end}, HB_labs(1));

%runner call...Ip
bd_HB = coco(prob_HB, 'Flutter_boundary', [], 1, {'U', 'alp'},...
    {[12.5, 35],[0.8, 5]*pi/180-0.8*pi/180});

%% .collect alpha-U results...

HBEquib = coco_bd_read('Flutter_boundary');

U_flut = coco_bd_col(HBEquib, 'U');
alp_flut = coco_bd_col(HBEquib, 'alp');
q_flut = coco_bd_col(HBEquib, 'x');

%get displacements...
for i=1:length(q_flut(1,:))
    [x,y,z] = run.getDisplField(q_flut(:,i), run.geom.L,...
        'beamModel');
    tipDisp_flut(i) = z(1,1)/run.geom.L;
end

plot(tipDispAx, U_flut, tipDisp_flut, 'r-', 'linewidth', 2);
