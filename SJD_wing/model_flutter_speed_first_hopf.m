function Y = model_flutter_speed_first_hopf(X,P)
%% Title section - Flutter speed (first Hopf point) versus AoA and E, G
%{
--------------------------------------------------------
Comments:
* The model computes the flutter speed (i.e., the first Hopf point) as a function of uncertain variables (i.e., E, G)
* We use a UQLab style notation, i.e., X: vector of uncertain variables; P: vector of deterministic parameters; Y: vector of model outputs
--------------------------------------------------------
Input:
* X: multiplicative scaling factor for the uncertain variables
    * X(1)      : Young's modulus E
    * X(2)      : torsional modulus G    
* P: deterministic parameters
    * P(1)      : angle of attack (rad)
--------------------------------------------------------
Output:
* Y: vector of quantities of interest (QIs)
    * Y(1)      : flutter speed (m/s)
--------------------------------------------------------
%}
%%
load('run_ONERA.mat'); run = run_ONERA;

%set initial parameters - note mach =0 for SJD wing
run = run.setPars('EI_1', X(1), 'EI_2', X(1), 'G', X(2), 'alpha0', 0*pi/180, 'alpha', 0, 'g', 9.81, 'mach', 0);

run.dispFlag = false; %this suppresses unnecessary message displaying

qstr0 = fsolve(@(q_all)(run.structDeriv(q_all,'alpha', 0, 'alpha0', 0)),...
    run.q0_struct);

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
U0 =  16; %initial speed...

%collect and initial solution for U=11 using fsolve...
q0 = fsolve(@(q_all)(run.aero_structDeriv(q_all, 'U', 11,...
    'alpha', P(1))),...
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
    {'U', 'alp'}, [U0, P(1)]); %continuation variables and initial values...

dt_date_time = datetime('now');
fmt_date_time = "HHmmssSSS";
str_date_time = string(dt_date_time, fmt_date_time);

%string for run id....
equib_run_id{1} = ['equilibrium_branch_AoA=',...
    % num2str(P(1)*180/pi,3)];
    str_date_time];

%call coco.....
bd_TrvEquib = coco(prob_equib, equib_run_id{1}, [],...
    1, {'U'}, [15, 35]); %run continuuation between 10<U<21

%.read coco solutions....see coco examples...
equib = coco_bd_read(equib_run_id{1});
q_stat{1} = coco_bd_col(equib, 'x');
U_stat{1} = coco_bd_col(equib, 'U');

%record flutter speeds...
hopf_idx = coco_bd_idxs(equib, 'HB');

if ~isempty(hopf_idx)
    U_1dim_HB{1} = U_stat{1}(hopf_idx);
else
    U_1dim_HB{1} = nan;
end

Y = U_1dim_HB{1}(1);
end