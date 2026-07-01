function [Y1, Y2, Y3, Y4] = model_W1_3(X)
%% Title section - QIs versus wing root pitch angle and lambda_L, alpha_L, EI, GJ, Sxx (W1.3 test - LCOs)
%{
--------------------------------------------------------
Comments:
* This model refers to the W1.3 test (i.e., wind tunnel test for LCOs)
* The model computes the amplitudes of the LCO responses (bending, torsion) and the frequencies of the LCO orbits versus the uncertain variables (i.e., lambda_L, alpha_L, EI, GJ, Sxx)
* We use the notation: X: vector of uncertain variables; Y1, ..., Yn: model outputs
--------------------------------------------------------
Input:
* X: multiplicative scaling factor for the uncertain variables
    * X(1)      : ONERA lift coefficient lambda_L
    * X(2)      : ONERA lift coefficient alpha_L
    * X(3)      : bending rigidity EI
    * X(4)      : torsional rigidity GJ   
    * X(5)      : rotational inertia m_{Sxx} of fairing segments
--------------------------------------------------------
Outputs:
--------------------------------------------------------
* Y1        : amplitudes of the LCO bending responses (max response-min response) for increasing airspeed 
* Y2        : amplitudes of the LCO torsional responses (max-min) for increasing airspeed 
* Y3        : frequencies of the LCO orbits for increasing airspeed 
* Y4        : vector of increasing airspeeds at which the LCOs are observed
%}
%%
%add folder below to path - contains functions relating to the experimental
%data and test-specific computations
addpath('WTTests'); 
load('run_ONERA.mat'); run = run_ONERA;

%set to work with modal form...
run = run.setTransform('modal', 14); 

ang = 1.0*pi/180; %wing root pitch angle to run..; transform from deg to rad

%required: COCO version 2020 installed and initialized in path.
%the continuation processes are called across two steps..

%STEP1: equilibirum continuation: same function 'W_statStab' with the
%second entry (isCOCO) set to true. Collect 'coco_log' from output.
%W_statStab: function in the WTTests folder, labelled based on the test, take in
%the 'run' object and returns the test-specific measurements
[statResp, Uf, beta_yf, beta_xf, coco_log] = W_statStab(run, true, ang, 0.275*X(1), 0.44*X(2), 'EI', X(3), 'GJ', X(4), 'Sxx', X(5), 'Szz', X(5));
%STEP2: call continuation function: pass in the 'coco_log' from STEP1
[lcoResp, coco_log] = W_LCO(run, coco_log);

max_beta_y_LCO = lcoResp.mx_beta_y;
min_beta_y_LCO = lcoResp.mn_beta_y;
Y1 = max_beta_y_LCO - min_beta_y_LCO; %peak-to-peak amplitude
max_beta_x_LCO = lcoResp.mx_beta_x;
min_beta_x_LCO = lcoResp.mn_beta_x;
Y2 = max_beta_x_LCO - min_beta_x_LCO; %peak-to-peak amplitude
Y3 = lcoResp.frq;
Y4 = lcoResp.U;

end