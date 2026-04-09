function Y = model_W1_2_Uf_all_unc(X)
%% Title section - Hopf bifurcations speeds versus lambda_L, alpha_L, EI, GJ, Sxx, Szz (W1.2 test)
%{
--------------------------------------------------------
Comments:
* This model refers to the W1.2 test (i.e., wind tunnel test)
* The model computes the Hopf bifurcations speeds as a function of uncertain variables (i.e., lambda_L, alpha_L, EI, GJ, Sxx, Szz)
* The wing root pitch angle is fixed to a value close to one experimental value used in the WT test
* We use the notation: X: vector of uncertain variables; Y: vector of model outputs
--------------------------------------------------------
Input:
* X: multiplicative scaling factor for the uncertain variables
    * X(1)      : ONERA lift coefficient lambda_L
    * X(2)      : ONERA lift coefficient alpha_L
    * X(3)      : bending rigidity EI
    * X(4)      : torsional rigidity GJ   
    * X(5)      : rotational inertia m_{Sxx} of fairing segments
--------------------------------------------------------
Output:
* Y: Hopf bifurcations speeds (m/s)
    * Y(1)      :first Hopf bifurcation speed (i.e., flutter speed)
    * Y(2)      :second Hopf bifurcation speed
--------------------------------------------------------
%}
%%
%add folder below to path - contains functions relating to the experimental
%data and test-specific computations
addpath('WTTests'); 
load('run_ONERA.mat'); run = run_ONERA;

%set to work with modal form...
run = run.setTransform('modal', 14); 

ang = 1.0*pi/180; %angle to run..; transform from deg to rad

%function in the WTTests folder, labelled based on the test, take in
%the 'run' object and returns the test-specific measurements
[~, Y, ~, ~] = W_statStab_parfor(run, false, ang, X(1), X(2), 'EI', X(3), 'GJ', X(4), 'Sxx', X(5), 'Szz', X(5));

end