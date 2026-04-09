function [Y1, Y2, Y3, Y4, Y5] = model_W1_W2(X)
%% Title section - QIs versus wing root pitch angle and lambda_L, alpha_L, EI, GJ, Sxx (W1.1, W1.2, W2 tests)
%{
--------------------------------------------------------
Comments:
* This model refers to the W1 and W2 tests (i.e., wind tunnel tests)
* The model computes the static response, Hopf bifurcations speeds, wind-on frequencies and damping ratios versus the uncertain variables (i.e., lambda_L, alpha_L, EI, GJ, Sxx)
* We use the notation: X: vector of uncertain variables; Y1, ..., Y5: model outputs
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
* Y1        : Hopf bifurcations speeds (m/s)
    * Y1(1)      : first Hopf bifurcation speed (i.e., flutter speed) 
    * Y1(2)      : second Hopf bifurcation speed 
* Y2        : bending strains - static response (beta_y) (1/m)
    * Y2(1)     : bending strain at the first Hopf bifurcation speed 
    * Y2(2)     : bending strain at the second Hopf bifurcation speed 
* Y3        : torsional strains - static response (beta_x) (1/m)
    * Y3(1)     : torsional strain at the first Hopf bifurcation speed 
    * Y3(2)     : torsional strain at the second Hopf bifurcation speed 
* Y4        : wind-on frequencies (Hz) for increasing airspeed 
    * Y4(1, :)  : modal freq. for the first out-of-plane bending mode (OOP1)
    * Y4(2, :)  : modal freq. for the first in-plane bending mode (IP1)
    * Y4(3, :)  : modal freq. for the second out-of-plane bending mode (OOP2)
    * Y4(4, :)  : modal freq. for the first torsional mode (TOR1)
* Y5        : wind-on damping ratios (-) for increasing airspeed 
    * Y5(1, :)  : damping ratios for the first out-of-plane bending mode (OOP1)
    * Y5(2, :)  : damping ratios for the first in-plane bending mode (IP1)
    * Y5(3, :)  : damping ratios for the second out-of-plane bending mode (OOP2)
    * Y5(4, :)  : damping ratios for the first torsional mode (TOR1)
--------------------------------------------------------
%}
%%
%add folder below to path - contains functions relating to the experimental
%data and test-specific computations
addpath('WTTests'); 
load('run_ONERA.mat'); run = run_ONERA;

%set to work with modal form...
run = run.setTransform('modal', 14); 

ang = 1.0*pi/180; %wing root pitch angle to run..; transform from deg to rad

%function in the WTTests folder, labelled based on the test, take in
%the 'run' object and returns the test-specific measurements
[statResp, Y1, Y2, Y3] = W_statStab(run, false, ang, 0.275*X(1), 0.44*X(2), 'EI', X(3), 'GJ', X(4), 'Sxx', X(5), 'Szz', X(5));

Y4 = statResp.frqs;
Y5 = statResp.damp;

end