function Y = model_G1_1_tip_deflection_LE(X,angl)
%% Title section - Tip deflections versus wing root pitch angles and EI, GJ (G1_1 test)
%{
--------------------------------------------------------
Comments:
* This model refers to the G1_1 test (i.e., ground static test, clean wing)
* The model computes the tip deflection (leading edge (LE)) as a function of uncertain variables (i.e., EI, GJ)
* We use the UQLab notation: X: vector of uncertain variables; angl: vector of deterministic parameters; Y: vector of model outputs
--------------------------------------------------------
Input:
* X: multiplicative scaling factor for the uncertain variables
    * X(1)      : bending rigidity EI
    * X(2)      : torsional rigidity GJ    
* angl: row vector of n deterministic parameters (please arrange in increasing order of magnitude)
    * angl(1)   : smallest wing root pitch angle (deg) (e.g., angl(1)=0)
    ...
    * angl(n)   : largest wing root pitch angle (deg) (e.g., angl(n)=60)
--------------------------------------------------------
Output:
* Y: 1xn row vector of quantities of interest (QIs)
    * Y(1, :)      :static tip deflection (LE) (m)
--------------------------------------------------------
%}
%%
%add folder below to path - contains functions relating to the experimental
%data and test-specific computations
addpath('groundTests'); 
load('run_ONERA.mat'); run = run_ONERA;

%set to work with modal form...
run = run.setTransform('modal', 14); 

angl = angl*pi/180; %transform from deg to rad

%function in the groundTests folder, labelled based on the test, take in
%the 'run' object and returns the test-specific measurements
[delta_LE, ~, ~, ~] = G1_1(run, angl, 'EI', X(1), 'GJ', X(2));

Y = zeros(1, length(angl));
Y(1, :) = delta_LE;

end