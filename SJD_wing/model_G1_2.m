function Y = model_G1_2(X,angl)
%% Title section - Tip deflections and strains versus wing root pitch angles and EI, GJ (G1_2 test)
%{
--------------------------------------------------------
Comments:
* This model refers to the G1_2 test (i.e., ground static test with wing tip weight)
* The model computes the tip deflection (leading edge (LE) and trailing edge (TE)) as a function of uncertain variables (i.e., EI, GJ)
* The model also computes the strains (bending (beta_y) and torsional (beta_x)) as a function of uncertain variables (i.e., EI, GJ)
* We use the notation: X: vector of uncertain variables; angl: vector of deterministic parameters; Y: vector of model outputs
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
* Y: 4xn matrix of quantities of interest (QIs)
    * Y(1, :)      :static tip deflection (LE) (m)
    * Y(2, :)      :static tip deflection (TE) (m)
    * Y(3, :)      :bending strain (beta_y) (1/m)
    * Y(4, :)      :torsional strain (beta_x) (1/m)
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
[delta_LE, delta_TE, beta_y, beta_x] = G1_2(run, angl, 'EI', X(1), 'GJ', X(2));

Y = zeros(4, length(angl));
Y(1, :) = delta_LE;
Y(2, :) = delta_TE;
Y(3, :) = beta_y;
Y(4, :) = beta_x;

end