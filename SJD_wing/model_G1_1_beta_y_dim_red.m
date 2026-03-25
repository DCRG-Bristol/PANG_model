function Y = model_G1_1_beta_y_dim_red(X)
%% Title section - Strain versus wing root pitch angles and EI (G1_1 test)
%{
--------------------------------------------------------
Comments:
* This model refers to the G1_1 test (i.e., ground static test, clean wing)
* The model computes the strain (bending (beta_y)) as a function of uncertain variables (i.e., EI)
* The variable GJ is fixed here (dimension reduction), as a result of sensitivity analysis
* The wing root pitch angles are fixed to the experimental values used in the ground static test
* We use the UQLab notation: X: vector of uncertain variables; Y: vector of model outputs
--------------------------------------------------------
Input:
* X: multiplicative scaling factor for the uncertain variables
    * X(1)      : bending rigidity EI
--------------------------------------------------------
Output:
* Y: bending strains corresponding to the different wing root angles (the angles arranged in increasing order of magnitude)
    * Y(1)      :bending strain (beta_y) (1/m) at the smallest root angle
    * Y(2)      :bending strain (beta_y) (1/m) at the 2nd smallest root angle
    ...
    * Y(5)      :bending strain (beta_y) (1/m) at the largest root angle
--------------------------------------------------------
%}
%%
%add folder below to path - contains functions relating to the experimental
%data and test-specific computations
addpath('groundTests'); 
load('run_ONERA.mat'); run = run_ONERA;
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');

%set to work with modal form...
run = run.setTransform('modal', 14); 

angl = exprData{1}.rootAngl(2:end);
angl = angl*pi/180; %transform from deg to rad

%function in the groundTests folder, labelled based on the test, take in
%the 'run' object and returns the test-specific measurements
[Y] = G1_1_beta_y(run, angl, 'EI', X(1), 'GJ', 1);

end