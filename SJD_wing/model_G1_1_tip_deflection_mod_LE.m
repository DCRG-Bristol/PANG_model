function Y = model_G1_1_tip_deflection_mod_LE(X)
%% Title section - Tip deflections versus wing root pitch angles and EI, GJ (G1_1 test)
%{
--------------------------------------------------------
Comments:
* This model refers to the G1_1 test (i.e., ground static test, clean wing)
* The model computes the tip deflection (leading edge (LE)) as a function of uncertain variables (i.e., EI, GJ)
* The wing root pitch angles are fixed to the experimental values used in the ground static test
* We use the UQLab notation: X: vector of uncertain variables; Y: vector of model outputs
--------------------------------------------------------
Input:
* X: multiplicative scaling factor for the uncertain variables
    * X(1)      : bending rigidity EI
    * X(2)      : torsional rigidity GJ    
--------------------------------------------------------
Output:
* Y: tip deflections corresponding to the different wing root angles (the angles arranged in increasing order of magnitude)
    * Y(1)      :static tip deflection (LE) at the smallest root angle (m)
    * Y(2)      :static tip deflection (LE) at the 2nd smallest root angle (m)
    ...
    * Y(5)      :static tip deflection (LE) at the largest root angle (m)
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

angl = exprData{1}.rootAngl;
angl = angl*pi/180; %transform from deg to rad

%function in the groundTests folder, labelled based on the test, take in
%the 'run' object and returns the test-specific measurements
[Y] = G1_1_tip_deflection(run, angl, 'EI', X(1), 'GJ', X(2));
Y = Y-Y(1);         % shift so that the deflection is zero at alpha=0

end