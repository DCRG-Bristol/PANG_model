function Y = model_G2_damp_4_dim_red(X)
%% Title section - Fourth damping ratio versus wing root pitch angles and GJ, Sxx, Szz (G2 test)
%{
--------------------------------------------------------
Comments:
* This model refers to the G2 test (i.e., ground vibration test)
* The model computes the 4th damping ratio (TOR1 torsion) as a function of uncertain variables (i.e., EI, GJ, Sxx, Szz)
* We use the notation: X: vector of uncertain variables; angl: vector of deterministic parameters; Y: vector of model outputs
--------------------------------------------------------
Input:
* X: multiplicative scaling factor for the uncertain variables
    * X(1)      : torsional rigidity GJ   
    * X(2)      : rotational inertia m_{Sxx} of fairing segments
--------------------------------------------------------
Output:
* Y: 4th damping ratio corresponding to the different wing root angles (the angles arranged in increasing order of magnitude)
    * Y(1)      :4th damping ratio at the smallest root angle (Hz)
    * Y(2)      :4th damping ratio at the 2nd smallest root angle (Hz)
    ...
    * Y(6)      :4th damping ratio at the largest root angle (Hz)
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

angl = exprData{1, 3}.rootAngl;
angl = angl*pi/180; %transform from deg to rad

%function in the groundTests folder, labelled based on the test, take in
%the 'run' object and returns the test-specific measurements
[Y] = G2_damp_4(run, angl, 'EI', 1, 'GJ', X(1), 'Sxx', X(2), 'Szz', X(2));

end