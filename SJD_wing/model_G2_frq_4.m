function Y = model_G2_frq_4(X)
%% Title section - Fourth frequency versus wing root pitch angles and EI, GJ, Sxx, Szz (G2 test)
%{
--------------------------------------------------------
Comments:
* This model refers to the G2 test (i.e., ground vibration test)
* The model computes the 4th modal frequency (torsion) as a function of uncertain variables (i.e., EI, GJ, Sxx, Szz)
* The wing root pitch angles are fixed to the experimental values used in the GVT test
* We use the notation: X: vector of uncertain variables; Y: vector of model outputs
--------------------------------------------------------
Input:
* X: multiplicative scaling factor for the uncertain variables
    * X(1)      : bending rigidity EI
    * X(2)      : torsional rigidity GJ   
    * X(3)      : rotational inertia m_{Sxx} of fairing segments
    * X(4)      : rotational inertia m_{Szz} of fairing segments
--------------------------------------------------------
Output:
* Y: 4th modal frequency corresponding to the different wing root angles (the angles arranged in increasing order of magnitude)
    * Y(1)      :4th modal frequency at the smallest root angle (Hz)
    * Y(2)      :4th modal frequency at the 2nd smallest root angle (Hz)
    ...
    * Y(6)      :4th modal frequency at the largest root angle (Hz)
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
[Y] = G2_frq_4(run, angl, 'EI', X(1), 'GJ', X(2), 'Sxx', X(3), 'Szz', X(4));

end