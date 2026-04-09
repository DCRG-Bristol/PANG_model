function Y = model_G2(X,angl)
%% Title section - Frequencies versus wing root pitch angles and EI, GJ, Sxx, Szz (G2 test)
%{
--------------------------------------------------------
Comments:
* This model refers to the G2 test (i.e., ground vibration test)
* The model computes the lowest four modal frequencies as a function of uncertain variables (i.e., EI, GJ, Sxx, Szz)
* We use the notation: X: vector of uncertain variables; angl: vector of deterministic parameters; Y: vector of model outputs
--------------------------------------------------------
Input:
* X: multiplicative scaling factor for the uncertain variables
    * X(1)      : bending rigidity EI
    * X(2)      : torsional rigidity GJ   
    * X(3)      : rotational inertia m_{Sxx} of fairing segments
    * X(4)      : rotational inertia m_{Szz} of fairing segments
* angl: row vector of n deterministic parameters (please arrange in increasing order of magnitude)
    * angl(1)   : smallest wing root pitch angle (deg) (e.g., angl(1)=0)
    ...
    * angl(n)   : largest wing root pitch angle (deg) (e.g., angl(n)=60)
--------------------------------------------------------
Output:
* Y: 4xn matrix of quantities of interest (QIs)
    * Y(1, :)      :modal freq. for the first out-of-plane bending mode (Hz)
    * Y(2, :)      :modal freq. for the first in-plane bending mode (Hz)
    * Y(3, :)      :modal freq. for the second out-of-plane bending mode (Hz)
    * Y(4, :)      :modal freq. for the first torsion mode (Hz)
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
[frqs] = G2(run, angl, 'EI', X(1), 'GJ', X(2), 'Sxx', X(3), 'Szz', X(4));

Y = zeros(4, length(angl));
Y(1, :) = frqs(1, :);
Y(2, :) = frqs(2, :);
Y(3, :) = frqs(3, :);
Y(4, :) = frqs(4, :);

end