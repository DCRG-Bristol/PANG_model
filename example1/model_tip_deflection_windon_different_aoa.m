function Y = model_tip_deflection_windon_different_aoa(X, P)
%% Title section - Wind-on tip deflections versus airspeed and E, G
%{
--------------------------------------------------------
Comments:
* The model computes the wind-on tip deflection as a function of uncertain variables (i.e., E, G)
* We use a UQLab style notation, i.e., X: vector of uncertain variables; P: vector of deterministic parameters; Y: vector of model outputs
--------------------------------------------------------
Input:
* X: uncertain variables
    * X(1)      : Young's modulus E for OOP bending (Pa)
    * X(2)      : Young's modulus E for IP bending (Pa)
    * X(3)      : torsional modulus G    
* P: deterministic parameters
    * P(1)      : airspeed (m/s)
    * P(2)      : angle of attack (deg) 
--------------------------------------------------------
Output:
* Y: vector of quantities of interest (QIs)
    * Y(1)      : tip LE displacement in z direction (OOP)
--------------------------------------------------------
%}
%%
%written: SANUJA JAYATILAKE, sj17967@bristol.ac.uk, 13th MAY 2025.

%this is an example of calling the analysis module with embedded aerodynamics model.
%this was created during the 'example_build.m', which created and saved the
% model-specific analysis model (variable 'run_ONERA') in 'run_ONERA.mat';

%load and assign the saved module to the variable 'run'
load('run_ONERA.mat'); run = run_ONERA;

%ASSIGN BASELINE VALUES TO THE PARAMETERS
run = run.setPars('EI_1', X(1), 'EI_2', X(2), 'G', X(3), 'alpha0', 0*pi/180, 'g', 9.81, 'mach', 0.78);

%this sets the baseline parameters - this has to be done for all user
%defined parameters (defined in buildSystem.buildBase.par) -in this case
% there is only 'E', the youn's modulus.

% The embedded parameters will have default values but can be set. These are
% >'U': the airspeed, 
% >alpha': the angle (attitude) between the far filed flow and the O(x,y) plane in the beam model
% >'alpha0' smalled angle between the horizontal and the O(x,y) plane (for gravity)
% >'g': gravitational constant, 9.81 m/s^2
% >'rho'; air density
%> 'mach', mach number

% The parameters set through setPars will be the default/baseline values -
% but can be appended as the model is called - will be shown later

%_______________________________________________________________________________________________________________________________

%SET TRANSFIRMATION FOR STRUCTURAL COORDINATES: set.Transform(type, ORDER)

%the structural degrees of freedom can either be tracked using
% > modal coordinates (type = 'modal') or,
% > using the shape-functional coordinates (type='default'), with the first ORDER (integer) modes.

run = run.setTransform('default'); %....this askes for modelling with functional coorinates 

%_______________________________________________________________________________________________________________________________

% SET STRCTURAL DAMPING MATRIX

%this is to be ssigned to the dampMatr property of the inherited analysis
%module (in this case run).
% 
% In this example, we use the ralaigh-damping model, where the damping
% ratio is a weighed sum of the linear stiffness and mass matrices...

%to do this we first need to get the linearised system matrices: we use the
%method 'getStructModes(q)' whihc computes the strutural system's linear properties 

q0 = run.q0_struct; %this releases a zero-state vector with the correct size - for a structure-only problem (i.e no aerodynamics..)
[shp, evals, Kmat, Cmat, Mmat] = run.getStructModes(q0); 

%in the above, the 'q0_struct' methor releases a zero-state vector
%correct size for a structure-only problem. The 'getStructModes' method
%releases the eigenvector 'shp', the eigenvalues 'evals', and the linear
%structural matrices. Note the default dmaping matrix is zero. 

%NOTE: 'getStructModes' releases eigenvectors based on the asigned
%transformation (see previous step). The damping matrix always needs to be
%in the shape-functional coordinates when assigned, hence 'getStructmodes'
%was called here with type='default' for  set.Transform(type);

%we now synthesise a damping matrix 'Dmat' and assign to teh damping
%property
Dmat = (1e-5)*Mmat + (1e-5)*Kmat;
run.dampMatr = Dmat; 

%_________________________________________________________________________________________________________
%SOME OTHER BITS TO SUPPORT THE PROCEEDING WORK...

%now that the damping model is set, we can now switch to work with the modal
%coordinates, 
run = run.setTransform('modal', 6); %...using the first 6 modes

%NOTE: when setting the above, it uses the eigenmode basis of the system
%with the baseline physical propertues set through 'setPars()'

%and it'll be useful to have zero state vectors of correct sizes to access for structure only and aero-structural 
%problem...

q0 = run.q0_struct; %structure-only problems...
qAero0 = run.q0_aeroStruct; %aeroelastic problems (includes embedded aerodynamic states)

%the below is for plotting work.. we get a mesh of beam-lengthwise
%coordinates. Note the 'geom' class is brought forward as a property in the
%analysis modules, allowing you to access the beam length L.

xMesh = linspace(0,run.geom.L,50);

%__________________________________________________________________________________________________________

%% ANALYSIS

%there're thre fundamental methods to use here...

% > 'qdt = structDeriv(q, 'parName1', val1, 'parName2', val2,...)' 
% computes the time-derivetives of the structural states (FOR A STRUCRURE
% ONLY/WIND_OFF PROBPLEM), with the state vector q (of size similar to q0_struct). You can
%append the parameters from the default/baseline values as you call this. 


% > 'qdt = aero_structDeriv(q, 'parName1', val1, 'parName2', val2,...)' 
% computes the time-derivetives of the aeroelastic states (FOR AN AEROELASTIC PROBLEM),
% with the state vector q (of size similar to q0_aeroStruct). You can
%append the parameters from the default/baseline values as you call this. 

% > [x, y, z] = run.getDisplField(q, xMesh ,frame)'
%gets the displacement field of the wing along arc-lethwise station along
%the beam in xMesh. 
% 
% x = [x_LeadingEdge, x_TrailingEdge], with as many rows as the length of
% xMesh (and similar for y and z outouts).
%
% frame = 'aircraft' releases displacemnets in the aircraft frame (standard
% flight mechanics frame), frame='beamModel' in the O(x,y,z) frame in the
% beam model (difference between the two beeing sweep!)

%% EXAMPLE STATIC SOLUTION - AEROELASTIC PROBLEM>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%static aeroelastic solutions under varying airspeed U, the a time-marched
%simulation of a gust about a computed aeroelastic equilibria

q = qAero0; %initieslised zero vector for the aeroelastic problem..

pitch = P(2)*pi/180; %pitch setting..

%we do not want to keep changing these in the loop so update the default
% value to match the fligth condition of interest..
run = run.setPars('alpha', pitch, 'alpha0', pitch);

q = fsolve(@(q_all)(run.aero_structDeriv(q_all, 'U', P(1))), q); %note E is not set here... the same basline value assigned in line 8 is used
[x, y, z] = run.getDisplField(q, xMesh(end) ,'beamModel'); %get displacements at teh beam tip.. in the beam model frame
Y = z(1,:); % deflection at the leaidng edge

end