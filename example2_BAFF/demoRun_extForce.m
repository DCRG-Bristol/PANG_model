function demoRun_extForce

clear all; close all

%written: SANUJA JAYATILAKE, sj17967@bristol.ac.uk, 13th MAY 2025.

%this is an example of calling the analysis module with embedded aerodynamics model.
%this was created during the 'example_build.m', which created and saved the
% model-specific analysis model (variable 'run_ONERA') in 'run_ONERA.mat';

%load and assign the saved module to the variable 'run'
load('run_ext_a320.mat'); run = run_ext;

%% INITIALISATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%ASSIGN BASELINE VALUES TO THE PARAMETERS
run = run.setPars('alpha0', 0*pi/180, 'g', 9.81, 'fuel_frac', 1);

%this sets the baseline parameters - this has to be done for all user
%defined parameters (defined in buildSystem.buildBase.par) -in this case
% there is only 'E', the youn's modulus.

% The embedded parameters will have default values but can be set. These are
% >'U': the airspeed, 
% >alpha': the angle (attitude) between the far filed flow and the O(x,y) plane in the beam model
% >'alpha0' smalled angle between the horizontal and the O(x,y) plane (for gravity)
% >'g': gravitational constant, 9.81 m/s^2
% >'rho'; air density

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

%now that the damping moel is set, we can now switch to work with the modal
%coordinates, 
run = run.setTransform('modal', 6); %...using the first 6 modes

%NOTE: when setting the above, it uses the eigenmode basis of the system
%with the baseline physical propertues set through 'setPars()'

%and it'll be useful to have zero state vectors of correct sizes to access for structure only and aero-structural 
%problem...

q0 = run.q0_struct; %structure-only problems...
qAero0 = run.q0_aeroStruct; %aeroelastic problems (includes embedded aerodynamic states)

%__________________________________________________________________________________________________________

%% ANALYSIS

%there're thre fundamental methods to use here...

% > 'qdt = structDeriv(q, 'parName1', val1, 'parName2', val2,...)' 
% computes the time-derivetives of the structural states (FOR A STRUCRURE
% ONLY/WIND_OFF PROBPLEM), with the state vector q (of size similar to q0_struct). You can
%append the parameters from the default/baseline values as you call this. 


% > 'qdt = aero_structDeriv(q, 'parName1', L, D, M, ang, val1, 'parName2', val2,...)' 
% computes the time-derivetives of the aeroelastic states (FOR AN AEROELASTIC PROBLEM),
% with the state vector q (of size similar to q0_aeroStruct). 
% > L, D and M are coloumn vectors of spanwise follower lift drag and moments (PER UNIT SPAN)
% to be applied at the collocation stations along the arc-lngth of the beam 
% -the collocation points are at the mid-points of the aerodynamic grid
% defined in basis.xi (NOTE: there are length(basis.xi)-1 collocation
% stations) during model building. lift and drag are applied at the
% respectve quater chord positions. 
% The applied loads are rotated by 'ang' - a column vector of length length(basis.xi)-1 
%you can append the parameters from the default/baseline values as you call this.

%for context, the collocation stations are at:
xColloc = 0.5*(run.basis.xi(1:end-1) + run.basis.xi(2:end));

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

% > [x, y, z] = run.getMesh(q, type)
%same as getDisplField, but only computes displacemnets at the aerodynamic
%mesh on the beam as defined in basis.xi in 'example_build.m'


%% apply increasing lift forces....

Lfctr = [0,  5,  10]*10000; %lift per unit span (equally distributed) to apply
figure; subplot(1,1,1); 
clr = {'c', 'b', 'r'};

[L0, D0, M0] = run.initLoads; %this releases a set of `ones' vectors of the correct size for the lift drag and moments (per unit span!)
ang = ones(size(L0))*0;
q = run.q0_aeroStruct;

for i=1:length(Lfctr)
    q = fsolve(@(q_all)(run.aero_structDeriv(q_all, L0*Lfctr(i), D0, M0, ang)), q);

    %plot LE and TE of the aerodynamic mesh printed...
    [x, y, z] = run.getMesh(q, 'aircraft');
    plot3(x(1,:), y(1,:), z(1,:), clr{i}, 'marker', 'x'); hold on;
    plot3(x(2,:), y(2,:), z(2,:), clr{i}, 'marker', 'x'); hold on;

    %plot the quater-chord points of the collocation points where the forces
    %are applied..
    [xC, yC, zC] = run.getDisplField(q, xColloc ,'aircraft');
    xQtr = 0.75*xC(1,:) + 0.25*xC(2,:); 
    yQtr = 0.75*yC(1,:) + 0.25*yC(2,:); 
    zQtr = 0.75*zC(1,:) + 0.25*zC(2,:); 
    plot3(xQtr, yQtr, zQtr, 'o--', 'color', 'k', 'linewidth', 1); hold on;
end

set(gca, 'dataAspectRatio', [1,1,1], 'zDir', 'reverse'); drawnow;
xlabel('x, [m]'); ylabel('y, [m]'); zlabel('z, [m]');

end