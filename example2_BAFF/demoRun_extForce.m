function demoRun_extForce

clear all; close all

%written: SANUJA JAYATILAKE, sj17967@bristol.ac.uk, 19th MAY 2025.

%this is an example of calling the analysis module with an external aerodynamics model.
%this was created during the 'example_build.m', which created and saved the
% model-specific analysis model (variable 'run_ext_a320') in 'run_ONERA.mat';

%load and assign the saved module to the variable 'run'
load('run_ext_a320.mat'); run = run_ext;

%% INITIALISATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%ASSIGN BASELINE VALUES TO THE PARAMETERS
%the parameters to set comes from a set of default parameters that always
%exist by default = any others that were user-specified when buidling the
%model..

run = run.setPars('alpha0', 0*pi/180, 'g', 9.81, 'fuel_frac', 0.72);

%this sets the baseline parameters - this has to be done for all user
%specified  parameters (defined in buildSystem.buildBase.par) -in this case
% there is only 'fuel_frac', the fuel fraction (embedded by default if a BAFF model is used!).

% The embedded parameters will have default values but can be set. These are
% >'U': the airspeed, 
% >alpha': the angle (attitude) between the far filed flow and the O(x,y) plane in the beam model
% >'alpha0' smallest angle between the horizontal and the O(x,y) plane (for gravity)
% >'g': gravitational constant, 9.81 m/s^2
% >'rho'; air density

%Note: when using an external aerodynamic model, 'rho', 'U', 'alpha' makes
%no difference!

% The parameters set through setPars will be the default/baseline values -
% but can be appended as the model is called - will be shown later

%_______________________________________________________________________________________________________________________________

%SET TRANSFIRMATION FOR STRUCTURAL COORDINATES: 
% setTransform(type, ORDER)

%the structural degrees of freedom can either be tracked using
% > modal coordinates (type = 'modal') with the first ORDER (integer) modes, or,
% > using the shape-functional coordinates (type='default'), .

run = run.setTransform('default'); %....this asks for modelling with functional coordinates 

%_______________________________________________________________________________________________________________________________

% SET STRCTURAL DAMPING MATRIX

%this is to be assigned to the dampMatr property of the inherited analysis
%module ('run').
% 
% In this example, we use the raleigh-damping model, where the damping
% ratio is a weighed sum of the linear stiffness and mass matrices...

%to do this we first need to get the linearised mass and stiffness matrices: use the
%method 'getStructModes(q)' which computes the structural system's linear properties 

q0 = run.q0_struct; %this releases a zero-state vector with the correct size - for a structure-only problem (i.e no aerodynamics..)
[shp, evals, Kmat, Cmat, Mmat] = run.getStructModes(q0); 

%in the above, the 'q0_struct' method releases a zero-state vector
%correct size for a structure-only problem. The 'getStructModes' method
%releases the eigenvector 'shp', the eigenvalues 'evals', and the linear
%structural matrices. Note the default damping matrix ('Cmat') is zero. 

%NOTE: 'getStructModes' releases eigenvectors based on the asigned
%transformation (see previous step). The damping matrix always needs to be
%in the shape-functional coordinates when assigned, hence 'getStructmodes'
%was called here with type='default' for  setTransform(type);

% now synthesise a damping matrix 'Dmat' and assign to teh damping
%property
Dmat = (1e-5)*Mmat + (1e-5)*Kmat;
run.dampMatr = Dmat; 

%_________________________________________________________________________________________________________
%SOME OTHER BITS TO SUPPORT THE PROCEEDING WORK...

%now that the damping moel is set, we can now switch to work with the modal
%coordinates, 
run = run.setTransform('modal', 10); %...using the first 10 modes

%NOTE: when setting the above, it uses the eigenmode basis of the system
%with the baseline physical propertues set through 'setPars()'

%and it'll be useful to have zero state vectors of correct sizes to access for structure only and aero-structural 
%problems...note this will be updated whenever a coordinate transformation is set

q0 = run.q0_struct; %structure-only problems...
qAero0 = run.q0_aeroStruct; %aeroelastic problems (includes embedded aerodynamic states)

%note in this specific instance (working with external aero models),
%q0_aeroStruct is the same as q0_struct (no additional aero sttes are
%handled in the PANG model)

%you can check the wind-off mode shapes with the method below...
run.plotStructModes(q0);
%__________________________________________________________________________________________________________

%% ANALYSIS

%there're three fundamental methods to use here...

% > 'qdt = structDeriv(q, 'parName1', val1, 'parName2', val2,...)' 
% releases the time-derivetives of the structural states (FOR A STRUCRURE
% ONLY/WIND_OFF PROBPLEM), with the state vector q (of size similar to q0_struct). You can
%append the parameters from the default/baseline values as you call this. 

% > 'qdt = aero_structDeriv(q, 'parName1', L, D, M, ang, val1, 'parName2', val2,...)' 
% computes the time-derivetives of the aeroelastic states (FOR AN AEROELASTIC PROBLEM),
% with the state vector q (of size similar to q0_aeroStruct). 
% > L, D and M are coloumn vectors of spanwise follower lift, drag and
% moments for each panel to be applied at the collocation stations along the beam 
% -the collocation points are at the mid-points of the aerodynamic grid
% defined in basis.xi (thus, there are length(basis.xi)-1 collocation
% stations) during model building. Lift and drag are applied at the
% respective quarter chord positions. 
% The applied loads are rotated by 'ang' - a column vector of length length(basis.xi)-1 
%you can append the parameters from the default/baseline values as you call
%this. (NOTE: this rotation is performed on top of the twist). If lift and
%drag are calculated perpendicular and parallel to the Cl=0 line in the
%aerofoil the ang can be a vector of zeros.

% > [x, y, z] = run.getDisplField(q, xMesh ,frame)
%gets the LE/TE displacement field of the wing along arc-lengthwise
%stations of the elastic axis of the beam in xMesh. 
% 
% x = [x_LeadingEdge; x_TrailingEdge], with as many columns as the length of
% xMesh (and similar for y and z outouts).
%
% frame = 'aircraft' releases displacements in the aircraft frame (standard
% flight mechanics frame), frame='beamModel' in the O(x,y,z) frame in the
% beam model (difference between the two being sweep and teh direction of the z axis)

% > [x, y, z] = run.getMesh(q, type)
%same as getDisplField, but only computes displacemnets at the aerodynamic
%mesh on the beam as defined in basis.xi in 'example_build.m'. Always use
%this to extract the mesh for computing aerodynamic loads!

% > [x, y, z] = run.getColloc(q, type)
%this shoots out the collocation points used to apply the loads...basically
%the quarter chord point in each panel...

%% apply increasing lift forces....

%this is an example of applying a uniformly dirtibuted lift on the wing (no
%moments and drag) in increments as the equilibria are calculated..
%note: aero_structDeriv gives the veclocities of the states so can be used
%for time marching, the example below only shows a static solution
%sequence.

Lfctr = [0,  2.5,  5]*10000; %increments lift (equally distributed) to apply

figure; subplot(1,1,1); 
clr = {'k', 'b', 'r'}; %some colors to plot the step incrmenets

%initialise loading vectors and the rotation vectors...
[L0, D0, M0] = run.initLoads; %this releases a set of 'ones vectors' of the correct size for the lift, drag and moments
ang = ones(size(L0))*0; %ang has the same size (i.e. number of collocation stations)
q = qAero0; %set of zero state vectors of the correct size...

for i=1:length(Lfctr)
    q = fsolve(@(q_all)(run.aero_structDeriv(q_all, L0*Lfctr(i), D0*0, M0*0, ang)), q); %call the time-derivative function as the function to solve for the static equilibria

    %plot LE and TE of the aerodynamic mesh printed...
    [x, y, z] = run.getMesh(q, 'aircraft');
    plot3(x(1,:), y(1,:), z(1,:), clr{i}, 'marker', 'x'); hold on;
    plot3(x(2,:), y(2,:), z(2,:), clr{i}, 'marker', 'x'); hold on;

    %plot the collocation points at which the loads are applied...
    [xQtr, yQtr, zQtr] = run.getColloc(q, 'aircraft');
    plot3(xQtr, yQtr, zQtr, 'o--', 'color', clr{i}, 'linewidth', 1,...
        'markerSize', 4); hold on;
end

set(gca, 'dataAspectRatio', [1,1,1], 'zDir', 'reverse'); drawnow;
xlabel('x, [m]'); ylabel('y, [m]'); zlabel('z, [m]');

% plot the last load case to show the panels with more detail...............
panNum = length(run.basis.xi)-1;

figure;
for pan=1:panNum
    xPan = [x(1,pan), x(2,pan), x(2,pan+1), x(1,pan+1), x(1,pan)]; %a trapeziod with LE-TE-TE-LE-LE points...
    yPan = [y(1,pan), y(2,pan), y(2,pan+1), y(1,pan+1), y(1,pan)]; 
    zPan = [z(1,pan), z(2,pan), z(2,pan+1), z(1,pan+1), z(1,pan)]; 

    plot3(xPan, yPan, zPan, 'x-', 'Color',clr{end}, 'linewidth', 1); hold on; %plot the panel
end
plot3(xQtr, yQtr, zQtr, 'o--', 'color', clr{end}, 'linewidth', 1,...
    'markerSize', 4); hold on;
set(gca, 'dataaspectratio', [1,1,1], 'ZDir', 'reverse');
grid minor;
xlabel('X'); ylabel('Y'); zlabel('Z'); grid minor;


end