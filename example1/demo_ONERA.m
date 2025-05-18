function demo_ONERA

clear all; close all;

%written: SANUJA JAYATILAKE, sj17967@bristol.ac.uk, 13th MAY 2025.

%this is an example of calling the analysis module with embedded aerodynamics model.
%this was created during the 'example_build.m', which created and saved the
% model-specific analysis model (variable 'run_ONERA') in 'run_ONERA.mat';

%load and assign the saved module to the variable 'run'
load('run_ONERA.mat'); run = run_ONERA;

%% INITIALISATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%ASSIGN BASELINE VALUES TO THE PARAMETERS
run = run.setPars('E', 70e9, 'alpha0', 0*pi/180, 'g', 9.81, 'mach', 0.78);

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

%% EXAMPLE STATIC SOLUTION - AERO OFF PROBLEM>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%this is an example involving rottaing the wing root pitch/attitude between
%0 and 90 degrees from the horizontal and get static structural
%deflections. This is demonstarted for three values of Young's moduli 'E'
%(uuser assined parameter - see example_build.m')

q = q0; %zero state vector for structure only problem..
alp = linspace(pi/2,0,15); %ranges of pitch angles...
Er = [60, 70, 80]*1e9; %values of youn's moduli
clr = {'b', 'r', 'k'}; %some colors for plotting 

figure('Name', 'static structural solution sequence');

subplot(1,2,1); xlabel('Pitch angle, [deg]'); ylabel('Tip Deflection, [m]'); hold on;
subplot(1,2,2); xlabel('Pitch angle, [deg]'); ylabel('Modal frequencies, [Hz]'); hold on;

for k_idx = 1:length(Er) %loop through young's moduli
    par_label{k_idx} = ['E = ',num2str(Er(k_idx)/(1e9)), 'GPa'];

    for a_idx = 1:length(alp) %loop through pitch angles

        q = fsolve(@(qstr)(run.structDeriv(qstr,'alpha0', alp(a_idx), 'E', Er(k_idx))), q); 

        [x, y, z] = run.getDisplField(q, xMesh(end) ,'beamModel'); %get displacements at x
        tipDefl_str{k_idx}(a_idx) = z(1,:); %z of size [2, length(xMeasureStat)] hold the z displacements in z(1,:) and TE in z(2,:)

        %eigenvalues....
        [~, e] = run.getStructModes(q, 'E', Er(k_idx), 'alpha0', alp(a_idx)); %this calculates the structural modes with the reference parameters...
        EigVals{k_idx}(:,a_idx) = e;
    end

    subplot(1,2,1); plot(alp*180/pi, tipDefl_str{k_idx}, 'x-', 'color', clr{k_idx}); hold on; 
    subplot(1,2,2); 
    for mode=1:4
        plot(alp*180/pi, abs(EigVals{k_idx}(mode,:))./(2*pi), 'x-', 'color', clr{k_idx}); hold on;
    end
    drawnow;
end
subplot(1,2,1); legend(par_label)

%% EXAMPLE STATIC SOLUTION - AEROELASTIC PROBLEM>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%static aeroelastic solutions under varying airspeed U, the a time-marched
%simulation of a gust about a computed aeroelastic equilibria

Ur = linspace(50,175,10); %airspeed range
q = qAero0; %initieslised zero vector for the aeroelastic problem..

pitch = 5*pi/180; %pitch setting..

%we do not want to keep changing these in the loop so update the default
% value to match the fligth condition of interest..
run = run.setPars('alpha', pitch, 'alpha0', pitch);

for u_idx=1:length(Ur)
    U = Ur(u_idx);
    q = fsolve(@(q_all)(run.aero_structDeriv(q_all, 'U', U)), q); %note E is not set here... the same basline value assigned in line 8 is used

    [x, y, z] = run.getDisplField(q, xMesh(end) ,'beamModel'); %get displacements at teh beam tip.. in the beam model frame
    tipDefl(u_idx) = z(1,:); % deflection at the leaidng edge
end

figure; plot(Ur,tipDefl, 'x-'); xlabel('Airspeed, [m/s]'); ylabel('tip LE vertical deflection, [m]');

%a comparison of the displacemnets in the beam and aircaft frames...using
%the last computed equilirbrium

figure('Name', 'Static aeroelastic equilibria'); 

[x, y, z] = run.getDisplField(q, xMesh ,'beamModel'); %get displacements at x
X = [x(1,:), flip(x(2,:))];
Y = [y(1,:), flip(y(2,:))];
Z = [z(1,:), flip(z(2,:))];
subplot(1,2,1)
plot3(X, Y, Z, 'x-', 'linewidth', 2); title('beam frame')
xlabel('x, [m]'); ylabel('y, [m]'); zlabel('z, [m]');
set(gca, 'dataaspectratio', [1,1,1])

[x, y, z] = run.getDisplField(q, xMesh ,'aircraft'); %get displacements at x
X = [x(1,:), flip(x(2,:))];
Y = [y(1,:), flip(y(2,:))];
Z = [z(1,:), flip(z(2,:))];
subplot(1,2,2)
plot3(X, Y, Z, 'x-', 'linewidth', 2); title('aircraft frame')
set(gca, 'dataaspectratio', [1,1,1])
set(gca, 'zDir', 'reverse');
xlabel('x, [m]'); ylabel('y, [m]'); zlabel('z, [m]');

% ODE SOLUTION SEQUENCE OF A GUST>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

gustSeq_aoa = @(t)(0.5*(1-cos(2*pi*t/0.5)).*heaviside(0.5-t)); %gust input..

%run simulation....'alpha' is appended with the gust function as the
%time-derivatve function is called...

[T,Q] = ode45(@(t,x)(run.aero_structDeriv(x, 'U', U, 'alpha', pitch + (5*pi/180)*gustSeq_aoa(t))),...
    linspace(0,1,100), q);

%animate the response...........................................
figure; 
[x, y, z] = run.getDisplField(q, xMesh ,'aircraft'); %get displacements at x
X = [x(1,:), flip(x(2,:))]; Y = [y(1,:), flip(y(2,:))]; Z = [z(1,:), flip(z(2,:))];
plot3(X, Y, Z, '-', 'linewidth', 1, 'color', 'k'); hold on;
animObj = plot3(X, Y, Z, '-', 'linewidth', 2, 'color', 'r');
set(gca,'zDir', 'reverse', 'dataaspectratio', [1,1,1]); grid minor;
view([-120,20])
ylim([0, 1.1*run.geom.L]);
zlim([-0.3, 0.3]*run.geom.L);

drawnow;

set(gca, 'dataaspectratio', [1,1,1])
for t_idx = 2:length(T)
    [x, y, z] = run.getDisplField(Q(t_idx,:)', xMesh ,'aircraft'); %get displacements at x
    X = [x(1,:), flip(x(2,:))]; Y = [y(1,:), flip(y(2,:))]; Z = [z(1,:), flip(z(2,:))];
    pause(0.05)
    set(animObj, 'xdata', X, 'ydata', Y, 'zdata', Z);  

    tipDisp_sim(t_idx-1) = z(1,end);
end

figure('name', 'Tip displacement repsosne....');
plot(T(2:end), tipDisp_sim/run.geom.L, 'k-', 'linewidth', 1);
xlabel('t, [s]'); ylabel('LE tip deflection, [-]');

end