function demoRun_extForce

clear all; close all

load('run_ext.mat');
run = run_ext;

run = run.setPars('E', 70e9, 'alpha0', 0*pi/180, 'g', 9.81); %this sets the baseline parameters - these numbers will be used if properties are not specified when calling further functions

run = run.setTransform('modal', 6); %the model maybe called either with functional coordinates or with model coordinates, this askes for modal coordinates with the first 12 modes
run = run.setTransform('default'); %....this askes for modelling with functional coorinates 

%lets work with the default transformation for now..

q0 = run.q0_struct; %this releases a zero-state vector with the correct size - for a structure-only problem (i.e no aerodynamics..)
[shp, evals, Kmat, Cmat, Mmat] = run.getStructModes(q0); %this does an eigenvalue analysis (structure only) and returns linearired matrices.. note here its done with the baseline parameters..

%we can create a damping matrix, not Cmat above is just zero by default. 
Dmat = (1e-5)*Mmat + (1e-5)*Kmat;
run.dampMatr = Dmat; %this creates a dmping matrix... note this has to be for the functional coorinates..the 'getStructModes' done before was first made with the defalt transformation

%now lets work with the modal basis, only using the first siz modes...
run = run.setTransform('modal', 6);

%% apply increasing lift forces....

Lfctr = [0, 2,  5,  7.5, 10]*10000;
figure; subplot(1,1,1); 
clr = {'r', 'b', 'k', 'r', 'b', 'k'};

[L0, D0, M0] = run.initLoads; %this releases a set of `ones' vectors of the correct size for the lift drag and moments (per unit span!)

ang = ones(run.basis.Ngam,1)*0;
q = run.q0_aeroStruct;

for i=1:length(Lfctr)
    q = fsolve(@(q_all)(run.aero_structDeriv(q_all, L0*Lfctr(i), D0, M0, ang)), q);
    %the arguments after q_all are respectivley follow lift drag and momeents at the collocation stations. the lift and drag can be rotated by 'ang' if desired

    tic
    [x, y, z] = run.getMesh(q, 'aircraft');
    toc
    plot3(x(1,:), y(1,:), z(1,:), clr{i}, 'marker', 'x'); hold on;
    plot3(x(2,:), y(2,:), z(2,:), clr{i}, 'marker', 'x'); hold on;
    set(gca, 'dataAspectRatio', [1,1,1]); drawnow;
end

end