%% example for building the model...

%note: 'tbx' fodler must be added to the path

clear all; close all;

buildOb = buildSystem.buildBase; %call bass class for model building...
buildOb.fileLocation = fileparts(mfilename('fullpath')); %set the file location where a folder with a model-specific items could be created...
buildOb.name = 'SJD_Wing'; %a name for the model

%% some example 

L = 0.655; %semi-span

%% planform geometry...

geom = buildSystem.geom;
geom.L = L;
geom.a = @(x)(0*ones(size(x))); %position of the elastic axis... theodorsen definition
geom.b = @(x)(0.075*ones(size(x))); %semi-chord distribution
geom.sweep = 0*pi/180; %sweep (of beam line!)
geom.twist = @(x)(0*ones(size(x))); %wing twist distribution
geom.WTWidth = 0.8;

buildOb.geom = geom; %write to buildBass class...

%% basis properties

basis = buildSystem.basis;
basis.Nw = 6; %size of shape basis..out of plane bending
basis.Nv = 4; %.. chord wise bending
basis.Nthet = 4;%...torsion
basis.xi = linspace(0.005, L, 21); %aerodynamic grid (spanwise)... collocation points assigned at the mid-points of xi

buildOb.basis = basis; %write to buildBass class...

%% elastic properties

% this is an example of a case where we want to change the Young's modulus
% when calling the buit model as a parameter... note that the rigidities
% EI1, etc depends in this..

elas(1) = buildSystem.structure.elasBase; %call elastic property class

%stiffness properties... matrices with these rigidities will be multiplied
%by a user dfined parameters...
elas(1).EI1 = @(x)(2.268*ones(size(x)));
elas(1).EI2 = @(x)(59.136*ones(size(x)));
elas(1).EI12 = @(x)(0);
elas(1).GJ = @(x)(0); %GJ handled seperately as we want this to be independednt of E
elas(1).name = 'K_E'; %namee for this matrix..
elas(1).fctrId = 'EI'; %this property, which must be one of that defined among user parameters, will scale this matrix computed for elas(1)

elas(2).EI1 = @(x)(0);
elas(2).EI2 = @(x)(0);
elas(2).EI12 = @(x)(0);
elas(2).GJ = @(x)(3.3843); %GJ handled seperately as we want this to be independednt of E
elas(2).name = 'K_G'; %namee for this matrix..
elas(2).fctrId = 'GJ'; %this property, which must be one of that defined among user parameters, will scale this matrix computed for elas(1)

buildOb.elas = elas; %write to buildBass class...

%% add inertia

inertia(1) = buildSystem.structure.inertiaBase; %inertia property class
inertia(1).e = @(x)(0); %position of mass axis.. theodorsen 
inertia(1).m = @(x)(0.1242); %mass per unit length [kg/m]
inertia(1).mxx = @(x)(2.7400e-06); %twisting inertia about chordwise about the CoM [kgm]
inertia(1).mzz = @(x)(2.6500e-06); %chord-plane rotation inertia [kgm]
inertia(1).myy = @(x)(0); %rotation inertia..bending rotation [kgm]

for i=1:10
    elem(i) = buildSystem.structure.descrElem; %descrete inertia for the engine....handled as a point mass
    %all properties for this follow the same definitions as that for
    %inertia...but NOT PER UNIT LENGTH AND ARE NOT FUNCTION HANDLES
    elem(i).m = 0.0613;
    elem(i).mxx = 0;
    elem(i).mzz = 0;
    elem(i).myy = 1.8399e-05;
    elem(i).e = 0.05;
    elem(i).xp = 0.0375 + (i-1)*0.0650; %attachment potision...say 30% of the semi-span
end

%guide-masses
xpLoc = 0.0375 + ([2,4,6,8]-1)*0.0650;
for i=11:14
    elem(i) = buildSystem.structure.descrElem; %descrete inertia for the engine....handled as a point mass
    %all properties for this follow the same definitions as that for
    %inertia...but NOT PER UNIT LENGTH AND ARE NOT FUNCTION HANDLES
    elem(i).m = 0.01;
    elem(i).mxx = 5.7704e-06;
    elem(i).mzz = 5.6943e-06;
    elem(i).myy = 2.4889e-07;
    elem(i).e = -0.4;
    elem(i).xp = xpLoc(i-10); %attachment potision...say 30% of the semi-span
end

inertia(1).elem = elem;
inertia(1).name = 'massModel';
inertia(1).fctrId = []; %no scalling is imposed on the mass matrix.. however can be implemented if needed

%%

inertia(2) = buildSystem.structure.inertiaBase; %inertia property class
inertia(2).e = @(x)(0); %position of mass axis.. theodorsen 
inertia(2).m = @(x)(0); %mass per unit length [kg/m]
inertia(2).mxx = @(x)(0); %twisting inertia about chordwise about the CoM [kgm]
inertia(2).mzz = @(x)(0); %chord-plane rotation inertia [kgm]
inertia(2).myy = @(x)(0); %rotation inertia..bending rotation [kgm]

for i=1:10
    elem_xx(i) = buildSystem.structure.descrElem; %descrete inertia for the engine....handled as a point mass
    %all properties for this follow the same definitions as that for
    %inertia...but NOT PER UNIT LENGTH AND ARE NOT FUNCTION HANDLES
    elem_xx(i).m = 0;
    elem_xx(i).mxx = 1.1796e-04;
    elem_xx(i).mzz = 0;
    elem_xx(i).myy = 0;
    elem_xx(i).e = 0.05;
    elem_xx(i).xp = 0.0375 + (i-1)*0.0650; %attachment potision...say 30% of the semi-span
end

inertia(2).elem = elem_xx;
inertia(2).name = 'SxxMatr';
inertia(2).fctrId = 'Sxx'; %no scalling is imposed on the mass matrix.. however can be implemented if needed

%%

inertia(3) = buildSystem.structure.inertiaBase; %inertia property class
inertia(3).e = @(x)(0); %position of mass axis.. theodorsen 
inertia(3).m = @(x)(0); %mass per unit length [kg/m]
inertia(3).mxx = @(x)(0); %twisting inertia about chordwise about the CoM [kgm]
inertia(3).mzz = @(x)(0); %chord-plane rotation inertia [kgm]
inertia(3).myy = @(x)(0); %rotation inertia..bending rotation [kgm]

for i=1:10
    elem_zz(i) = buildSystem.structure.descrElem; %descrete inertia for the engine....handled as a point mass
    %all properties for this follow the same definitions as that for
    %inertia...but NOT PER UNIT LENGTH AND ARE NOT FUNCTION HANDLES
    elem_zz(i).m = 0;
    elem_zz(i).mxx = 0;
    elem_zz(i).mzz = 1.3197e-04;
    elem_zz(i).myy = 0;
    elem_zz(i).e = 0.05;
    elem_zz(i).xp = 0.0375 + (i-1)*0.0650; %attachment potision...say 30% of the semi-span
end

inertia(3).elem = elem_zz;
inertia(3).name = 'SzzMatr';
inertia(3).fctrId = 'Szz'; %no scalling is imposed on the mass matrix.. however can be implemented if needed


buildOb.inertia = inertia;

%% gravity..

grav = buildOb.inertia2grav(inertia(1)); % this method in the buildBass class automatically generates a class with gravitational/weigth properties using the assigned inertia properties
buildOb.grav = grav; %write to buildBass class...

%% create the model...

buildOb.prepFolder; %this method sets up the folder and the necesary sub-content to 
buildOb.writeModel;

%% derive analysis models... the buildBass instance is passed in 

run_ONERA = analysis.oneraBase(buildOb); %analysis module foe the embedded-aerodynamic model
run_ext = analysis.extAeroBase(buildOb); %analysis module foe the external aerodynamic model

%add problem-specific aerodynamic data....
run_ONERA.Cl_grad = @(alp, U)(aeroCurves.clGradFcn(alp,U));
run_ONERA.Cl = @(alp, U)(aeroCurves.clFcn(alp,U));
run_ONERA.Cm = @(alp, U)(aeroCurves.cmFcn(alp,U));
run_ONERA.Cd = @(alp, U)(aeroCurves.cdFcn(alp,U));
run_ONERA.lossFactor = 0.05;
run_ONERA.ML = 0.44;
run_ONERA.lam = 0.275;

%%
run_ONERA = run_ONERA.setPars('EI', 1, 'GJ', 1, 'Sxx', 1, 'Szz', 1);
run_ext = run_ext.setPars('EI', 1, 'GJ', 1, 'Sxx', 1, 'Szz', 1);

%set damping matrices...
qstr0 = fsolve(@(q_all)(run_ONERA.structDeriv(q_all,'alpha', 0, 'alpha0', 0)),...
    run_ONERA.q0_struct);

%get linear modes and matrices....
[shp, evals, Kmat, Cmat, Mmat] = run_ONERA.getStructModes(qstr0);

%add damping (values from paper)
Dmat = 0.071*Mmat + (10e-5)*Kmat;
run_ONERA.dampMatr = Dmat;

qstr0 = fsolve(@(q_all)(run_ext.structDeriv(q_all,'alpha', 0, 'alpha0', 0)),...
    run_ext.q0_struct);

%get linear modes and matrices....
[shp, evals, Kmat, Cmat, Mmat] = run_ext.getStructModes(qstr0);

%add damping (values from paper)
Dmat = 0.071*Mmat + (10e-5)*Kmat;
run_ext.dampMatr = Dmat;


%save these analysis modules for calling later...
save('run_ONERA.mat', 'run_ONERA')
save('run_ext.mat', 'run_ext')

%%