%% example for building the model...
%note: 'toolbox' fodler must be added to the path

%example parameters for 777-like wing (stipulated!)

clear all; close all;

buildOb = buildSystem.buildBase; %call bass class for model building...
buildOb.fileLocation = fileparts(mfilename('fullpath')); %set the file location where a folder with a model-specific items could be created...
buildOb.name = '777_likeWing'; %a name for the model
buildOb.par = {'E'}; %you can define some parameters (elastic and mass) that can be varied when calling a built model.. they should be defined here...
%% some example 

L = 29.4; %semi-span

%% planform geometry...

geom = buildSystem.geom;
geom.L = L;
geom.a = @(x)(0.8-1); %position of the elastic axis... theodorsen definition
geom.b = @(x)(7 + (1.5-7)*x/L); %semi-chord distribution
geom.sweep = 35*pi/180; %sweep (of beam line!)
geom.twist = @(x)(7*pi/180 + (pi/180)*(-4-7)*x/L); %wing twist distribution

buildOb.geom = geom; %write to buildBass class...

%% basis properties

basis = buildSystem.basis;
basis.Nw = 4; %size of shape basis..out of plane bending
basis.Nv = 4; %.. chord wise bending
basis.Nthet = 4;%...torsion
basis.xi = linspace(0, L, 21); %aerodynamic grid (spanwise)... collocation points assigned at the mid-points of xi

buildOb.basis = basis; %write to buildBass class...

%% elastic properties

% this is an example of a case where we want to change the Young's modulus
% when calling the buit model as a parameter... note that the rigidities
% EI1, etc depends in this..

elas(1) = buildSystem.structure.elasBase; %call elastic property class

%stiffness properties... matrices with these rigidities will be multiplied
%by a user dfined parameters...
elas(1).EI1 = @(x)(100*exp(18.8-7.3*x/L)/(70e9));
elas(1).EI2 = @(x)(5*100*exp(18.8-7.3*x/L)/(70e9));
elas(1).EI12 = @(x)(0);
elas(1).GJ = @(x)(0); %GJ handled seperately as we want this to be independednt of E
elas(1).name = 'K_E'; %namee for this matrix..
elas(1).fctrId = 'E'; %this property, which must be one of that defined among user parameters, will scale this matrix computed for elas(1)

%a stiffnes matrix for torsional rigidity...here we keep it constant (i.e.
%not dependednt on 

elas(2) = buildSystem.structure.elasBase; %call elastic property class
elas(2).EI1 = @(x)(0);
elas(2).EI2 = @(x)(0);
elas(2).EI12 = @(x)(0);
elas(2).GJ = @(x)(100*exp(18.3-6.3*x/L)); %we want to use a fixed GJ model here... 
elas(2).name = 'K_G';
elas(2).fctrId = []; %..hence the fctrId property is left empty...

buildOb.elas = elas; %write to buildBass class...

%% add inertia

inertia = buildSystem.structure.inertiaBase; %inertia property class
inertia.e = @(x)(0.8-1); %position of mass axis.. theodorsen 
inertia.m = @(x)(1e3*(0.7 + (0.025-0.7)*x/L)); %mass per unit length [kg/m]
inertia.mxx = @(x)(0.25*(1e3*(0.7 + (0.025-0.7)*x/L)).*(7 + (1.5-7)*x/L).^2); %twisting inertia about chordwise about the CoM [kgm]
inertia.mzz = @(x)(0.25*(1e3*(0.7 + (0.025-0.7)*x/L)).*(7 + (1.5-7)*x/L).^2); %chord-plane rotation inertia [kgm]
inertia.myy = @(x)(0); %rotation inertia..bending rotation [kgm]

elem(1) = buildSystem.structure.descrElem; %descrete inertia for the engine....handled as a point mass
%all properties for this follow the same definitions as that for
%inertia...but NOT PER UNIT LENGTH AND ARE NOT FUNCTION HANDLES
elem(1).m = 7800;
elem(1).mxx = 0;
elem(1).mzz = 0;
elem(1).myy = 0;
elem(1).e = -0.1;
elem(1).xp = 0.3*L; %attachment potision...say 30% of the semi-span
inertia(1).elem = elem;
inertia(1).name = 'massModel';
inertia(1).fctrId = []; %no scalling is imposed on the mass matrix.. however can be implemented if needed

buildOb.inertia = inertia; %write to buildBass class...

%% gravity..

grav = buildOb.inertia2grav(inertia); % this method in the buildBass class automatically generates a class with gravitational/weigth properties using the assigned inertia properties
buildOb.grav = grav; %write to buildBass class...

%% create the model...

buildOb.prepFolder; %this method sets up the folder and the necesary sub-content to 
buildOb.writeModel;

%% derive analysis models... the buildBass instance is passed in 

run_ONERA = analysis.oneraBase(buildOb); %analysis module foe the embedded-aerodynamic model
run_ext = analysis.extAeroBase(buildOb); %analysis module foe the external aerodynamic model

%save these analysis modules for calling later...
save('run_ONERA.mat', 'run_ONERA')
save('run_ext.mat', 'run_ext')

%%