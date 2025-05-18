%% example for building the model...

clear all;

%load baff object_____________
ADP = load('UB321_baseline_simple.mat');
ADP = ADP(1).ADP.BuildBaff;
Baff_mdl = ADP.Baff;

%%
buildOb = buildSystem.buildBase; %call bass class for model building...
buildOb.fileLocation = fileparts(mfilename('fullpath')); %set the file location where a folder with a model-specific items could be created...
buildOb.name = 'a320_likeWing'; %a name for the model

%% basis properties

basis = buildSystem.basis;
basis.Nw = 4; %size of shape basis..out of plane bending
basis.Nv = 4; %.. chord wise bending
basis.Nthet = 4;%...torsion

buildOb.basis = basis; %write to buildBass class...

buildOb = buildOb.baff2PANG(Baff_mdl, 'Wing_RHS'); %this builds a PANG model definition from a baff
buildOb.basis.xi = linspace(0, buildOb.geom.L, 20); %re assign a grid for aero: beam arc-lengthwise points

%% create the model...

buildOb.prepFolder; %this method sets up the folder and the necesary sub-content to 
buildOb.writeModel; %develop system

%% derive analysis models... the buildBass instance is passed in 

run_ONERA = analysis.oneraBase(buildOb); %analysis module foe the embedded-aerodynamic model
run_ext = analysis.extAeroBase(buildOb); %analysis module foe the external aerodynamic model

%save these analysis modules for calling later...
save('run_ONERA_a320.mat', 'run_ONERA')
save('run_ext_a320.mat', 'run_ext')

%%