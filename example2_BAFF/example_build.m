%% example for building the model...

%clear all;

%this script require a 'BAFF' object as an input - assigned to the Baff_mdl variable. See https://github.com/dcrg-bristol/baff

% ADP = load('UB321_baseline_simple.mat'); 
% ADP = ADP(1).ADP.BuildBaff;

Baff_mdl = ADP.Baff; % see 'Basic Model description' in https://github.com/DCRG-Bristol/Sprint0626

%%
buildOb = buildSystem.buildBase; %base class for model building...
buildOb.fileLocation = fileparts(mfilename('fullpath')); %set the file location where a folder with a model-specific items could be created...can be changed...
buildOb.name = 'a320_likeWing'; %a name for the model

%% basis properties

if ADP.HingeEta==1
    buildOb = buildOb.baff2PANG(Baff_mdl, 'Wing_RHS'); %this builds a PANG model definition from a baff. 'Wing_RHS' is the name of the item in the baff input to use for generating the model in this case..
else
    buildOb = buildOb.baff2PANG(Baff_mdl, 'Wing_RHS', 'FFWT_RHS'); %add wing tip if needed
end

%OPTIONAL - assign system size...
buildOb.basis.Nw = 4; %# of out of plane bending fcns..
buildOb.basis.Nv = 4; %.. chord wise bending
buildOb.basis.Nthet = 4;%...torsion

%set the arc-length wise stations (easured along the beam's elastic axis)
%to be used for printing the aerodynamic mesh (calculated at the
%corresponding LE/TE pairs for each entry in xi). The loads are applied at
%the mid-points along this mesh, at the quarter chord position.
buildOb.basis.xi = linspace(0, buildOb.geom.L, 20);

%% create the model...

buildOb.prepFolder; %this method sets up the folder for saving model-specific contents.. 
buildOb.writeModel; %thus generates the model 

%% derive analysis models... the buildBass instance is passed in 

run_ONERA = analysis.oneraBase(buildOb); %analysis module for the embedded-aerodynamic model
run_ext = analysis.extAeroBase(buildOb); %analysis module for the external aerodynamic model (e.g. external UVLM module)

%save these analysis modules for calling later.... once tehse are made, 
save('run_ONERA_a320.mat', 'run_ONERA')
save('run_ext_a320.mat', 'run_ext')

%%