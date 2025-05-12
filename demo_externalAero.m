
%% example
clear all;
close all;

buildOb = buildSystem.buildBase;
%buildOb.fileLocation = cd;
buildOb.name = 'a320_likeWing';
buildOb.par = {'E'}; %DECLARE ALL USER-EMBEDDED VARIABLES...

%% 
L = 17.05;
cr = 5.82;
ct = 1.40;

%% planform geometry...

geom = buildSystem.geom;
geom.L = L;

%%
geom.a = @(x)(-0.05);
geom.b = @(x)(0.5*(wing.cr - (cr-ct)*(x/L)));
geom.sweep = 25*pi/180;
geom.twist = @(x)((2*pi/180) - ((2*pi/180 + 2*pi/180)*(x/L)));

buildOb.geom = geom;

%% basis properties

basis = buildSystem.basis;
basis.Nw = 4;
basis.Nv = 4;
basis.Nthet = 4;
basis.xi = linspace(0, geom.L, 21);

buildOb.basis = basis;

%% elastic properties

elas(1) = buildSystem.structure.elasBase;
elas(1).EI1 = @(x)(0.12   - (0.12-0.01)*(x/L));
elas(1).EI2 = @(x)(0.08   - (0.08-0.005)*(x/L));
elas(1).EI12 = @(x)(0);
elas(1).GJ = @(x)(0);
elas(1).name = 'K_E';
elas(1).fctrId = 'E';

elas(2) = buildSystem.structure.elasBase;
elas(2).EI1 = @(x)(0);
elas(2).EI2 = @(x)(0);
elas(2).EI12 = @(x)(0);
elas(2).GJ = @(x)(26e9*(0.05-(0.05-0.002)*(x/L)));
elas(2).name = 'K_G';
elas(2).fctrId = [];

buildOb.elas = elas;

%% add inertia

inertia = buildSystem.structure.inertiaBase;
inertia.e = @(x)(0.1); %CoG nondim position [-]
inertia.m = @(x)(300  - (300-100)*(x/L)); %mass per unit length [kg/m]
inertia.mxx = @(x)(20   - (20-2)*(x/L)); %twisting inertia, CoG [kgm]
inertia.mzz = @(x)(20   - (20-2)*(x/L)); %choord-plane rotation inertia [kgm]
inertia.myy = @(x)(0);

elem(1) = buildSystem.structure.descrElem;
elem(1).m = 2300;
elem(1).mxx = 0;
elem(1).mzz = 0;
elem(1).myy = 0;
elem(1).e = -0.07;
elem(1).xp = 0.3*L;
inertia(1).elem = elem;
inertia(1).name = 'massModel';
inertia(1).fctrId = [];

buildOb.inertia = inertia;

%% gravity..

grav = buildOb.inertia2grav(inertia); %make gravity property with inertia properties
buildOb.grav = grav;

%%
 
buildOb.prepFolder;
buildOb.writeModel;
run = analysis.extAeroBase(buildOb);
save('run.mat', 'run')



%%