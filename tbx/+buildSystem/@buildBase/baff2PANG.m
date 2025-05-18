%script to build PANG model with baff file

function buildObj = baff2PANG(buildObj, baffObj, itmName)

%inputs:
%baffObj: baff model for aircraft..(baff class)
%itmName: name of wing object in baff (string)

%find the matchig item in baff model...
id = find(strcmp([baffObj.Wing.Name],itmName));
baffItm = baffObj.Wing(id); %baff item to extract properties from....

%get elastic axis coordinates.. aircraft frame...
EA0 = baffItm.GetGlobalPos([baffItm.Stations.Eta], [0;0;0]); %elastic axis in BAFF coordinates
EA0 = EA0-EA0(:,1); %offset wing root
EA_r = [abs(EA0(2,:)); -EA0(1,:)]; %PANG frame....

sweep = asin(-EA_r(2,end)/EA_r(1,end)); %estimate beam sweep

% beam loci...beam arc length coordinates...
s_statn = EA_r(:,end)'*EA_r/norm(EA_r(:,end)); %beam stations.. project to beam: arc length

%aero dynamic/planfome__________________________________________________________

%get aero aerometry at beam stations...
xi = interp1([baffItm.Stations.Eta], s_statn, [baffItm.AeroStations.Eta]); %aero mesh...

%LE/TE offsets...
cLE = [[baffItm.AeroStations.Chord].*[baffItm.AeroStations.BeamLoc]]*cos(sweep);
cTE = [[baffItm.AeroStations.Chord].*(1-[baffItm.AeroStations.BeamLoc])]*cos(sweep);
twistDat = [baffItm.AeroStations.Twist];

%x-coordinates of LE and TE....
LE_x = interp1([baffItm.Stations.Eta], EA_r(1,:), [baffItm.AeroStations.Eta]) + cLE*sin(sweep);
TE_x = interp1([baffItm.Stations.Eta], EA_r(1,:), [baffItm.AeroStations.Eta]) - cTE*sin(sweep);

%find place where TE_x=0 to clip...
x_TEClip = interp1(TE_x, xi, 0, 'linear','extrap');
sel_idx = find((xi>x_TEClip));

if x_TEClip(1)>0
    cTE = [0, interp1(xi, cTE, x_TEClip), cTE(sel_idx)];
    cLE = [interp1(xi, cLE, 0, 'linear','extrap'), interp1(xi, cLE, x_TEClip, 'linear','extrap'), cLE(sel_idx)];
    twistDat = [interp1(xi, twistDat, 0, 'linear','extrap'), interp1(xi, twistDat, x_TEClip, 'linear','extrap'), twistDat(sel_idx)];
    xi = [0, x_TEClip, xi(sel_idx)];
    [xi, odr] = sort(xi);
    cLE = cLE(odr);
    cTE = cTE(odr);
    twistDat = twistDat(odr);
else
    cTE = [interp1(xi, cTE, x_TEClip), cTE(sel_idx)];
    cLE = [interp1(xi, cLE, x_TEClip, 'linear','extrap'), cLE(sel_idx)];
    twistDat = [interp1(xi, twistDat, x_TEClip, 'linear','extrap'), twistDat(sel_idx)];
    xi = [x_TEClip, xi(sel_idx)];
    [xi, odr] = sort(xi);
    cLE = cLE(odr);
    cTE = cTE(odr);
    twistDat = twistDat(odr);
end

b_dat = 0.5*(cTE+cLE);

cLEf = @(x)(interp1(xi, cLE, x));

b = @(x)(interp1(xi, b_dat, x));
a = @(x)(cLEf(x)./b(x) - 1);

%set class property: geometry>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
geom = buildSystem.geom;
geom.a = a; %a; %position of the elastic axis... theodorsen definition
geom.b = b; %semi-chord distribution
geom.sweep = sweep; %sweep (of beam line!)
geom.L = s_statn(end);
geom.twist = @(x)(interp1(xi, twistDat, x));
buildObj.geom = geom;

%set class property: basis for aero >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if isempty(buildObj.basis)
    basis = buildSystem.basis;
    basis.xi = xi;
    buildObj.basis = basis;
else
    buildObj.basis.xi = xi;
end

%beam properties__________________________________________________________
for stn=1:length(s_statn)
    I_matr = baffItm.Stations(stn).I;
    J = baffItm.Stations(stn).J;
    rho = baffItm.Stations(stn).Mat.rho;
    E = baffItm.Stations(stn).Mat.E;
    G = baffItm.Stations(stn).Mat.G;

    EI1(stn) = E*I_matr(2,2);
    EI2(stn) = E*I_matr(1,1);
    EI12(stn) = E*I_matr(1,2);
    GJ(stn) = G*J;

    m(stn) = rho*baffItm.Stations(stn).A;
    mxx(stn) = rho*(I_matr(1,1) + I_matr(2,2));
    mzz(stn) = rho*(I_matr(1,1));
    myy(stn) = rho*(I_matr(2,2));
end

%set class property: elastic >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
elas(1) = buildSystem.structure.elasBase; %call elastic property class
elas.GJ = @(x)interp1(s_statn, GJ, x);
elas.EI1 = @(x)interp1(s_statn, EI1, x);
elas.EI2 = @(x)interp1(s_statn, EI2, x);
elas.EI12 = @(x)interp1(s_statn, EI12, x);
elas(1).name = 'baff_struct'; %namee for this matrix..
elas(1).fctrId = [];

%set class property: inertia: fixed >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
inertia(1) = buildSystem.structure.inertiaBase; %inertia property class
inertia(1).e = a; %position of mass axis.. theodorsen
inertia(1).m = @(x)interp1(s_statn, m, x); %mass per unit length [kg/m]
inertia(1).mxx = @(x)interp1(s_statn, mxx, x); %twisting inertia about chordwise about the CoM [kgm]
inertia(1).mzz = @(x)interp1(s_statn, mzz, x); %chord-plane rotation inertia [kgm]
inertia(1).myy = @(x)interp1(s_statn, myy, x); %rotation inertia..bending rotation [kgm]
inertia(1).name = 'baff_struct'; %namee for this matrix..
inertia(1).fctrId = [];

%go through children and group types..
fuel_idx = 0;
mass_idx = 0;

for i=1:length(baffItm.Children)
    type = class(baffItm.Children(i));

    switch type
        case 'baff.Mass'
            mass_idx = mass_idx+1;
            mass_elem(mass_idx) = buildSystem.structure.descrElem;
            mass_elem(mass_idx).m = baffItm.Children(i).mass;
            mass_elem(mass_idx).xp =...
                interp1([baffItm.Stations.Eta], s_statn, baffItm.Children(i).Eta);
            mass_elem(mass_idx).e =...
                -baffItm.Children(i).Offset(1)/b(mass_elem(mass_idx).xp)-1;
            mass_elem(mass_idx).mxx = 0;
            mass_elem(mass_idx).myy = 0;
            mass_elem(mass_idx).mzz = 0;

        case 'baff.DraggableBluffBody'
            mass_idx = mass_idx+1;
            mass_elem(mass_idx) = buildSystem.structure.descrElem;
            mass_elem(mass_idx).m = sum(baffItm.Children(i).Children.mass);
            mass_elem(mass_idx).xp =...
                interp1([baffItm.Stations.Eta], s_statn, baffItm.Children(i).Eta);
            mass_elem(mass_idx).e =...
                -baffItm.Children(i).Offset(1)/b(mass_elem(mass_idx).xp)-1;
            mass_elem(mass_idx).mxx = 0;
            mass_elem(mass_idx).myy = 0;
            mass_elem(mass_idx).mzz = 0;

        case 'baff.Fuel'
            fuel_idx = fuel_idx+1;
            fuel_elem(fuel_idx) = buildSystem.structure.descrElem;
            fuel_elem(fuel_idx).m = sum(baffItm.Children(i).mass);
            fuel_elem(fuel_idx).xp =...
                interp1([baffItm.Stations.Eta], s_statn, baffItm.Children(i).Eta);
            fuel_elem(fuel_idx).e =...
                -baffItm.Children(i).Offset(1)/b(fuel_elem(fuel_idx).xp)-1;
            fuel_elem(fuel_idx).mxx = 0;
            fuel_elem(fuel_idx).myy = 0;
            fuel_elem(fuel_idx).mzz = 0; 
    end
end

%masses compute equivalent distributed
if mass_idx>0
    inertia(1).elem = mass_elem;
end

if fuel_idx>0
    inertia(2) = buildSystem.structure.inertiaBase; %inertia property class
    inertia(2).e = @(x)0; %position of mass axis.. theodorsen
    inertia(2).m = @(x)0; %mass per unit length [kg/m]
    inertia(2).mxx = @(x)0; %twisting inertia about chordwise about the CoM [kgm]
    inertia(2).mzz = @(x)0; %chord-plane rotation inertia [kgm]
    inertia(2).myy = @(x)0; %rotation inertia..bending rotation [kgm]
    inertia(2).name = 'fuel'; %namee for this matrix..
    inertia(2).fctrId = 'fuel_frac';
    inertia(2).elem = fuel_elem;

    fprintf('fuel masses found: ...embedding user parameter fuel_frac\n')
end

%set class property: embed to buildBasis class >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
buildObj.inertia = inertia;
buildObj.elas = elas;
grav = buildObj.inertia2grav(inertia);
buildObj.grav = grav;

end

