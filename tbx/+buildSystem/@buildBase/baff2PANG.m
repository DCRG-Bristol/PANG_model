%script to build PANG model with baff file

function buildObj = baff2PANG(buildObj, baffObj, varargin)

%inputs:
%baffObj: baff model for aircraft..(baff class)
%itmName: name of wing object in baff (string)

EA0 = [];
eta_beam = [];
globLgth =0;

%find the matchig item in baff model...
for i=1:length(varargin)
    id(i) = find(strcmp([baffObj.Wing.Name],varargin{i}));
    baffItm(i) = baffObj.Wing(id(i)); %baff item to extract properties from....

    %collect all ETAS and  combine
    EA0 = [EA0, [baffItm(i).GetGlobalPos([baffItm(i).Stations.Eta], [0;0;0])]]; %elastic axis in BAFF coordinates
    eta_beam = [eta_beam, (i-1)+[baffItm(i).Stations.Eta]];
    globLgth = globLgth + baffItm(i).EtaLength;
end

%remove duplicates...
[eta_beam, idx] = unique(eta_beam);
EA0 = EA0(:,idx);
beam_idx = idx; %unique beam indices to keep..used later

%get elastic axis coordinates.. aircraft frame...
%EA0 = baffItm.GetGlobalPos([baffItm.Stations.Eta], [0;0;0]); %elastic axis in BAFF coordinates
r0 = EA0(:,1); %record of wing root in weird baff frame
EA0 = EA0-EA0(:,1); %offset wing root
EA_r = [abs(EA0(2,:)); -EA0(1,:)]; %PANG frame....

sweep = asin(-EA_r(2,end)/EA_r(1,end)); %estimate beam sweep

% beam loci...beam arc length coordinates...
s_statn = EA_r(:,end)'*EA_r/norm(EA_r(:,end)); %beam stations.. project to beam: arc length

%aero dynamic/planfome__________________________________________________________

%get aero aerometry at beam stations...
xi = [];
cLE =[];
cTE = [];
LE_x = [];
TE_x = [];
twistDat = [];
for i=1:length(varargin)
    xi = [xi, ...
        interp1(eta_beam, s_statn, i-1+[baffItm(i).AeroStations.Eta])]; %aero mesh...

    %LE/TE offsets...
    cLE = [cLE,...
        ([baffItm(i).AeroStations.Chord].*[baffItm(i).AeroStations.BeamLoc])*cos(sweep)];
    cTE = [cTE,...
        ([baffItm(i).AeroStations.Chord].*(1-[baffItm(i).AeroStations.BeamLoc]))*cos(sweep)];
    twistDat = [twistDat, [baffItm(i).AeroStations.Twist]];

    %x-coordinates of LE and TE....
    LE_x = [LE_x, interp1(eta_beam, EA_r(1,:), i-1+[baffItm(i).AeroStations.Eta])]; 
    TE_x = [TE_x, interp1(eta_beam, EA_r(1,:), i-1+[baffItm(i).AeroStations.Eta])]; 
end
LE_x = LE_x  + cLE*sin(sweep);
TE_x = TE_x - cTE*sin(sweep);

%remoe duplicates...
[xi, idx] = unique(xi);
cLE = cLE(idx);
cTE = cTE(idx);
LE_x = LE_x(idx);
TE_x = TE_x(idx);
twistDat = twistDat(idx);

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
    basis.Nw = 4;
    basis.Nv = 4;
    basis.Nthet = 4;
    basis.xi = xi;
    buildObj.basis = basis;
else
    buildObj.basis.xi = xi;
end

%beam properties__________________________________________________________
counter = 1;
for i=1:length(varargin)
    for stn=1:length([baffItm(i).Stations])
        I_matr = baffItm(i).Stations(stn).I;
        J = baffItm(i).Stations(stn).J;
        rho = baffItm(i).Stations(stn).Mat.rho;
        E = baffItm(i).Stations(stn).Mat.E;
        G = baffItm(i).Stations(stn).Mat.G;

        EI1(counter) = E*I_matr(2,2);
        EI2(counter) = E*I_matr(3,3);
        EI12(counter) = E*I_matr(2,3);
        GJ(counter) = G*J;

        m(counter) = rho*baffItm(i).Stations(stn).A;
        mxx(counter) = rho*(I_matr(2,2) + I_matr(3,3));
        mzz(counter) = rho*(I_matr(3,3));
        myy(counter) = rho*(I_matr(2,2));
        counter = counter + 1;
    end
end

%only keep the unique ones...
idx = beam_idx;
EI1 = EI1(idx);
EI2 = EI2(idx);
EI12 = EI12(idx);
GJ = GJ(idx);
m = m(idx);
mxx = mxx(idx);
mzz = mzz(idx);
myy=  myy(idx);

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

for j=1:length(varargin)
    for i=1:length(baffItm(j).Children)
        type = class(baffItm(j).Children(i));

        switch type
            case 'baff.Mass'
                mass_idx = mass_idx+1;
                mass_elem(mass_idx) = buildSystem.structure.descrElem;
                mass_elem(mass_idx).m = baffItm(j).Children(i).mass;
                mass_elem(mass_idx).xp =...
                    interp1(eta_beam, s_statn, j-1+baffItm(j).Children(i).Eta);

                %find offsets from nearest beam...
                rawPostn = baffItm(j).GetGlobalPos(baffItm(j).Children(i).Eta, [0;0;0])-r0;
                for dim=1:3
                    beamPostn(dim) = interp1(eta_beam(:)', EA0(dim,:), j-1+baffItm(j).Children(i).Eta);
                end
                yOfst = -(rawPostn(1) - beamPostn(1))*cos(sweep);
                mass_elem(mass_idx).e = -yOfst/b(mass_elem(mass_idx).xp)+a(mass_elem(mass_idx).xp);
                mass_elem(mass_idx).mxx = 0;
                mass_elem(mass_idx).myy = 0;
                mass_elem(mass_idx).mzz = 0;

            case 'cast.drag.DraggableBluffBody'
                for mm=1:length(baffItm(j).Children(i).Children)
                    mass_idx = mass_idx+1;
                    mass_elem(mass_idx) = buildSystem.structure.descrElem;
                    mass_elem(mass_idx).m = (baffItm(j).Children(i).Children(mm).mass);
                    mass_elem(mass_idx).xp =...
                        interp1(eta_beam, s_statn, j-1+baffItm(j).Children(i).Eta);

                    %find offsets from nearest beam...
                    rawPostn = baffItm(j).Children(i).Children(mm).GetGlobalPos(baffItm(j).Children(i).Children(mm).Eta, [0;0;0])-r0;
                    for dim=1:3
                        beamPostn(dim) = interp1(eta_beam(:)', EA0(dim,:), j-1+baffItm(j).Children(i).Children(mm).Eta);
                    end
                    yOfst = -(rawPostn(1) - beamPostn(1))*cos(sweep);

                    mass_elem(mass_idx).e = -yOfst/b(mass_elem(mass_idx).xp)+a(mass_elem(mass_idx).xp);
                    mass_elem(mass_idx).mxx = 0;
                    mass_elem(mass_idx).myy = 0;
                    mass_elem(mass_idx).mzz = 0;
                end

            case 'baff.Fuel'
                fuel_idx = fuel_idx+1;
                fuel_elem(fuel_idx) = buildSystem.structure.descrElem;
                fuel_elem(fuel_idx).m = sum(baffItm(j).Children(i).mass);
                fuel_elem(fuel_idx).xp =...
                    interp1(eta_beam, s_statn, j-1+baffItm(j).Children(i).Eta);

                %find offsets from nearest beam...
                rawPostn = baffItm(j).GetGlobalPos(baffItm(j).Children(i).Eta, [0;0;0])-r0;
                for dim=1:3
                    beamPostn(dim) = interp1(eta_beam(:)', EA0(dim,:), j-1+baffItm(j).Children(i).Eta);
                end
                yOfst = -(rawPostn(1) - beamPostn(1))*cos(sweep);

                fuel_elem(fuel_idx).e = -yOfst/b(fuel_elem(fuel_idx).xp)+a(fuel_elem(fuel_idx).xp);
                fuel_elem(fuel_idx).mxx = 0;
                fuel_elem(fuel_idx).myy = 0;
                fuel_elem(fuel_idx).mzz = 0;
        end
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

