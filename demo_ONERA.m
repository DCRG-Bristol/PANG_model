function demo_ONERA

clear all;

load('run_ONERA.mat');
run = run_ONERA;

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

%the below is re-implemented to account for the new sizing..
q0 = run.q0_struct; %this releases a zero-state vector with the correct size - for a structure-only problem (i.e no aerodynamics..)
qAero0 = run.q0_aeroStruct; %this releases a zero-state vector with the correct size - for a probelm with aero...

%if getStructModes is called now, theeigenvalectors will basically be an identity matrix!!

%% static solution... this can be identically handled with an analysis module derived with the extAeroBase

q = q0; %structure only problem so assign teh correctly sized initial vector
alp = linspace(pi/2,0,15);
Er = [60, 70, 80]*1e9;

clr = {'b', 'r', 'k'};

figure;
xMeasureStat = run.geom.L; %station to get displacements for plotting: the planform geometery class is inherited in the run modules for convinience.. here this is used to get the length

for k_idx = 1:length(Er)
    E = Er(k_idx);
    par_label{k_idx} = ['E = ',num2str(E/(1e9)), 'GPa'];

    for a_idx = 1:length(alp)
        pitch = alp(a_idx);
        q = fsolve(@(qstr)(run.structDeriv(qstr,'alpha0', pitch, 'E', E, 'alpha', pitch)), q); 

        [x, y, z] = run.getDisplField(q, xMeasureStat ,'beamModel'); %get displacements at x
        tipDefl(a_idx) = z(1,:); %z of size [2, length(xMeasureStat)] hold the z displacements in z(1,:) and TE in z(2,:)

        %eigenvalues....
        [~, e] = run.getStructModes(q, 'E', E, 'alpha0', pitch); %this calculates the structural modes with the reference parameters...
        EigVals(:,a_idx) = e;
    end
    subplot(1,2,1); plot(alp*180/pi, tipDefl, 'x-', 'color', clr{k_idx}); hold on; 
    subplot(1,2,2); 
    for mode=1:4
        plot(alp*180/pi, abs(EigVals(mode,:))./(2*pi), 'x-', 'color', clr{k_idx}); hold on;
    end
    drawnow;
end

subplot(1,2,1); xlabel('Pitch angle, [deg]'); ylabel('Tip Deflection, [m]');
legend(par_label)
subplot(1,2,2); xlabel('Pitch angle, [deg]'); ylabel('Modal frequencies, [Hz]')

%%
clc;
Ur = linspace(75,200,15);
q = qAero0; %initieslised zero vector for the aeroelastic problem..

pitch=5*pi/180;
for u_idx=1:length(Ur)
    U = Ur(u_idx);
    q = fsolve(@(q_all)(run.aero_structDeriv(q_all, 'alpha0', pitch,'alpha', pitch + 2.5*pi/180, 'U', U)), q); %note E is not set here... the same basline value assigned in line 8 is used

    [x, y, z] = run.getDisplField(q, xMeasureStat ,'beamModel'); %get displacements at x
    tipDefl(u_idx) = z(1,:); %z of size [2, length(xMeasureStat)] hold the z displacements in z(1,:) and TE in z(2,:)
end

figure; plot(Ur,tipDefl, 'x-'); xlabel('Airspeed, [m/s]'); ylabel('tip LE vertical deflection, [m]');

% plot the last static solution in aircraft and beam frames...
figure; 
xStat = linspace(0, run.geom.L, 50);

[x, y, z] = run.getDisplField(q, xStat ,'beamModel'); %get displacements at x
X = [x(1,:), flip(x(2,:))];
Y = [y(1,:), flip(y(2,:))];
Z = [z(1,:), flip(z(2,:))];
subplot(1,2,1)
plot3(X, Y, Z, 'x-', 'linewidth', 2); title('beam frame')
xlabel('x, [m]'); ylabel('y, [m]'); zlabel('z, [m]');
set(gca, 'dataaspectratio', [1,1,1])

[x, y, z] = run.getDisplField(q, xStat ,'aircraft'); %get displacements at x
X = [x(1,:), flip(x(2,:))];
Y = [y(1,:), flip(y(2,:))];
Z = [z(1,:), flip(z(2,:))];
subplot(1,2,2)
plot3(X, Y, Z, 'x-', 'linewidth', 2); title('aircraft frame')
set(gca, 'dataaspectratio', [1,1,1])
set(gca, 'zDir', 'reverse');
xlabel('x, [m]'); ylabel('y, [m]'); zlabel('z, [m]');


end