%% Title section - Plots for comparing numerical predictions with the experimental data for tests G1.1, G1.2, and G2
%{
Comments:
* The experimental data 'exprData' is read from a .mat file created by Sanuja (wing root pitch angles included).
* Numerical predictions can be performed using the functions G1_1(), G1_2(), and G2().
* The inputs for the numerical predictions are the deterministic pitch angle and the uncertain parameters (EI, GJ, Sxx, Szz).
* Finally, the numerical predictions and the experimental data are plotted together (for comparison).
%}

clear all; close all;

%add folder below to path - contains functions relating to the experimental
%data and test-specific computations
addpath('groundTests'); 
load('run_ONERA.mat'); run = run_ONERA;

%set to work with modal form...
run = run.setTransform('modal', 14); 

%% load the data folder..

%check the structrued data set in the .mat file -should be self
%explanatory, test cases are referenced as per the document
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');

%% ground test case..G1.1

angl = linspace(0, 90, 15)*pi/180;

%function in the groundTests folder, labelled based on the test, take in
%the 'run' object and returns the test-specific measurements
%{
The uncertain parameters can be varied using multiplicative scaling factors, as follows:
* 'EI': scaling factor for out-of-plane EI ('EI'=1 corresponds to the tuned EI_{yy} value in Sanuja's paper)
* 'GJ': scaling factor for GJ ('GJ'=1 corresponds to the tuned GJ value in Sanuja's paper)
%}
[delta_LE, delta_TE, beta_y, beta_x] = G1_1(run, angl, 'EI', 1, 'GJ', 1);

tip_dispFig = figure('Name', 'tip deflections'); 
plot(angl*180/pi, delta_LE, 'r-'); hold on;
plot(angl*180/pi, delta_TE, 'r--'); hold on;
plot(exprData{1}.rootAngl, exprData{1}.delta_LE, 'rs', 'MarkerFaceColor',...
    'r'); hold on;
plot(exprData{1}.rootAngl, exprData{1}.delta_TE, 'rs'); hold on;
xlabel('Root angle, [deg]'); ylabel('\delta_{\bullet}, [m]'); hold on;
grid minor;

rootLoadFig = figure('Name', 'scaled load measurments');
subplot(1,2,1)
plot(angl*180/pi, beta_y, 'r-'); hold on;
plot(exprData{1}.rootAngl, exprData{1}.beta_y, 'rs'); hold on;
xlabel('Root angle, [deg]'); ylabel('\beta_y, [1/m]'); hold on;
grid minor;

subplot(1,2,2)
plot(angl*180/pi, beta_x, 'r-'); hold on;
plot(exprData{1}.rootAngl, exprData{1}.beta_x, 'rs'); hold on;
xlabel('Root angle, [deg]'); ylabel('\beta_x, [1/m]'); hold on;
grid minor;

drawnow;

%% G1.2..
%{
The uncertain parameters can be varied using multiplicative scaling factors, as follows:
* 'EI': scaling factor for out-of-plane EI ('EI'=1 corresponds to the tuned EI_{yy} value in Sanuja's paper)
* 'GJ': scaling factor for GJ ('GJ'=1 corresponds to the tuned GJ value in Sanuja's paper)
%}
[delta_LE, delta_TE, beta_y, beta_x] = G1_2(run, angl, 'EI', 1, 'GJ', 1);

set(0, 'CurrentFigure', tip_dispFig)
plot(angl*180/pi, delta_LE, 'b-'); hold on;
plot(angl*180/pi, delta_TE, 'b--'); hold on;
plot(exprData{2}.rootAngl, exprData{2}.delta_LE, 'bs', 'MarkerFaceColor',...
    'b'); hold on;
plot(exprData{2}.rootAngl, exprData{2}.delta_TE, 'bs'); hold on;
leg_stat = {'G1.1', [''], [''], [''], 'G1.2', [''], [''], ['']};
legend(leg_stat)

set(0, 'CurrentFigure', rootLoadFig)
subplot(1,2,1)
plot(angl*180/pi, beta_y, 'b-'); hold on;
plot(exprData{2}.rootAngl, exprData{2}.beta_y, 'bs'); hold on;
leg_statload = {'G1.1', [''], 'G1.2', ['']};
legend(leg_statload)

subplot(1,2,2)
plot(angl*180/pi, beta_x, 'b-'); hold on;
plot(exprData{2}.rootAngl, exprData{2}.beta_x, 'bs'); hold on;
legend(leg_statload)

drawnow;


%% G2
%{
The uncertain parameters can be varied using multiplicative scaling factors, as follows:
* 'EI': scaling factor for out-of-plane EI ('EI'=1 corresponds to the tuned EI_{yy} value in Sanuja's paper)
* 'GJ': scaling factor for GJ ('GJ'=1 corresponds to the tuned GJ value in Sanuja's paper)
* 'Sxx': scaling factor for m_{Sxx} ('Sxx'=1 corresponds to the tuned m_{Sxx} value in Sanuja's paper)
* 'Szz': scaling factor for m_{Szz} ('Szz'=1 corresponds to the tuned m_{Szz} value in Sanuja's paper)
The prediction corresponds to modal frequencies for the following modes of vibration
(see (Jayatilake), Supplementary ground test results of the highly flexible wing demonstrator (Appendix E)):
frqs[1]: first out-of-plane bending mode
frqs[2]: first in-plane bending mode
frqs[3]: second out-of-plane bending mode
frqs[4]: first torsion mode
%}
[frqs] = G2(run, angl, 'EI', 1, 'GJ', 1, 'Sxx', 1, 'Szz', 1);

figure;
for mode=1:4
    plot(angl*180/pi, frqs(mode,:), 'k-'); hold on;
    plot(exprData{3}.rootAngl, exprData{3}.frequencies(mode,:), 'bx'); hold on;
end
grid minor;
xlabel('Root Angle, [deg]');
ylabel('\omega, [Hz]')