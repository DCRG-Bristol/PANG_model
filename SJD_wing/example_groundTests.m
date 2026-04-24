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

angl = linspace(0, 90, 30)*pi/180;

%function in the groundTests folder, labelled based on the test, take in
%the 'run' object and returns the test-specific measurements
%{
The uncertain parameters can be varied using multiplicative scaling factors, as follows:
* 'EI': scaling factor for out-of-plane EI ('EI'=1 corresponds to the tuned EI_{yy} value in Sanuja's paper)
* 'GJ': scaling factor for GJ ('GJ'=1 corresponds to the tuned GJ value in Sanuja's paper)
%}
[delta_LE, delta_TE, beta_y, beta_x] = G1_1(run, angl, 'EI', 1.0143, 'GJ', 1.0789);

tip_dispFig = figure('Name', 'tip deflections'); 
plot(angl*180/pi, delta_LE-delta_LE(1), 'r-'); hold on;
plot(angl*180/pi, delta_TE-delta_TE(1), 'r--'); hold on;
errorbar(exprData{1}.rootAngl, exprData{1}.delta_LE-exprData{1}.delta_LE(1), 2e-3, 'rs', 'MarkerFaceColor',...
    'r'); hold on;
errorbar(exprData{1}.rootAngl, exprData{1}.delta_TE-exprData{1}.delta_TE(1), 2e-3, 'rs'); hold on;
xlabel('Root angle, [deg]'); ylabel('\delta_{\bullet}, [m]'); hold on;
grid minor;

load('ground_strains_lower_and_upper_bounds.mat');
rootLoadFig = figure('Name', 'scaled load measurments');
subplot(1,2,1)
plot(angl*180/pi, beta_y, 'r-'); hold on;
errorbar(exprData{1}.rootAngl, exprData{1}.beta_y, exprData{1}.beta_y-[0, beta_y_g1_1_lb], [0, beta_y_g1_1_ub]-exprData{1}.beta_y, 'rs'); hold on;
xlabel('Root angle, [deg]'); ylabel('\beta_y, [1/m]'); hold on;
grid minor;

subplot(1,2,2)
plot(angl*180/pi, beta_x, 'r-'); hold on;
errorbar(exprData{1}.rootAngl, exprData{1}.beta_x, exprData{1}.beta_x-[0, beta_x_g1_1_lb], [0, beta_x_g1_1_ub]-exprData{1}.beta_x, 'rs'); hold on;
xlabel('Root angle, [deg]'); ylabel('\beta_x, [1/m]'); hold on;
grid minor;

drawnow;

%% G1.2
%{
The uncertain parameters can be varied using multiplicative scaling factors, as follows:
* 'EI': scaling factor for out-of-plane EI ('EI'=1 corresponds to the tuned EI_{yy} value in Sanuja's paper)
* 'GJ': scaling factor for GJ ('GJ'=1 corresponds to the tuned GJ value in Sanuja's paper)
%}
[delta_LE_, delta_TE_, beta_y, beta_x] = G1_2(run, angl, 'EI', 1.0143, 'GJ', 1.0789);

set(0, 'CurrentFigure', tip_dispFig)
plot(angl*180/pi, delta_LE_-delta_LE(1), 'b-'); hold on;
plot(angl*180/pi, delta_TE_-delta_TE(1), 'b--'); hold on;
errorbar(exprData{2}.rootAngl, exprData{2}.delta_LE-exprData{1}.delta_LE(1), 2e-3, 'bs', 'MarkerFaceColor',...
    'b'); hold on;
errorbar(exprData{2}.rootAngl, exprData{2}.delta_TE-exprData{1}.delta_TE(1), 2e-3, 'bs'); hold on;
leg_stat = {'G1.1', '', '', '', '', '','', '', '','', '', '', 'G1.2'};
legend(leg_stat)

set(0, 'CurrentFigure', rootLoadFig)
subplot(1,2,1)
plot(angl*180/pi, beta_y, 'b-'); hold on;
errorbar(exprData{2}.rootAngl, exprData{2}.beta_y, exprData{2}.beta_y-[0, beta_y_g1_2_lb], [0, beta_y_g1_2_ub]-exprData{2}.beta_y, 'bs'); hold on;
leg_statload = {'G1.1', [''], 'G1.2', ['']};
legend(leg_statload)

subplot(1,2,2)
plot(angl*180/pi, beta_x, 'b-'); hold on;
errorbar(exprData{2}.rootAngl, exprData{2}.beta_x, exprData{2}.beta_x-[0, beta_x_g1_2_lb], [0, beta_x_g1_2_ub]-exprData{2}.beta_x, 'bs'); hold on;
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
[frqs] = G2(run, angl, 'EI', 1.0143, 'GJ', 1.0789, 'Sxx', 1.1157, 'Szz', 1.1157);
load('GVT_freq_lower_and_upper_bounds.mat');

figure;
plot(angl*180/pi, frqs(1,:), 'k-'); hold on;
errorbar(exprData{3}.rootAngl, exprData{3}.frequencies(1,:), exprData{3}.frequencies(1,:)-oop1_mode_freq_lb, oop1_mode_freq_ub-exprData{3}.frequencies(1,:), 'bx'); hold on;
plot(angl*180/pi, frqs(2,:), 'k-'); hold on;
errorbar(exprData{3}.rootAngl, exprData{3}.frequencies(2,:), exprData{3}.frequencies(2,:)-ip1_mode_freq_lb, ip1_mode_freq_ub-exprData{3}.frequencies(2,:), 'bx'); hold on;
plot(angl*180/pi, frqs(3,:), 'k-'); hold on;
errorbar(exprData{3}.rootAngl, exprData{3}.frequencies(3,:), exprData{3}.frequencies(3,:)-oop2_mode_freq_lb, oop2_mode_freq_ub-exprData{3}.frequencies(3,:), 'bx'); hold on;
plot(angl*180/pi, frqs(4,:), 'k-'); hold on;
errorbar(exprData{3}.rootAngl, exprData{3}.frequencies(4,:), exprData{3}.frequencies(4,:)-tor1_mode_freq_lb, tor1_mode_freq_ub-exprData{3}.frequencies(4,:), 'bx'); hold on;
grid minor;
xlabel('Root Angle, [deg]');
ylabel('\omega, [Hz]')