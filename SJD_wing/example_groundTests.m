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
[delta_LE, delta_TE, beta_y, beta_x] = G1_1(run, angl, 'g', 9.81);

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

[delta_LE, delta_TE, beta_y, beta_x] = G1_2(run, angl, 'g', 9.81);

set(0, 'CurrentFigure', tip_dispFig)
plot(angl*180/pi, delta_LE, 'b-'); hold on;
plot(angl*180/pi, delta_TE, 'b--'); hold on;
plot(exprData{2}.rootAngl, exprData{2}.delta_LE, 'bs', 'MarkerFaceColor',...
    'b'); hold on;
plot(exprData{2}.rootAngl, exprData{2}.delta_TE, 'bs'); hold on;

set(0, 'CurrentFigure', rootLoadFig)
subplot(1,2,1)
plot(angl*180/pi, beta_y, 'b-'); hold on;
plot(exprData{2}.rootAngl, exprData{2}.beta_y, 'bs'); hold on;

subplot(1,2,2)
plot(angl*180/pi, beta_x, 'b-'); hold on;
plot(exprData{2}.rootAngl, exprData{2}.beta_x, 'bs'); hold on;

drawnow;

%% G2

[frqs] = G2(run, angl, 'g', 9.81, 'Sxx', 1);

figure;
for mode=1:4
    plot(angl*180/pi, frqs(mode,:), 'k-'); hold on;
    plot(exprData{3}.rootAngl, exprData{3}.frequencies(mode,:), 'bx'); hold on;
end
grid minor;
xlabel('Root Angle, [deg]');
ylabel('\omega, [Hz]')