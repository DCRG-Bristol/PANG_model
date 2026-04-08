clear all; close all;

%add folder below to path - contains functions relating to the experimental
%data and test-specific computations
addpath('WTTests'); 
load('run_ONERA.mat'); run = run_ONERA;

%set to work with modal form...
run = run.setTransform('modal', 14); 


%% run for a selected angle

ang = 1.0*pi/180; %angle to run..

%function to run static respnses and stability checks...
[statResp, Uf, beta_yf, beta_xf] = W_statStab(run, false, ang, 0.275, 0.44, 'EI', 1, 'GJ', 1, 'Sxx', 1, 'Szz', 1);
%inputs as below>>>>
%run - run object
%ang - wing root pitch angle
%true/false - to use continuation or not..
%'0.275' - lambda_L
%'0.44' - alpha_L
%followed by..'string'-value pairs for structural parameters..

% see plotting exmaple below for output detailss>>>>>>>>>>>>>>>>>>>>>>>>
%static repsonse...
figure;
plot(statResp.U, statResp.beta_y, 'k-'); hold on;
plot(Uf, beta_yf, 'ro', 'markerfacecolor', 'r');
xlabel('U, [m/s]'); ylabel('\beta_y, [1/m]');
legend({'Static response', 'Hopf Bifurcations'})

figure;
for mode=1:4
    subplot(2,1,1);
    plot(statResp.U, statResp.damp(mode,:)); hold on;

    subplot(2,1,2);
    plot(statResp.U, statResp.frqs(mode,:)); hold on;
end
subplot(2,1,1); xline(Uf); xlabel('U, [m/s]'); ylabel('\zeta, [-]');
set(gca, 'xAxisLocation', 'origin')
subplot(2,1,2); xline(Uf); xlabel('U, [m/s]'); ylabel('\omega, [Hz]');

%% work in progress - potential option with COCO...

% %run for a fixed angle case...
% [statResp, Uf, beta_yf, beta_xf, coco_log] = W_statStab(run, true, ang, 0.275, 0.44, 'Sxx', 1);
% 
% %compute entire boundary...
% [U_bd, beta_ybd, beta_xbd] = W_stabEnvel(run, coco_log);
