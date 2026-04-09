clear all; close all;

%add folder below to path - contains functions relating to the experimental
%data and test-specific computations
addpath('WTTests'); 
load('run_ONERA.mat'); run = run_ONERA;

%set to work with modal form...
run = run.setTransform('modal', 14); 


%% run for a selected angle

ang = 1.6*pi/180; %angle to run..
ang_expr = 1.9; %note angles to access experimental data in degrees

%function to retrieve processes experimental data...
[ang_true, exp_statResp, exp_Uf, exp_beta_yf, exp_beta_xf] = expr_statStab(ang_expr);
%ang_true - actual angle (closest to ang_expr) where experimnetal data is
%available and is accessed) - all other outputs are similar to model
%values...AVAILABLE ANGLES: [0, 0.3, 0.6, 1.1, 1.4, 1.9, 2.5] degrees

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
%plot numerical...
plot(statResp.U, statResp.beta_y, 'k-'); hold on;
plot(Uf, beta_yf, 'ro', 'markerfacecolor', 'r'); hold on;
%plot experimental...
plot(exp_statResp.U, exp_statResp.beta_y, 'bx-'); hold on;
plot(exp_Uf, exp_beta_yf, 'r^'); hold on;
xlabel('U, [m/s]'); ylabel('\beta_y, [1/m]');
legend({'Static response (model)', 'Hopf Bifurcations (model)', 'Static response (expr)', 'Hopf Bifurcations (expr.)'})

exp_mkrs = {'mx', 'ko', 'b^', 'rs'};
clrs = {'m', 'k', 'b', 'r'};
figure;
for mode=1:4
    subplot(2,1,1);
    plot(statResp.U, statResp.damp(mode,:), 'color', clrs{mode}); hold on;
    plot(exp_statResp.U, exp_statResp.damp(mode,:), exp_mkrs{mode},...
        'markerFaceColor', clrs{mode}); hold on;

    subplot(2,1,2);
    plot(statResp.U, statResp.frqs(mode,:), clrs{mode}); hold on;
    plot(exp_statResp.U, exp_statResp.frqs(mode,:), exp_mkrs{mode},...
        'markerFaceColor', clrs{mode}); hold on;
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
