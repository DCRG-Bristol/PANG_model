clear all; close all;

%add folder below to path - contains functions relating to the experimental
%data and test-specific computations
addpath('WTTests'); 
load('run_ONERA.mat'); run = run_ONERA;

%set to work with modal form...
run = run.setTransform('modal', 14); 


%% run for a selected angle

ang = 1.2*pi/180; %angle to run..
ang_expr = 1.4; %note angles to access experimental data in degrees

%function to retrieve processes experimental data...
[ang_true, exp_statResp, exp_Uf, exp_beta_yf, exp_beta_xf] = expr_statStab(ang_expr);
%ang_true - actual angle (closest to ang_expr) where experimnetal data is
%available and is accessed) - all other outputs are similar to model
%values...AVAILABLE ANGLES: [0, 0.3, 0.6, 1.1, 1.4, 1.9, 2.5] degrees

load('particleswarm_mod_lse_optim_rand_strain_unc.mat', 'x_opt_lse')
%function to run static respnses and stability checks...
for ii=1:size(x_opt_lse, 1)
    [statResp(ii), Uf(ii, :), beta_yf(ii, :), beta_xf(ii, :)] = W_statStab(run, false, ang, 0.275, 0.44, 'EI', x_opt_lse(ii, 1), 'GJ', x_opt_lse(ii, 2), 'Sxx', x_opt_lse(ii, 3), 'Szz', x_opt_lse(ii, 3));
end
%inputs as below>>>>
%run - run object
%ang - wing root pitch angle
%true/false - to use continuation or not..
%'0.275' - lambda_L
%'0.44' - alpha_L
%followed by..'string'-value pairs for structural parameters..

%%
% see plotting exmaple below for output detailss>>>>>>>>>>>>>>>>>>>>>>>>
%static repsonse...
load('WT_strains_lower_and_upper_bounds.mat')
figure;
hold on

for ii = 1:size(x_opt_lse, 1)
    plot(statResp(1, ii).U, statResp(1, ii).beta_y, 'k-'); hold on;
    plot(Uf(ii, :), beta_yf(ii, :), 'ro', 'markerfacecolor', 'r'); hold on;
end    
plot(exp_statResp.U, exp_statResp.beta_y, 'bx-');
errorbar(exp_statResp.U, exp_statResp.beta_y, exp_statResp.beta_y-[0; beta_y_wt_experimental_data_set_angle_5_lb(2:end)], [0; beta_y_wt_experimental_data_set_angle_5_ub(2:end)]-exp_statResp.beta_y, 0.3, 0.3, 'bx-');
errorbar(exp_Uf, exp_beta_yf, [], [], 0.3, 0.3, 'r^');
xlabel('U, [m/s]'); ylabel('\beta_y, [1/m]');
legend({'Static response (model)', 'Hopf Bifurcations (model)', '', '','', '','', '','', '','', '','', '','', '','', '','', '','', '','', '','', '','', '','', '','', '','', '','', '','', '','', '','Static response (expr)', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'Hopf Bifurcations (expr.)'});    
hold off

%%

% see plotting exmaple below for output detailss>>>>>>>>>>>>>>>>>>>>>>>>
%static repsonse...
load('WT_strains_lower_and_upper_bounds.mat')
figure;
%plot numerical...
plot(statResp.U, statResp.beta_y, 'k-'); hold on;
plot(Uf, beta_yf, 'ro', 'markerfacecolor', 'r'); hold on;
%plot experimental...
if ang_expr == 1.4
    plot(exp_statResp.U, exp_statResp.beta_y, 'bx-'); hold on;
    errorbar(exp_statResp.U, exp_statResp.beta_y, exp_statResp.beta_y-[0; beta_y_wt_experimental_data_set_angle_5_lb(2:end)], [0; beta_y_wt_experimental_data_set_angle_5_ub(2:end)]-exp_statResp.beta_y, 0.3, 0.3, 'bx-'); hold on;
    errorbar(exp_Uf, exp_beta_yf, [], [], 0.3, 0.3, 'r^'); hold on;
    xlabel('U, [m/s]'); ylabel('\beta_y, [1/m]');
    legend({'Static response (model)', 'Hopf Bifurcations (model)', 'Static response (expr)', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'Hopf Bifurcations (expr.)'})
elseif ang_expr == 0.3
    plot(exp_statResp.U, exp_statResp.beta_y, 'bx-'); hold on;
    errorbar(exp_statResp.U, exp_statResp.beta_y, exp_statResp.beta_y-[0; beta_y_wt_experimental_data_set_angle_2_lb(2:end)], [0; beta_y_wt_experimental_data_set_angle_2_ub(2:end)]-exp_statResp.beta_y, 0.3, 0.3, 'bx-'); hold on;
    errorbar(exp_Uf, exp_beta_yf, [], [], 0.3, 0.3, 'r^'); hold on;
    xlabel('U, [m/s]'); ylabel('\beta_y, [1/m]');
    legend({'Static response (model)', 'Hopf Bifurcations (model)', '', 'Static response (expr)', '', '', '', '', '', '', '', '', '', '', '', 'Hopf Bifurcations (expr.)'})
elseif ang_expr == 0.6
    plot(exp_statResp.U, exp_statResp.beta_y, 'bx-'); hold on;
    errorbar(exp_statResp.U, exp_statResp.beta_y, exp_statResp.beta_y-[0; beta_y_wt_experimental_data_set_angle_3_lb(2:end)], [0; beta_y_wt_experimental_data_set_angle_3_ub(2:end)]-exp_statResp.beta_y, 0.3, 0.3, 'bx-'); hold on;
    errorbar(exp_Uf, exp_beta_yf, [], [], 0.3, 0.3, 'r^'); hold on;
    xlabel('U, [m/s]'); ylabel('\beta_y, [1/m]');
    legend({'Static response (model)', 'Hopf Bifurcations (model)', '', 'Static response (expr)', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'Hopf Bifurcations (expr.)'})
elseif ang_expr == 1.1
    plot(exp_statResp.U, exp_statResp.beta_y, 'bx-'); hold on;
    errorbar(exp_statResp.U, exp_statResp.beta_y, exp_statResp.beta_y-[0; beta_y_wt_experimental_data_set_angle_4_lb(2:end)], [0; beta_y_wt_experimental_data_set_angle_4_ub(2:end)]-exp_statResp.beta_y, 0.3, 0.3, 'bx-'); hold on;
    errorbar(exp_Uf, exp_beta_yf, [], [], 0.3, 0.3, 'r^'); hold on;
    xlabel('U, [m/s]'); ylabel('\beta_y, [1/m]');
    legend({'Static response (model)', 'Hopf Bifurcations (model)', '', 'Static response (expr)', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'Hopf Bifurcations (expr.)'})
else
    plot(exp_statResp.U, exp_statResp.beta_y, 'bx-'); hold on;
    plot(exp_Uf, exp_beta_yf, 'r^'); hold on;
    xlabel('U, [m/s]'); ylabel('\beta_y, [1/m]');
    legend({'Static response (model)', 'Hopf Bifurcations (model)', 'Static response (expr)', 'Hopf Bifurcations (expr.)'})
end    

%%
load('OMA_freq_lower_and_upper_bounds.mat');
load('OMA_damp_lower_and_upper_bounds.mat');
exp_mkrs = {'mx', 'ko', 'b^', 'rs'};
clrs = {'m', 'k', 'b', 'r'};
figure;
for mode=1:4
    subplot(2,1,1);
    for ii = 1:size(x_opt_lse, 1)
        plot(statResp(1, ii).U, statResp(1, ii).damp(mode,:), 'color', clrs{mode}); hold on;
    end
    plot(exp_statResp.U, exp_statResp.damp(mode,:), exp_mkrs{mode},...
            'markerFaceColor', clrs{mode}); hold on;
    if ang_expr == 1.4 & mode == 3
        errorbar(exp_statResp.U([1,5,6,7,8,9,10]), exp_statResp.damp(mode,[1,5,6,7,8,9,10]), exp_statResp.damp(mode,[1,5,6,7,8,9,10])-oop2_mode_angle_5_damp_lb, oop2_mode_angle_5_damp_ub-exp_statResp.damp(mode,[1,5,6,7,8,9,10]), 'bx'); hold on;
    elseif ang_expr == 1.4 & mode == 4
        errorbar(exp_statResp.U([1,5,6,7,8,9,10,11,12,18]), exp_statResp.damp(mode,[1,5,6,7,8,9,10,11,12,18]), exp_statResp.damp(mode,[1,5,6,7,8,9,10,11,12,18])-tor1_mode_angle_5_damp_lb, tor1_mode_angle_5_damp_ub-exp_statResp.damp(mode,[1,5,6,7,8,9,10,11,12,18]), 'rx'); hold on;
    end

    subplot(2,1,2);
    for ii = 1:size(x_opt_lse, 1)
        plot(statResp(1, ii).U, statResp(1, ii).frqs(mode,:), clrs{mode}); hold on;
    end
    plot(exp_statResp.U, exp_statResp.frqs(mode,:), exp_mkrs{mode},...
        'markerFaceColor', clrs{mode}); hold on;
    if ang_expr == 1.4 & mode == 3
        errorbar(exp_statResp.U([1,5,6,7,8,9,10]), exp_statResp.frqs(mode,[1,5,6,7,8,9,10]), exp_statResp.frqs(mode,[1,5,6,7,8,9,10])-oop2_mode_angle_5_freq_lb, oop2_mode_angle_5_freq_ub-exp_statResp.frqs(mode,[1,5,6,7,8,9,10]), 'bx'); hold on;
    elseif ang_expr == 1.4 & mode == 4
        errorbar(exp_statResp.U([1,5,6,7,8,9,10,11,12,18]), exp_statResp.frqs(mode,[1,5,6,7,8,9,10,11,12,18]), exp_statResp.frqs(mode,[1,5,6,7,8,9,10,11,12,18])-tor1_mode_angle_5_freq_lb, tor1_mode_angle_5_freq_ub-exp_statResp.frqs(mode,[1,5,6,7,8,9,10,11,12,18]), 'rx'); hold on;
    end
end
for ii = 1:size(x_opt_lse, 1)
    subplot(2,1,1); xline(Uf(ii, :)); xlabel('U, [m/s]'); ylabel('\zeta, [-]');
end
set(gca, 'xAxisLocation', 'origin')
for ii = 1:size(x_opt_lse, 1)
    subplot(2,1,2); xline(Uf(ii, :)); xlabel('U, [m/s]'); ylabel('\omega, [Hz]');
end
hold off

%%
load('OMA_freq_lower_and_upper_bounds.mat');
load('OMA_damp_lower_and_upper_bounds.mat');
exp_mkrs = {'mx', 'ko', 'b^', 'rs'};
clrs = {'m', 'k', 'b', 'r'};
figure;
for mode=1:4
    subplot(2,1,1);
    plot(statResp.U, statResp.damp(mode,:), 'color', clrs{mode}); hold on;
    plot(exp_statResp.U, exp_statResp.damp(mode,:), exp_mkrs{mode},...
        'markerFaceColor', clrs{mode}); hold on;
    if ang_expr == 1.4 & mode == 3
        errorbar(exp_statResp.U([1,5,6,7,8,9,10]), exp_statResp.damp(mode,[1,5,6,7,8,9,10]), exp_statResp.damp(mode,[1,5,6,7,8,9,10])-oop2_mode_angle_5_damp_lb, oop2_mode_angle_5_damp_ub-exp_statResp.damp(mode,[1,5,6,7,8,9,10]), 'bx'); hold on;
    elseif ang_expr == 1.4 & mode == 4
        errorbar(exp_statResp.U([1,5,6,7,8,9,10,11,12,18]), exp_statResp.damp(mode,[1,5,6,7,8,9,10,11,12,18]), exp_statResp.damp(mode,[1,5,6,7,8,9,10,11,12,18])-tor1_mode_angle_5_damp_lb, tor1_mode_angle_5_damp_ub-exp_statResp.damp(mode,[1,5,6,7,8,9,10,11,12,18]), 'rx'); hold on;
    end

    subplot(2,1,2);
    plot(statResp.U, statResp.frqs(mode,:), clrs{mode}); hold on;
    plot(exp_statResp.U, exp_statResp.frqs(mode,:), exp_mkrs{mode},...
        'markerFaceColor', clrs{mode}); hold on;
    if ang_expr == 1.4 & mode == 3
        errorbar(exp_statResp.U([1,5,6,7,8,9,10]), exp_statResp.frqs(mode,[1,5,6,7,8,9,10]), exp_statResp.frqs(mode,[1,5,6,7,8,9,10])-oop2_mode_angle_5_freq_lb, oop2_mode_angle_5_freq_ub-exp_statResp.frqs(mode,[1,5,6,7,8,9,10]), 'bx'); hold on;
    elseif ang_expr == 1.4 & mode == 4
        errorbar(exp_statResp.U([1,5,6,7,8,9,10,11,12,18]), exp_statResp.frqs(mode,[1,5,6,7,8,9,10,11,12,18]), exp_statResp.frqs(mode,[1,5,6,7,8,9,10,11,12,18])-tor1_mode_angle_5_freq_lb, tor1_mode_angle_5_freq_ub-exp_statResp.frqs(mode,[1,5,6,7,8,9,10,11,12,18]), 'rx'); hold on;
    end
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
