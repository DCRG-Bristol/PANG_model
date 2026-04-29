%% WT exp.angle 1.4, speed 0m/s, damp TOR1
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_1 = modalData.stabPls{5, 1};
figure; plot(modal_data_idx_angle_5_idx_speed_1, 'k.')

r_min = 16;
r_max = 17;

tor1_mode_angle_5_idx_speed_1 = modal_data_idx_angle_5_idx_speed_1(abs(modal_data_idx_angle_5_idx_speed_1) >= r_min & abs(modal_data_idx_angle_5_idx_speed_1) <= r_max);
tor1_mode_freq_angle_5_idx_speed_1 = abs(tor1_mode_angle_5_idx_speed_1);
tor1_mode_damp_angle_5_idx_speed_1 = -real(tor1_mode_angle_5_idx_speed_1)./abs(tor1_mode_angle_5_idx_speed_1);
tor1_mode_damp_freq_angle_5_idx_speed_1 = [tor1_mode_damp_angle_5_idx_speed_1(1:end), tor1_mode_freq_angle_5_idx_speed_1(1:end)];

mu_tor1_mode_damp_freq_angle_5_idx_speed_1 = mean(tor1_mode_damp_freq_angle_5_idx_speed_1, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_5_idx_speed_1 = cov(tor1_mode_damp_freq_angle_5_idx_speed_1);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_5_idx_speed_1 = mvnpdf(tor1_mode_damp_freq_angle_5_idx_speed_1, mu_tor1_mode_damp_freq_angle_5_idx_speed_1, Sigma_tor1_mode_damp_freq_angle_5_idx_speed_1);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_5_idx_speed_1(:,1), tor1_mode_damp_freq_angle_5_idx_speed_1(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_5_idx_speed_1);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_5_idx_speed_1 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_5_idx_speed_1;
    plot(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_1(:,1), ellipse_tor1_mode_damp_freq_angle_5_idx_speed_1(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_5_idx_speed_1_lb = min(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_1(:, 2));
tor1_mode_freq_angle_5_idx_speed_1_ub = max(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_1(:, 2));
tor1_mode_damp_angle_5_idx_speed_1_lb = min(tor1_mode_damp_freq_angle_5_idx_speed_1(:, 1));
tor1_mode_damp_angle_5_idx_speed_1_ub_ = sort(tor1_mode_damp_freq_angle_5_idx_speed_1(:, 1), 'descend');
tor1_mode_damp_angle_5_idx_speed_1_ub = tor1_mode_damp_angle_5_idx_speed_1_ub_(6);

%% WT exp.angle 1.4, speed 0m/s, damp OOP2
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_1 = modalData.stabPls{5, 1};
figure; plot(modal_data_idx_angle_5_idx_speed_1, 'k.')

r_min = 11;
r_max = 12;

oop2_mode_angle_5_idx_speed_1 = modal_data_idx_angle_5_idx_speed_1(abs(modal_data_idx_angle_5_idx_speed_1) >= r_min & abs(modal_data_idx_angle_5_idx_speed_1) <= r_max);
oop2_mode_freq_angle_5_idx_speed_1 = abs(oop2_mode_angle_5_idx_speed_1);
oop2_mode_damp_angle_5_idx_speed_1 = -real(oop2_mode_angle_5_idx_speed_1)./abs(oop2_mode_angle_5_idx_speed_1);
oop2_mode_damp_freq_angle_5_idx_speed_1 = [oop2_mode_damp_angle_5_idx_speed_1(1:end), oop2_mode_freq_angle_5_idx_speed_1(1:end)];

mu_oop2_mode_damp_freq_angle_5_idx_speed_1 = mean(oop2_mode_damp_freq_angle_5_idx_speed_1, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_5_idx_speed_1 = cov(oop2_mode_damp_freq_angle_5_idx_speed_1);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_5_idx_speed_1 = mvnpdf(oop2_mode_damp_freq_angle_5_idx_speed_1, mu_oop2_mode_damp_freq_angle_5_idx_speed_1, Sigma_oop2_mode_damp_freq_angle_5_idx_speed_1);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_5_idx_speed_1(:,1), oop2_mode_damp_freq_angle_5_idx_speed_1(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_5_idx_speed_1);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_5_idx_speed_1 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_5_idx_speed_1;
    plot(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_1(:,1), ellipse_oop2_mode_damp_freq_angle_5_idx_speed_1(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_5_idx_speed_1_lb = min(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_1(:, 2));
oop2_mode_freq_angle_5_idx_speed_1_ub = max(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_1(:, 2));
oop2_mode_damp_angle_5_idx_speed_1_lb = min(oop2_mode_damp_freq_angle_5_idx_speed_1(:, 1));
oop2_mode_damp_angle_5_idx_speed_1_ub = max(oop2_mode_damp_freq_angle_5_idx_speed_1(:, 1));

%% WT exp.angle 1.4, speed 12.5m/s, damp TOR1
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_5 = modalData.stabPls{5, 5};
figure; plot(modal_data_idx_angle_5_idx_speed_5, 'k.')

r_min = 15;
r_max = 16;

tor1_mode_angle_5_idx_speed_5 = modal_data_idx_angle_5_idx_speed_5(abs(modal_data_idx_angle_5_idx_speed_5) >= r_min & abs(modal_data_idx_angle_5_idx_speed_5) <= r_max);
tor1_mode_freq_angle_5_idx_speed_5 = abs(tor1_mode_angle_5_idx_speed_5);
tor1_mode_damp_angle_5_idx_speed_5 = -real(tor1_mode_angle_5_idx_speed_5)./abs(tor1_mode_angle_5_idx_speed_5);
tor1_mode_damp_freq_angle_5_idx_speed_5 = [tor1_mode_damp_angle_5_idx_speed_5(1:end), tor1_mode_freq_angle_5_idx_speed_5(1:end)];

mu_tor1_mode_damp_freq_angle_5_idx_speed_5 = mean(tor1_mode_damp_freq_angle_5_idx_speed_5, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_5_idx_speed_5 = cov(tor1_mode_damp_freq_angle_5_idx_speed_5);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_5_idx_speed_5 = mvnpdf(tor1_mode_damp_freq_angle_5_idx_speed_5, mu_tor1_mode_damp_freq_angle_5_idx_speed_5, Sigma_tor1_mode_damp_freq_angle_5_idx_speed_5);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_5_idx_speed_5(:,1), tor1_mode_damp_freq_angle_5_idx_speed_5(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_5_idx_speed_5);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_5_idx_speed_5 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_5_idx_speed_5;
    plot(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_5(:,1), ellipse_tor1_mode_damp_freq_angle_5_idx_speed_5(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_5_idx_speed_5_lb = min(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_5(:, 2));
tor1_mode_freq_angle_5_idx_speed_5_ub = max(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_5(:, 2));
tor1_mode_damp_angle_5_idx_speed_5_lb = min(tor1_mode_damp_freq_angle_5_idx_speed_5(:, 1));
tor1_mode_damp_angle_5_idx_speed_5_ub = max(tor1_mode_damp_freq_angle_5_idx_speed_5(:, 1));

%% WT exp.angle 1.4, speed 12.5m/s, damp OOP2
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_5 = modalData.stabPls{5, 5};
figure; plot(modal_data_idx_angle_5_idx_speed_5, 'k.')

r_min = 11;
r_max = 12;

oop2_mode_angle_5_idx_speed_5 = modal_data_idx_angle_5_idx_speed_5(abs(modal_data_idx_angle_5_idx_speed_5) >= r_min & abs(modal_data_idx_angle_5_idx_speed_5) <= r_max);
oop2_mode_freq_angle_5_idx_speed_5 = abs(oop2_mode_angle_5_idx_speed_5);
oop2_mode_damp_angle_5_idx_speed_5 = -real(oop2_mode_angle_5_idx_speed_5)./abs(oop2_mode_angle_5_idx_speed_5);
oop2_mode_damp_freq_angle_5_idx_speed_5 = [oop2_mode_damp_angle_5_idx_speed_5(1:end), oop2_mode_freq_angle_5_idx_speed_5(1:end)];

mu_oop2_mode_damp_freq_angle_5_idx_speed_5 = mean(oop2_mode_damp_freq_angle_5_idx_speed_5, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_5_idx_speed_5 = cov(oop2_mode_damp_freq_angle_5_idx_speed_5);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_5_idx_speed_5 = mvnpdf(oop2_mode_damp_freq_angle_5_idx_speed_5, mu_oop2_mode_damp_freq_angle_5_idx_speed_5, Sigma_oop2_mode_damp_freq_angle_5_idx_speed_5);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_5_idx_speed_5(:,1), oop2_mode_damp_freq_angle_5_idx_speed_5(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_5_idx_speed_5);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_5_idx_speed_5 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_5_idx_speed_5;
    plot(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_5(:,1), ellipse_oop2_mode_damp_freq_angle_5_idx_speed_5(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_5_idx_speed_5_lb = min(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_5(:, 2));
oop2_mode_freq_angle_5_idx_speed_5_ub = max(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_5(:, 2));
oop2_mode_damp_angle_5_idx_speed_5_lb = min(oop2_mode_damp_freq_angle_5_idx_speed_5(:, 1));
oop2_mode_damp_angle_5_idx_speed_5_ub_ = sort(oop2_mode_damp_freq_angle_5_idx_speed_5(:, 1), 'descend');
oop2_mode_damp_angle_5_idx_speed_5_ub = oop2_mode_damp_angle_5_idx_speed_5_ub_(3);

%% WT exp.angle 1.4, speed 15m/s, damp TOR1
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_6 = modalData.stabPls{5, 6};
figure; plot(modal_data_idx_angle_5_idx_speed_6, 'k.')

r_min = 15;
r_max = 15.5;

tor1_mode_angle_5_idx_speed_6 = modal_data_idx_angle_5_idx_speed_6(abs(modal_data_idx_angle_5_idx_speed_6) >= r_min & abs(modal_data_idx_angle_5_idx_speed_6) <= r_max);
tor1_mode_freq_angle_5_idx_speed_6 = abs(tor1_mode_angle_5_idx_speed_6);
tor1_mode_damp_angle_5_idx_speed_6 = -real(tor1_mode_angle_5_idx_speed_6)./abs(tor1_mode_angle_5_idx_speed_6);
tor1_mode_damp_freq_angle_5_idx_speed_6 = [tor1_mode_damp_angle_5_idx_speed_6(1:end), tor1_mode_freq_angle_5_idx_speed_6(1:end)];

mu_tor1_mode_damp_freq_angle_5_idx_speed_6 = mean(tor1_mode_damp_freq_angle_5_idx_speed_6, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_5_idx_speed_6 = cov(tor1_mode_damp_freq_angle_5_idx_speed_6);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_5_idx_speed_6 = mvnpdf(tor1_mode_damp_freq_angle_5_idx_speed_6, mu_tor1_mode_damp_freq_angle_5_idx_speed_6, Sigma_tor1_mode_damp_freq_angle_5_idx_speed_6);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_5_idx_speed_6(:,1), tor1_mode_damp_freq_angle_5_idx_speed_6(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_5_idx_speed_6);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_5_idx_speed_6 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_5_idx_speed_6;
    plot(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_6(:,1), ellipse_tor1_mode_damp_freq_angle_5_idx_speed_6(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_5_idx_speed_6_lb = min(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_6(:, 2));
tor1_mode_freq_angle_5_idx_speed_6_ub = max(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_6(:, 2));
tor1_mode_damp_angle_5_idx_speed_6_lb = min(tor1_mode_damp_freq_angle_5_idx_speed_6(:, 1));
tor1_mode_damp_angle_5_idx_speed_6_ub = max(tor1_mode_damp_freq_angle_5_idx_speed_6(:, 1));

%% WT exp.angle 1.4, speed 15m/s, damp OOP2
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_6 = modalData.stabPls{5, 6};
figure; plot(modal_data_idx_angle_5_idx_speed_6, 'k.')

r_min = 11;
r_max = 12;

oop2_mode_angle_5_idx_speed_6 = modal_data_idx_angle_5_idx_speed_6(abs(modal_data_idx_angle_5_idx_speed_6) >= r_min & abs(modal_data_idx_angle_5_idx_speed_6) <= r_max);
oop2_mode_freq_angle_5_idx_speed_6 = abs(oop2_mode_angle_5_idx_speed_6);
oop2_mode_damp_angle_5_idx_speed_6 = -real(oop2_mode_angle_5_idx_speed_6)./abs(oop2_mode_angle_5_idx_speed_6);
oop2_mode_damp_freq_angle_5_idx_speed_6 = [oop2_mode_damp_angle_5_idx_speed_6(1:end), oop2_mode_freq_angle_5_idx_speed_6(1:end)];

mu_oop2_mode_damp_freq_angle_5_idx_speed_6 = mean(oop2_mode_damp_freq_angle_5_idx_speed_6, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_5_idx_speed_6 = cov(oop2_mode_damp_freq_angle_5_idx_speed_6);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_5_idx_speed_6 = mvnpdf(oop2_mode_damp_freq_angle_5_idx_speed_6, mu_oop2_mode_damp_freq_angle_5_idx_speed_6, Sigma_oop2_mode_damp_freq_angle_5_idx_speed_6);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_5_idx_speed_6(:,1), oop2_mode_damp_freq_angle_5_idx_speed_6(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_5_idx_speed_6);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_5_idx_speed_6 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_5_idx_speed_6;
    plot(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_6(:,1), ellipse_oop2_mode_damp_freq_angle_5_idx_speed_6(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_5_idx_speed_6_lb = min(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_6(:, 2));
oop2_mode_freq_angle_5_idx_speed_6_ub = max(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_6(:, 2));
oop2_mode_damp_angle_5_idx_speed_6_lb = min(oop2_mode_damp_freq_angle_5_idx_speed_6(:, 1));
oop2_mode_damp_angle_5_idx_speed_6_ub = max(oop2_mode_damp_freq_angle_5_idx_speed_6(:, 1));

%% WT exp.angle 1.4, speed 17m/s, damp TOR1
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_7 = modalData.stabPls{5, 7};
figure; plot(modal_data_idx_angle_5_idx_speed_7, 'k.')

r_min = 14.5;
r_max = 15;

tor1_mode_angle_5_idx_speed_7 = modal_data_idx_angle_5_idx_speed_7(abs(modal_data_idx_angle_5_idx_speed_7) >= r_min & abs(modal_data_idx_angle_5_idx_speed_7) <= r_max);
tor1_mode_freq_angle_5_idx_speed_7 = abs(tor1_mode_angle_5_idx_speed_7);
tor1_mode_damp_angle_5_idx_speed_7 = -real(tor1_mode_angle_5_idx_speed_7)./abs(tor1_mode_angle_5_idx_speed_7);
tor1_mode_damp_freq_angle_5_idx_speed_7 = [tor1_mode_damp_angle_5_idx_speed_7(1:end), tor1_mode_freq_angle_5_idx_speed_7(1:end)];

mu_tor1_mode_damp_freq_angle_5_idx_speed_7 = mean(tor1_mode_damp_freq_angle_5_idx_speed_7, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_5_idx_speed_7 = cov(tor1_mode_damp_freq_angle_5_idx_speed_7);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_5_idx_speed_7 = mvnpdf(tor1_mode_damp_freq_angle_5_idx_speed_7, mu_tor1_mode_damp_freq_angle_5_idx_speed_7, Sigma_tor1_mode_damp_freq_angle_5_idx_speed_7);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_5_idx_speed_7(:,1), tor1_mode_damp_freq_angle_5_idx_speed_7(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_5_idx_speed_7);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_5_idx_speed_7 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_5_idx_speed_7;
    plot(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_7(:,1), ellipse_tor1_mode_damp_freq_angle_5_idx_speed_7(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_5_idx_speed_7_lb = min(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_7(:, 2));
tor1_mode_freq_angle_5_idx_speed_7_ub = max(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_7(:, 2));
tor1_mode_damp_angle_5_idx_speed_7_lb = min(tor1_mode_damp_freq_angle_5_idx_speed_7(:, 1));
tor1_mode_damp_angle_5_idx_speed_7_ub_ = sort(tor1_mode_damp_freq_angle_5_idx_speed_7(:, 1), 'descend');
tor1_mode_damp_angle_5_idx_speed_7_ub = tor1_mode_damp_angle_5_idx_speed_7_ub_(2);

%% WT exp.angle 1.4, speed 17m/s, damp OOP2
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_7 = modalData.stabPls{5, 7};
figure; plot(modal_data_idx_angle_5_idx_speed_7, 'k.')

r_min = 11;
r_max = 12;

oop2_mode_angle_5_idx_speed_7 = modal_data_idx_angle_5_idx_speed_7(abs(modal_data_idx_angle_5_idx_speed_7) >= r_min & abs(modal_data_idx_angle_5_idx_speed_7) <= r_max);
oop2_mode_freq_angle_5_idx_speed_7 = abs(oop2_mode_angle_5_idx_speed_7);
oop2_mode_damp_angle_5_idx_speed_7 = -real(oop2_mode_angle_5_idx_speed_7)./abs(oop2_mode_angle_5_idx_speed_7);
oop2_mode_damp_freq_angle_5_idx_speed_7 = [oop2_mode_damp_angle_5_idx_speed_7(1:end), oop2_mode_freq_angle_5_idx_speed_7(1:end)];

mu_oop2_mode_damp_freq_angle_5_idx_speed_7 = mean(oop2_mode_damp_freq_angle_5_idx_speed_7, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_5_idx_speed_7 = cov(oop2_mode_damp_freq_angle_5_idx_speed_7);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_5_idx_speed_7 = mvnpdf(oop2_mode_damp_freq_angle_5_idx_speed_7, mu_oop2_mode_damp_freq_angle_5_idx_speed_7, Sigma_oop2_mode_damp_freq_angle_5_idx_speed_7);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_5_idx_speed_7(:,1), oop2_mode_damp_freq_angle_5_idx_speed_7(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_5_idx_speed_7);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_5_idx_speed_7 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_5_idx_speed_7;
    plot(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_7(:,1), ellipse_oop2_mode_damp_freq_angle_5_idx_speed_7(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_5_idx_speed_7_lb = min(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_7(:, 2));
oop2_mode_freq_angle_5_idx_speed_7_ub = max(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_7(:, 2));
oop2_mode_damp_angle_5_idx_speed_7_lb = min(oop2_mode_damp_freq_angle_5_idx_speed_7(:, 1));
oop2_mode_damp_angle_5_idx_speed_7_ub_ = sort(oop2_mode_damp_freq_angle_5_idx_speed_7(:, 1), 'descend');
oop2_mode_damp_angle_5_idx_speed_7_ub = oop2_mode_damp_angle_5_idx_speed_7_ub_(8);

%% WT exp.angle 1.4, speed 19m/s, damp TOR1
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_8 = modalData.stabPls{5, 8};
figure; plot(modal_data_idx_angle_5_idx_speed_8, 'k.')

r_min = 14;
r_max = 14.5;

tor1_mode_angle_5_idx_speed_8 = modal_data_idx_angle_5_idx_speed_8(abs(modal_data_idx_angle_5_idx_speed_8) >= r_min & abs(modal_data_idx_angle_5_idx_speed_8) <= r_max);
tor1_mode_freq_angle_5_idx_speed_8 = abs(tor1_mode_angle_5_idx_speed_8);
tor1_mode_damp_angle_5_idx_speed_8 = -real(tor1_mode_angle_5_idx_speed_8)./abs(tor1_mode_angle_5_idx_speed_8);
tor1_mode_damp_freq_angle_5_idx_speed_8 = [tor1_mode_damp_angle_5_idx_speed_8(1:end), tor1_mode_freq_angle_5_idx_speed_8(1:end)];

mu_tor1_mode_damp_freq_angle_5_idx_speed_8 = mean(tor1_mode_damp_freq_angle_5_idx_speed_8, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_5_idx_speed_8 = cov(tor1_mode_damp_freq_angle_5_idx_speed_8);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_5_idx_speed_8 = mvnpdf(tor1_mode_damp_freq_angle_5_idx_speed_8, mu_tor1_mode_damp_freq_angle_5_idx_speed_8, Sigma_tor1_mode_damp_freq_angle_5_idx_speed_8);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_5_idx_speed_8(:,1), tor1_mode_damp_freq_angle_5_idx_speed_8(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_5_idx_speed_8);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_5_idx_speed_8 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_5_idx_speed_8;
    plot(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_8(:,1), ellipse_tor1_mode_damp_freq_angle_5_idx_speed_8(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_5_idx_speed_8_lb = min(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_8(:, 2));
tor1_mode_freq_angle_5_idx_speed_8_ub = max(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_8(:, 2));
tor1_mode_damp_angle_5_idx_speed_8_lb = min(tor1_mode_damp_freq_angle_5_idx_speed_8(:, 1));
tor1_mode_damp_angle_5_idx_speed_8_ub = max(tor1_mode_damp_freq_angle_5_idx_speed_8(:, 1));

%% WT exp.angle 1.4, speed 19m/s, damp OOP2
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_8 = modalData.stabPls{5, 8};
figure; plot(modal_data_idx_angle_5_idx_speed_8, 'k.')

r_min = 11;
r_max = 12;

oop2_mode_angle_5_idx_speed_8 = modal_data_idx_angle_5_idx_speed_8(abs(modal_data_idx_angle_5_idx_speed_8) >= r_min & abs(modal_data_idx_angle_5_idx_speed_8) <= r_max);
oop2_mode_freq_angle_5_idx_speed_8 = abs(oop2_mode_angle_5_idx_speed_8);
oop2_mode_damp_angle_5_idx_speed_8 = -real(oop2_mode_angle_5_idx_speed_8)./abs(oop2_mode_angle_5_idx_speed_8);
oop2_mode_damp_freq_angle_5_idx_speed_8 = [oop2_mode_damp_angle_5_idx_speed_8(1:end), oop2_mode_freq_angle_5_idx_speed_8(1:end)];

mu_oop2_mode_damp_freq_angle_5_idx_speed_8 = mean(oop2_mode_damp_freq_angle_5_idx_speed_8, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_5_idx_speed_8 = cov(oop2_mode_damp_freq_angle_5_idx_speed_8);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_5_idx_speed_8 = mvnpdf(oop2_mode_damp_freq_angle_5_idx_speed_8, mu_oop2_mode_damp_freq_angle_5_idx_speed_8, Sigma_oop2_mode_damp_freq_angle_5_idx_speed_8);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_5_idx_speed_8(:,1), oop2_mode_damp_freq_angle_5_idx_speed_8(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_5_idx_speed_8);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_5_idx_speed_8 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_5_idx_speed_8;
    plot(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_8(:,1), ellipse_oop2_mode_damp_freq_angle_5_idx_speed_8(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_5_idx_speed_8_lb = min(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_8(:, 2));
oop2_mode_freq_angle_5_idx_speed_8_ub = max(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_8(:, 2));
oop2_mode_damp_angle_5_idx_speed_8_lb_ = sort(oop2_mode_damp_freq_angle_5_idx_speed_8(:, 1), 'ascend');
oop2_mode_damp_angle_5_idx_speed_8_lb = oop2_mode_damp_angle_5_idx_speed_8_lb_(4);
oop2_mode_damp_angle_5_idx_speed_8_ub_ = sort(oop2_mode_damp_freq_angle_5_idx_speed_8(:, 1), 'descend');
oop2_mode_damp_angle_5_idx_speed_8_ub = oop2_mode_damp_angle_5_idx_speed_8_ub_(4);

%% WT exp.angle 1.4, speed 21m/s, damp TOR1
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_9 = modalData.stabPls{5, 9};
figure; plot(modal_data_idx_angle_5_idx_speed_9, 'k.')

r_min = 13;
r_max = 14;

tor1_mode_angle_5_idx_speed_9 = modal_data_idx_angle_5_idx_speed_9(abs(modal_data_idx_angle_5_idx_speed_9) >= r_min & abs(modal_data_idx_angle_5_idx_speed_9) <= r_max);
tor1_mode_freq_angle_5_idx_speed_9 = abs(tor1_mode_angle_5_idx_speed_9);
tor1_mode_damp_angle_5_idx_speed_9 = -real(tor1_mode_angle_5_idx_speed_9)./abs(tor1_mode_angle_5_idx_speed_9);
tor1_mode_damp_freq_angle_5_idx_speed_9 = [tor1_mode_damp_angle_5_idx_speed_9(1:end), tor1_mode_freq_angle_5_idx_speed_9(1:end)];

mu_tor1_mode_damp_freq_angle_5_idx_speed_9 = mean(tor1_mode_damp_freq_angle_5_idx_speed_9, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_5_idx_speed_9 = cov(tor1_mode_damp_freq_angle_5_idx_speed_9);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_5_idx_speed_9 = mvnpdf(tor1_mode_damp_freq_angle_5_idx_speed_9, mu_tor1_mode_damp_freq_angle_5_idx_speed_9, Sigma_tor1_mode_damp_freq_angle_5_idx_speed_9);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_5_idx_speed_9(:,1), tor1_mode_damp_freq_angle_5_idx_speed_9(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_5_idx_speed_9);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_5_idx_speed_9 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_5_idx_speed_9;
    plot(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_9(:,1), ellipse_tor1_mode_damp_freq_angle_5_idx_speed_9(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_5_idx_speed_9_lb = min(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_9(:, 2));
tor1_mode_freq_angle_5_idx_speed_9_ub = max(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_9(:, 2));
tor1_mode_damp_angle_5_idx_speed_9_lb_ = sort(tor1_mode_damp_freq_angle_5_idx_speed_9(:, 1), 'ascend');
tor1_mode_damp_angle_5_idx_speed_9_lb = tor1_mode_damp_angle_5_idx_speed_9_lb_(1);
tor1_mode_damp_angle_5_idx_speed_9_ub_ = sort(tor1_mode_damp_freq_angle_5_idx_speed_9(:, 1), 'descend');
tor1_mode_damp_angle_5_idx_speed_9_ub = tor1_mode_damp_angle_5_idx_speed_9_ub_(4);

%% WT exp.angle 1.4, speed 21m/s, damp OOP2
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_9 = modalData.stabPls{5, 9};
figure; plot(modal_data_idx_angle_5_idx_speed_9, 'k.')

r_min = 11;
r_max = 12;

oop2_mode_angle_5_idx_speed_9 = modal_data_idx_angle_5_idx_speed_9(abs(modal_data_idx_angle_5_idx_speed_9) >= r_min & abs(modal_data_idx_angle_5_idx_speed_9) <= r_max);
oop2_mode_freq_angle_5_idx_speed_9 = abs(oop2_mode_angle_5_idx_speed_9);
oop2_mode_damp_angle_5_idx_speed_9 = -real(oop2_mode_angle_5_idx_speed_9)./abs(oop2_mode_angle_5_idx_speed_9);
oop2_mode_damp_freq_angle_5_idx_speed_9 = [oop2_mode_damp_angle_5_idx_speed_9(1:end), oop2_mode_freq_angle_5_idx_speed_9(1:end)];

mu_oop2_mode_damp_freq_angle_5_idx_speed_9 = mean(oop2_mode_damp_freq_angle_5_idx_speed_9, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_5_idx_speed_9 = cov(oop2_mode_damp_freq_angle_5_idx_speed_9);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_5_idx_speed_9 = mvnpdf(oop2_mode_damp_freq_angle_5_idx_speed_9, mu_oop2_mode_damp_freq_angle_5_idx_speed_9, Sigma_oop2_mode_damp_freq_angle_5_idx_speed_9);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_5_idx_speed_9(:,1), oop2_mode_damp_freq_angle_5_idx_speed_9(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_5_idx_speed_9);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_5_idx_speed_9 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_5_idx_speed_9;
    plot(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_9(:,1), ellipse_oop2_mode_damp_freq_angle_5_idx_speed_9(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_5_idx_speed_9_lb = min(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_9(:, 2));
oop2_mode_freq_angle_5_idx_speed_9_ub = max(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_9(:, 2));
oop2_mode_damp_angle_5_idx_speed_9_lb_ = sort(oop2_mode_damp_freq_angle_5_idx_speed_9(:, 1), 'ascend');
oop2_mode_damp_angle_5_idx_speed_9_lb = oop2_mode_damp_angle_5_idx_speed_9_lb_(2);
oop2_mode_damp_angle_5_idx_speed_9_ub_ = sort(oop2_mode_damp_freq_angle_5_idx_speed_9(:, 1), 'descend');
oop2_mode_damp_angle_5_idx_speed_9_ub = oop2_mode_damp_angle_5_idx_speed_9_ub_(8);

%% WT exp.angle 1.4, speed 22.5m/s, damp TOR1
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_10 = modalData.stabPls{5, 10};
figure; plot(modal_data_idx_angle_5_idx_speed_10, 'k.')

r_min = 12.5;
r_max = 13.5;

tor1_mode_angle_5_idx_speed_10 = modal_data_idx_angle_5_idx_speed_10(abs(modal_data_idx_angle_5_idx_speed_10) >= r_min & abs(modal_data_idx_angle_5_idx_speed_10) <= r_max);
tor1_mode_freq_angle_5_idx_speed_10 = abs(tor1_mode_angle_5_idx_speed_10);
tor1_mode_damp_angle_5_idx_speed_10 = -real(tor1_mode_angle_5_idx_speed_10)./abs(tor1_mode_angle_5_idx_speed_10);
tor1_mode_damp_freq_angle_5_idx_speed_10 = [tor1_mode_damp_angle_5_idx_speed_10(1:end), tor1_mode_freq_angle_5_idx_speed_10(1:end)];

mu_tor1_mode_damp_freq_angle_5_idx_speed_10 = mean(tor1_mode_damp_freq_angle_5_idx_speed_10, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_5_idx_speed_10 = cov(tor1_mode_damp_freq_angle_5_idx_speed_10);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_5_idx_speed_10 = mvnpdf(tor1_mode_damp_freq_angle_5_idx_speed_10, mu_tor1_mode_damp_freq_angle_5_idx_speed_10, Sigma_tor1_mode_damp_freq_angle_5_idx_speed_10);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_5_idx_speed_10(:,1), tor1_mode_damp_freq_angle_5_idx_speed_10(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_5_idx_speed_10);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_5_idx_speed_10 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_5_idx_speed_10;
    plot(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_10(:,1), ellipse_tor1_mode_damp_freq_angle_5_idx_speed_10(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_5_idx_speed_10_lb = min(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_10(:, 2));
tor1_mode_freq_angle_5_idx_speed_10_ub = max(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_10(:, 2));
tor1_mode_damp_angle_5_idx_speed_10_lb = min(tor1_mode_damp_freq_angle_5_idx_speed_10(:, 1));
tor1_mode_damp_angle_5_idx_speed_10_ub_ = sort(tor1_mode_damp_freq_angle_5_idx_speed_10(:, 1), 'descend');
tor1_mode_damp_angle_5_idx_speed_10_ub = tor1_mode_damp_angle_5_idx_speed_10_ub_(19);

%% WT exp.angle 1.4, speed 22.5m/s, damp OOP2
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_10 = modalData.stabPls{5, 10};
figure; plot(modal_data_idx_angle_5_idx_speed_10, 'k.')

r_min = 12.5;
r_max = 12.8;

oop2_mode_angle_5_idx_speed_10 = modal_data_idx_angle_5_idx_speed_10(abs(modal_data_idx_angle_5_idx_speed_10) >= r_min & abs(modal_data_idx_angle_5_idx_speed_10) <= r_max);
oop2_mode_freq_angle_5_idx_speed_10 = abs(oop2_mode_angle_5_idx_speed_10);
oop2_mode_damp_angle_5_idx_speed_10 = -real(oop2_mode_angle_5_idx_speed_10)./abs(oop2_mode_angle_5_idx_speed_10);
oop2_mode_damp_freq_angle_5_idx_speed_10 = [oop2_mode_damp_angle_5_idx_speed_10(1:end), oop2_mode_freq_angle_5_idx_speed_10(1:end)];

mu_oop2_mode_damp_freq_angle_5_idx_speed_10 = mean(oop2_mode_damp_freq_angle_5_idx_speed_10, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_5_idx_speed_10 = cov(oop2_mode_damp_freq_angle_5_idx_speed_10);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_5_idx_speed_10 = mvnpdf(oop2_mode_damp_freq_angle_5_idx_speed_10, mu_oop2_mode_damp_freq_angle_5_idx_speed_10, Sigma_oop2_mode_damp_freq_angle_5_idx_speed_10);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_5_idx_speed_10(:,1), oop2_mode_damp_freq_angle_5_idx_speed_10(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_5_idx_speed_10);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_5_idx_speed_10 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_5_idx_speed_10;
    plot(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_10(:,1), ellipse_oop2_mode_damp_freq_angle_5_idx_speed_10(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_5_idx_speed_10_lb = min(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_10(:, 2));
oop2_mode_freq_angle_5_idx_speed_10_ub = max(ellipse_oop2_mode_damp_freq_angle_5_idx_speed_10(:, 2));
oop2_mode_damp_angle_5_idx_speed_10_lb = min(oop2_mode_damp_freq_angle_5_idx_speed_10(:, 1));
oop2_mode_damp_angle_5_idx_speed_10_ub = max(oop2_mode_damp_freq_angle_5_idx_speed_10(:, 1));

%% WT exp.angle 1.4, speed 24m/s, damp TOR1
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_11 = modalData.stabPls{5, 11};
figure; plot(modal_data_idx_angle_5_idx_speed_11, 'k.')

r_min = 12;
r_max = 13;

tor1_mode_angle_5_idx_speed_11 = modal_data_idx_angle_5_idx_speed_11(abs(modal_data_idx_angle_5_idx_speed_11) >= r_min & abs(modal_data_idx_angle_5_idx_speed_11) <= r_max);
tor1_mode_freq_angle_5_idx_speed_11 = abs(tor1_mode_angle_5_idx_speed_11);
tor1_mode_damp_angle_5_idx_speed_11 = -real(tor1_mode_angle_5_idx_speed_11)./abs(tor1_mode_angle_5_idx_speed_11);
tor1_mode_damp_freq_angle_5_idx_speed_11 = [tor1_mode_damp_angle_5_idx_speed_11(1:end), tor1_mode_freq_angle_5_idx_speed_11(1:end)];

mu_tor1_mode_damp_freq_angle_5_idx_speed_11 = mean(tor1_mode_damp_freq_angle_5_idx_speed_11, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_5_idx_speed_11 = cov(tor1_mode_damp_freq_angle_5_idx_speed_11);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_5_idx_speed_11 = mvnpdf(tor1_mode_damp_freq_angle_5_idx_speed_11, mu_tor1_mode_damp_freq_angle_5_idx_speed_11, Sigma_tor1_mode_damp_freq_angle_5_idx_speed_11);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_5_idx_speed_11(:,1), tor1_mode_damp_freq_angle_5_idx_speed_11(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_5_idx_speed_11);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_5_idx_speed_11 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_5_idx_speed_11;
    plot(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_11(:,1), ellipse_tor1_mode_damp_freq_angle_5_idx_speed_11(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_5_idx_speed_11_lb = min(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_11(:, 2));
tor1_mode_freq_angle_5_idx_speed_11_ub = max(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_11(:, 2));
tor1_mode_damp_angle_5_idx_speed_11_lb = min(tor1_mode_damp_freq_angle_5_idx_speed_11(:, 1));
tor1_mode_damp_angle_5_idx_speed_11_ub = max(tor1_mode_damp_freq_angle_5_idx_speed_11(:, 1));

%% WT exp.angle 1.4, speed 24.7m/s, damp TOR1
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_12 = modalData.stabPls{5, 12};
figure; plot(modal_data_idx_angle_5_idx_speed_12, 'k.')

r_min = 12;
r_max = 13;

tor1_mode_angle_5_idx_speed_12 = modal_data_idx_angle_5_idx_speed_12(abs(modal_data_idx_angle_5_idx_speed_12) >= r_min & abs(modal_data_idx_angle_5_idx_speed_12) <= r_max);
tor1_mode_freq_angle_5_idx_speed_12 = abs(tor1_mode_angle_5_idx_speed_12);
tor1_mode_damp_angle_5_idx_speed_12 = -real(tor1_mode_angle_5_idx_speed_12)./abs(tor1_mode_angle_5_idx_speed_12);
tor1_mode_damp_freq_angle_5_idx_speed_12 = [tor1_mode_damp_angle_5_idx_speed_12(1:end), tor1_mode_freq_angle_5_idx_speed_12(1:end)];

mu_tor1_mode_damp_freq_angle_5_idx_speed_12 = mean(tor1_mode_damp_freq_angle_5_idx_speed_12, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_5_idx_speed_12 = cov(tor1_mode_damp_freq_angle_5_idx_speed_12);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_5_idx_speed_12 = mvnpdf(tor1_mode_damp_freq_angle_5_idx_speed_12, mu_tor1_mode_damp_freq_angle_5_idx_speed_12, Sigma_tor1_mode_damp_freq_angle_5_idx_speed_12);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_5_idx_speed_12(:,1), tor1_mode_damp_freq_angle_5_idx_speed_12(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_5_idx_speed_12);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_5_idx_speed_12 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_5_idx_speed_12;
    plot(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_12(:,1), ellipse_tor1_mode_damp_freq_angle_5_idx_speed_12(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_5_idx_speed_12_lb = min(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_12(:, 2));
tor1_mode_freq_angle_5_idx_speed_12_ub = max(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_12(:, 2));
tor1_mode_damp_angle_5_idx_speed_12_lb = min(tor1_mode_damp_freq_angle_5_idx_speed_12(:, 1));
tor1_mode_damp_angle_5_idx_speed_12_ub = max(tor1_mode_damp_freq_angle_5_idx_speed_12(:, 1));

%% WT exp.angle 1.4, speed 27.5m/s, damp TOR1
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_18 = modalData.stabPls{5, 18};
figure; plot(modal_data_idx_angle_5_idx_speed_18, 'k.')

r_min = 12;
r_max = 13;

tor1_mode_angle_5_idx_speed_18 = modal_data_idx_angle_5_idx_speed_18(abs(modal_data_idx_angle_5_idx_speed_18) >= r_min & abs(modal_data_idx_angle_5_idx_speed_18) <= r_max);
tor1_mode_freq_angle_5_idx_speed_18 = abs(tor1_mode_angle_5_idx_speed_18);
tor1_mode_damp_angle_5_idx_speed_18 = -real(tor1_mode_angle_5_idx_speed_18)./abs(tor1_mode_angle_5_idx_speed_18);
tor1_mode_damp_freq_angle_5_idx_speed_18 = [tor1_mode_damp_angle_5_idx_speed_18(1:end), tor1_mode_freq_angle_5_idx_speed_18(1:end)];

mu_tor1_mode_damp_freq_angle_5_idx_speed_18 = mean(tor1_mode_damp_freq_angle_5_idx_speed_18, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_5_idx_speed_18 = cov(tor1_mode_damp_freq_angle_5_idx_speed_18);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_5_idx_speed_18 = mvnpdf(tor1_mode_damp_freq_angle_5_idx_speed_18, mu_tor1_mode_damp_freq_angle_5_idx_speed_18, Sigma_tor1_mode_damp_freq_angle_5_idx_speed_18);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_5_idx_speed_18(:,1), tor1_mode_damp_freq_angle_5_idx_speed_18(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_5_idx_speed_18);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_5_idx_speed_18 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_5_idx_speed_18;
    plot(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_18(:,1), ellipse_tor1_mode_damp_freq_angle_5_idx_speed_18(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_5_idx_speed_18_lb = min(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_18(:, 2));
tor1_mode_freq_angle_5_idx_speed_18_ub = max(ellipse_tor1_mode_damp_freq_angle_5_idx_speed_18(:, 2));
tor1_mode_damp_angle_5_idx_speed_18_lb = min(tor1_mode_damp_freq_angle_5_idx_speed_18(:, 1));
tor1_mode_damp_angle_5_idx_speed_18_lb_ = sort(tor1_mode_damp_freq_angle_5_idx_speed_18(:, 1), 'descend');
tor1_mode_damp_angle_5_idx_speed_18_ub = tor1_mode_damp_angle_5_idx_speed_18_lb_(6);

%% bringing together freq UB and LB
tor1_mode_angle_5_freq_lb = [tor1_mode_freq_angle_5_idx_speed_1_lb, tor1_mode_freq_angle_5_idx_speed_5_lb, tor1_mode_freq_angle_5_idx_speed_6_lb, tor1_mode_freq_angle_5_idx_speed_7_lb, tor1_mode_freq_angle_5_idx_speed_8_lb, tor1_mode_freq_angle_5_idx_speed_9_lb, tor1_mode_freq_angle_5_idx_speed_10_lb, tor1_mode_freq_angle_5_idx_speed_11_lb, tor1_mode_freq_angle_5_idx_speed_12_lb, tor1_mode_freq_angle_5_idx_speed_18_lb];
tor1_mode_angle_5_freq_ub = [tor1_mode_freq_angle_5_idx_speed_1_ub, tor1_mode_freq_angle_5_idx_speed_5_ub, tor1_mode_freq_angle_5_idx_speed_6_ub, tor1_mode_freq_angle_5_idx_speed_7_ub, tor1_mode_freq_angle_5_idx_speed_8_ub, tor1_mode_freq_angle_5_idx_speed_9_ub, tor1_mode_freq_angle_5_idx_speed_10_ub, tor1_mode_freq_angle_5_idx_speed_11_ub, tor1_mode_freq_angle_5_idx_speed_12_ub, tor1_mode_freq_angle_5_idx_speed_18_ub];
oop2_mode_angle_5_freq_lb = [oop2_mode_freq_angle_5_idx_speed_1_lb, oop2_mode_freq_angle_5_idx_speed_5_lb, oop2_mode_freq_angle_5_idx_speed_6_lb, oop2_mode_freq_angle_5_idx_speed_7_lb, oop2_mode_freq_angle_5_idx_speed_8_lb, oop2_mode_freq_angle_5_idx_speed_9_lb, oop2_mode_freq_angle_5_idx_speed_10_lb];
oop2_mode_angle_5_freq_ub = [oop2_mode_freq_angle_5_idx_speed_1_ub, oop2_mode_freq_angle_5_idx_speed_5_ub, oop2_mode_freq_angle_5_idx_speed_6_ub, oop2_mode_freq_angle_5_idx_speed_7_ub, oop2_mode_freq_angle_5_idx_speed_8_ub, oop2_mode_freq_angle_5_idx_speed_9_ub, oop2_mode_freq_angle_5_idx_speed_10_ub];

%% bringing together damp UB and LB
tor1_mode_angle_5_damp_lb = [tor1_mode_damp_angle_5_idx_speed_1_lb, tor1_mode_damp_angle_5_idx_speed_5_lb, tor1_mode_damp_angle_5_idx_speed_6_lb, tor1_mode_damp_angle_5_idx_speed_7_lb, tor1_mode_damp_angle_5_idx_speed_8_lb, tor1_mode_damp_angle_5_idx_speed_9_lb, tor1_mode_damp_angle_5_idx_speed_10_lb, tor1_mode_damp_angle_5_idx_speed_11_lb, tor1_mode_damp_angle_5_idx_speed_12_lb, tor1_mode_damp_angle_5_idx_speed_18_lb];
tor1_mode_angle_5_damp_ub = [tor1_mode_damp_angle_5_idx_speed_1_ub, tor1_mode_damp_angle_5_idx_speed_5_ub, tor1_mode_damp_angle_5_idx_speed_6_ub, tor1_mode_damp_angle_5_idx_speed_7_ub, tor1_mode_damp_angle_5_idx_speed_8_ub, tor1_mode_damp_angle_5_idx_speed_9_ub, tor1_mode_damp_angle_5_idx_speed_10_ub, tor1_mode_damp_angle_5_idx_speed_11_ub, tor1_mode_damp_angle_5_idx_speed_12_ub, tor1_mode_damp_angle_5_idx_speed_18_ub];
oop2_mode_angle_5_damp_lb = [oop2_mode_damp_angle_5_idx_speed_1_lb, oop2_mode_damp_angle_5_idx_speed_5_lb, oop2_mode_damp_angle_5_idx_speed_6_lb, oop2_mode_damp_angle_5_idx_speed_7_lb, oop2_mode_damp_angle_5_idx_speed_8_lb, oop2_mode_damp_angle_5_idx_speed_9_lb, oop2_mode_damp_angle_5_idx_speed_10_lb];
oop2_mode_angle_5_damp_ub = [oop2_mode_damp_angle_5_idx_speed_1_ub, oop2_mode_damp_angle_5_idx_speed_5_ub, oop2_mode_damp_angle_5_idx_speed_6_ub, oop2_mode_damp_angle_5_idx_speed_7_ub, oop2_mode_damp_angle_5_idx_speed_8_ub, oop2_mode_damp_angle_5_idx_speed_9_ub, oop2_mode_damp_angle_5_idx_speed_10_ub];

%% saving freq UB and LB
save('OMA_freq_lower_and_upper_bounds.mat', 'tor1_mode_angle_5_freq_lb', 'tor1_mode_angle_5_freq_ub', 'oop2_mode_angle_5_freq_lb', 'oop2_mode_angle_5_freq_ub')

%% saving damp UB and LB
save('OMA_damp_lower_and_upper_bounds.mat', 'tor1_mode_angle_5_damp_lb', 'tor1_mode_angle_5_damp_ub', 'oop2_mode_angle_5_damp_lb', 'oop2_mode_angle_5_damp_ub');

%% preparation for error bars for strains (ang_expr = 1.4)
beta_y_g1_1_bending_moment_at_0 = -2.07;
beta_y_g1_1_at_0 = beta_y_g1_1_bending_moment_at_0/2.268;
addpath('WTTests'); 
ang_expr = 1.4; %note angles to access experimental data in degrees

%function to retrieve processes experimental data...
[ang_true, exp_statResp, exp_Uf, exp_beta_yf, exp_beta_xf] = expr_statStab(ang_expr);
beta_y_wt_experimental_data_set_angle_5 = exp_statResp.beta_y;

%% preparation for error bars for strains (ang_expr = 0.3)
ang_expr = 0.3; %note angles to access experimental data in degrees

%function to retrieve processes experimental data...
[ang_true, exp_statResp, exp_Uf, exp_beta_yf, exp_beta_xf] = expr_statStab(ang_expr);
beta_y_wt_experimental_data_set_angle_2 = exp_statResp.beta_y;

%% preparation for error bars for strains (ang_expr = 0.6)
ang_expr = 0.6; %note angles to access experimental data in degrees

%function to retrieve processes experimental data...
[ang_true, exp_statResp, exp_Uf, exp_beta_yf, exp_beta_xf] = expr_statStab(ang_expr);
beta_y_wt_experimental_data_set_angle_3 = exp_statResp.beta_y;

%% preparation for error bars for strains (ang_expr = 1.1)
ang_expr = 1.1; %note angles to access experimental data in degrees

%function to retrieve processes experimental data...
[ang_true, exp_statResp, exp_Uf, exp_beta_yf, exp_beta_xf] = expr_statStab(ang_expr);
beta_y_wt_experimental_data_set_angle_4 = exp_statResp.beta_y;

%% calculate the error bars for strains
beta_y_wt_experimental_data_set_angle_5_unscaled = beta_y_wt_experimental_data_set_angle_5+beta_y_g1_1_at_0;
beta_y_wt_experimental_data_set_angle_5_lb = beta_y_wt_experimental_data_set_angle_5-0.01*abs(beta_y_wt_experimental_data_set_angle_5_unscaled)-0.01*abs(beta_y_g1_1_at_0);
beta_y_wt_experimental_data_set_angle_5_ub = beta_y_wt_experimental_data_set_angle_5+0.01*abs(beta_y_wt_experimental_data_set_angle_5_unscaled)+0.01*abs(beta_y_g1_1_at_0);

beta_y_wt_experimental_data_set_angle_4_unscaled = beta_y_wt_experimental_data_set_angle_4+beta_y_g1_1_at_0;
beta_y_wt_experimental_data_set_angle_4_lb = beta_y_wt_experimental_data_set_angle_4-0.01*abs(beta_y_wt_experimental_data_set_angle_4_unscaled)-0.01*abs(beta_y_g1_1_at_0);
beta_y_wt_experimental_data_set_angle_4_ub = beta_y_wt_experimental_data_set_angle_4+0.01*abs(beta_y_wt_experimental_data_set_angle_4_unscaled)+0.01*abs(beta_y_g1_1_at_0);

beta_y_wt_experimental_data_set_angle_3_unscaled = beta_y_wt_experimental_data_set_angle_3+beta_y_g1_1_at_0;
beta_y_wt_experimental_data_set_angle_3_lb = beta_y_wt_experimental_data_set_angle_3-0.01*abs(beta_y_wt_experimental_data_set_angle_3_unscaled)-0.01*abs(beta_y_g1_1_at_0);
beta_y_wt_experimental_data_set_angle_3_ub = beta_y_wt_experimental_data_set_angle_3+0.01*abs(beta_y_wt_experimental_data_set_angle_3_unscaled)+0.01*abs(beta_y_g1_1_at_0);

beta_y_wt_experimental_data_set_angle_2_unscaled = beta_y_wt_experimental_data_set_angle_2+beta_y_g1_1_at_0;
beta_y_wt_experimental_data_set_angle_2_lb = beta_y_wt_experimental_data_set_angle_2-0.01*abs(beta_y_wt_experimental_data_set_angle_2_unscaled)-0.01*abs(beta_y_g1_1_at_0);
beta_y_wt_experimental_data_set_angle_2_ub = beta_y_wt_experimental_data_set_angle_2+0.01*abs(beta_y_wt_experimental_data_set_angle_2_unscaled)+0.01*abs(beta_y_g1_1_at_0);

save('WT_strains_lower_and_upper_bounds.mat', 'beta_y_wt_experimental_data_set_angle_2_lb', 'beta_y_wt_experimental_data_set_angle_2_ub', 'beta_y_wt_experimental_data_set_angle_3_lb', 'beta_y_wt_experimental_data_set_angle_3_ub', 'beta_y_wt_experimental_data_set_angle_4_lb', 'beta_y_wt_experimental_data_set_angle_4_ub', 'beta_y_wt_experimental_data_set_angle_5_lb', 'beta_y_wt_experimental_data_set_angle_5_ub');

%%
% load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
% modal_data_idx_angle_5_idx_speed_6 = modalData.stabPls{5, 6};
% figure; plot(modal_data_idx_angle_5_idx_speed_6, 'k.')
% 
% r_min = 1;
% r_max = 3;
% 
% oop1_mode = modal_data_idx_angle_5_idx_speed_6(abs(modal_data_idx_angle_5_idx_speed_6) >= r_min & abs(modal_data_idx_angle_5_idx_speed_6) <= r_max);
% oop1_mode_freq = abs(oop1_mode);
% oop1_mode_damp = -real(oop1_mode)./abs(oop1_mode);
% oop1_mode_damp_freq = [oop1_mode_damp, oop1_mode_freq];
% 
% mu_oop1_mode_damp_freq = mean(oop1_mode_damp_freq, 1);        % 1×2 mean vector
% Sigma_oop1_mode_damp_freq = cov(oop1_mode_damp_freq);         % 2×2 covariance matrix
% 
% pdf_vals_oop1_mode_damp_freq = mvnpdf(oop1_mode_damp_freq, mu_oop1_mode_damp_freq, Sigma_oop1_mode_damp_freq);
% 
% figure; hold on; axis equal
% scatter(oop1_mode_damp_freq(:,1), oop1_mode_damp_freq(:,2), 10, 'filled')
% 
% theta = linspace(0,2*pi,200);
% unit_circle = [cos(theta); sin(theta)];
% 
% [V,D] = eig(Sigma_oop1_mode_damp_freq);
% 
% for k = 1:3   % 1σ, 2σ, 3σ
%     ellipse = (k * V * sqrt(D) * unit_circle)' + mu_oop1_mode_damp_freq;
%     plot(ellipse(:,1), ellipse(:,2), 'LineWidth', 2)
% end
% 
% legend('Data','1σ','2σ','3σ')
% 
% 
% K = 2;   % number of Gaussians
% gm = fitgmdist(oop1_mode_damp_freq, K);
% idx = cluster(gm, oop1_mode_damp_freq);
% 
% 
% p = pdf(gm, oop1_mode_damp_freq);            % total mixture density
% p_k = posterior(gm, oop1_mode_damp_freq);    % N×K soft assignments
% 
% 
% figure; hold on; axis equal
% scatter(oop1_mode_damp_freq(:,1), oop1_mode_damp_freq(:,2), 8, 'k', 'filled')
% 
% theta = linspace(0,2*pi,200);
% circle = [cos(theta); sin(theta)];
% 
% for k = 1:K
%     [V,D] = eig(gm.Sigma(:,:,k));
%     ellipse = (V * sqrt(D) * circle)' + gm.mu(k,:);
%     plot(ellipse(:,1), ellipse(:,2), 'LineWidth', 2)
% end
% 
