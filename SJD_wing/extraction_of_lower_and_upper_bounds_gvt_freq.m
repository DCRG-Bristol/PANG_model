%% GVT angle 0, freq OOP1
load('gvt_poles.mat');
modal_data_angle_0 = pl_all{1, 6};
figure; plot(modal_data_angle_0, 'k.')

r_min = 1;
r_max = 3;

oop1_mode_angle_0 = modal_data_angle_0(abs(modal_data_angle_0) >= r_min & abs(modal_data_angle_0) <= r_max);
oop1_mode_freq_angle_0 = abs(oop1_mode_angle_0);
oop1_mode_damp_angle_0 = -real(oop1_mode_angle_0)./abs(oop1_mode_angle_0);
oop1_mode_damp_freq_angle_0 = [oop1_mode_damp_angle_0(3:end), oop1_mode_freq_angle_0(3:end)];

mu_oop1_mode_damp_freq_angle_0 = mean(oop1_mode_damp_freq_angle_0, 1);        % 1×2 mean vector
Sigma_oop1_mode_damp_freq_angle_0 = cov(oop1_mode_damp_freq_angle_0);         % 2×2 covariance matrix

pdf_vals_oop1_mode_damp_freq_angle_0 = mvnpdf(oop1_mode_damp_freq_angle_0, mu_oop1_mode_damp_freq_angle_0, Sigma_oop1_mode_damp_freq_angle_0);

figure; hold on; axis equal
scatter(oop1_mode_damp_freq_angle_0(:,1), oop1_mode_damp_freq_angle_0(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop1_mode_damp_freq_angle_0);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop1_mode_damp_freq_angle_0 = (k * V * sqrt(D) * unit_circle)' + mu_oop1_mode_damp_freq_angle_0;
    plot(ellipse_oop1_mode_damp_freq_angle_0(:,1), ellipse_oop1_mode_damp_freq_angle_0(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop1_mode_freq_angle_0_lb = min(ellipse_oop1_mode_damp_freq_angle_0(:, 2));
oop1_mode_freq_angle_0_ub = max(ellipse_oop1_mode_damp_freq_angle_0(:, 2));

%% GVT angle 0, freq IP1
load('gvt_poles.mat');
modal_data_angle_0 = pl_all{1, 6};
% figure; plot(modal_data_angle_0, 'k.')

r_min = 8;
r_max = 9;

ip1_mode_angle_0 = modal_data_angle_0(abs(modal_data_angle_0) >= r_min & abs(modal_data_angle_0) <= r_max);
ip1_mode_freq_angle_0 = abs(ip1_mode_angle_0);
ip1_mode_damp_angle_0 = -real(ip1_mode_angle_0)./abs(ip1_mode_angle_0);
ip1_mode_damp_freq_angle_0 = [ip1_mode_damp_angle_0(3:end), ip1_mode_freq_angle_0(3:end)];

mu_ip1_mode_damp_freq_angle_0 = mean(ip1_mode_damp_freq_angle_0, 1);        % 1×2 mean vector
Sigma_ip1_mode_damp_freq_angle_0 = cov(ip1_mode_damp_freq_angle_0);         % 2×2 covariance matrix

pdf_vals_ip1_mode_damp_freq_angle_0 = mvnpdf(ip1_mode_damp_freq_angle_0, mu_ip1_mode_damp_freq_angle_0, Sigma_ip1_mode_damp_freq_angle_0);

figure; hold on; axis equal
scatter(ip1_mode_damp_freq_angle_0(:,1), ip1_mode_damp_freq_angle_0(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_ip1_mode_damp_freq_angle_0);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_ip1_mode_damp_freq_angle_0 = (k * V * sqrt(D) * unit_circle)' + mu_ip1_mode_damp_freq_angle_0;
    plot(ellipse_ip1_mode_damp_freq_angle_0(:,1), ellipse_ip1_mode_damp_freq_angle_0(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

ip1_mode_freq_angle_0_lb = min(ellipse_ip1_mode_damp_freq_angle_0(:, 2));
ip1_mode_freq_angle_0_ub = max(ellipse_ip1_mode_damp_freq_angle_0(:, 2));

%% GVT angle 0, freq OOP2
load('gvt_poles.mat');
modal_data_angle_0 = pl_all{1, 6};
figure; plot(modal_data_angle_0, 'k.')

r_min =11.1;
r_max = 11.5;

oop2_mode_angle_0 = modal_data_angle_0(abs(modal_data_angle_0) >= r_min & abs(modal_data_angle_0) <= r_max);
oop2_mode_freq_angle_0 = abs(oop2_mode_angle_0);
oop2_mode_damp_angle_0 = -real(oop2_mode_angle_0)./abs(oop2_mode_angle_0);
oop2_mode_damp_freq_angle_0 = [oop2_mode_damp_angle_0(1:end), oop2_mode_freq_angle_0(1:end)];

mu_oop2_mode_damp_freq_angle_0 = mean(oop2_mode_damp_freq_angle_0, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_0 = cov(oop2_mode_damp_freq_angle_0);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_0 = mvnpdf(oop2_mode_damp_freq_angle_0, mu_oop2_mode_damp_freq_angle_0, Sigma_oop2_mode_damp_freq_angle_0);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_0(:,1), oop2_mode_damp_freq_angle_0(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_0);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_0 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_0;
    plot(ellipse_oop2_mode_damp_freq_angle_0(:,1), ellipse_oop2_mode_damp_freq_angle_0(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_0_lb = min(ellipse_oop2_mode_damp_freq_angle_0(:, 2));
oop2_mode_freq_angle_0_ub = max(ellipse_oop2_mode_damp_freq_angle_0(:, 2));

%% GVT angle 0, freq TOR1
load('gvt_poles.mat');
modal_data_angle_0 = pl_all{1, 6};
% figure; plot(modal_data_angle_0, 'k.')

r_min = 16;
r_max = 17;

tor1_mode_angle_0 = modal_data_angle_0(abs(modal_data_angle_0) >= r_min & abs(modal_data_angle_0) <= r_max);
tor1_mode_freq_angle_0 = abs(tor1_mode_angle_0);
tor1_mode_damp_angle_0 = -real(tor1_mode_angle_0)./abs(tor1_mode_angle_0);
tor1_mode_damp_freq_angle_0 = [tor1_mode_damp_angle_0(1:end), tor1_mode_freq_angle_0(1:end)];

mu_tor1_mode_damp_freq_angle_0 = mean(tor1_mode_damp_freq_angle_0, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_0 = cov(tor1_mode_damp_freq_angle_0);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_0 = mvnpdf(tor1_mode_damp_freq_angle_0, mu_tor1_mode_damp_freq_angle_0, Sigma_tor1_mode_damp_freq_angle_0);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_0(:,1), tor1_mode_damp_freq_angle_0(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_0);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_0 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_0;
    plot(ellipse_tor1_mode_damp_freq_angle_0(:,1), ellipse_tor1_mode_damp_freq_angle_0(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_0_lb = min(ellipse_tor1_mode_damp_freq_angle_0(:, 2));
tor1_mode_freq_angle_0_ub = max(ellipse_tor1_mode_damp_freq_angle_0(:, 2));


%% GVT angle 10, freq OOP1
load('gvt_poles.mat');
modal_data_angle_10 = pl_all{1, 5};
figure; plot(modal_data_angle_10, 'k.')

r_min = 1;
r_max = 3;

oop1_mode_angle_10 = modal_data_angle_10(abs(modal_data_angle_10) >= r_min & abs(modal_data_angle_10) <= r_max);
oop1_mode_freq_angle_10 = abs(oop1_mode_angle_10);
oop1_mode_damp_angle_10 = -real(oop1_mode_angle_10)./abs(oop1_mode_angle_10);
oop1_mode_damp_freq_angle_10 = [oop1_mode_damp_angle_10(2:end), oop1_mode_freq_angle_10(2:end)];

mu_oop1_mode_damp_freq_angle_10 = mean(oop1_mode_damp_freq_angle_10, 1);        % 1×2 mean vector
Sigma_oop1_mode_damp_freq_angle_10 = cov(oop1_mode_damp_freq_angle_10);         % 2×2 covariance matrix

pdf_vals_oop1_mode_damp_freq_angle_10 = mvnpdf(oop1_mode_damp_freq_angle_10, mu_oop1_mode_damp_freq_angle_10, Sigma_oop1_mode_damp_freq_angle_10);

figure; hold on; axis equal
scatter(oop1_mode_damp_freq_angle_10(:,1), oop1_mode_damp_freq_angle_10(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop1_mode_damp_freq_angle_10);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop1_mode_damp_freq_angle_10 = (k * V * sqrt(D) * unit_circle)' + mu_oop1_mode_damp_freq_angle_10;
    plot(ellipse_oop1_mode_damp_freq_angle_10(:,1), ellipse_oop1_mode_damp_freq_angle_10(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop1_mode_freq_angle_10_lb = min(ellipse_oop1_mode_damp_freq_angle_10(:, 2));
oop1_mode_freq_angle_10_ub = max(ellipse_oop1_mode_damp_freq_angle_10(:, 2));

%% GVT angle 10, freq IP1
load('gvt_poles.mat');
modal_data_angle_10 = pl_all{1, 5};
% figure; plot(modal_data_angle_10, 'k.')

r_min = 8;
r_max = 9;

ip1_mode_angle_10 = modal_data_angle_10(abs(modal_data_angle_10) >= r_min & abs(modal_data_angle_10) <= r_max);
ip1_mode_freq_angle_10 = abs(ip1_mode_angle_10);
ip1_mode_damp_angle_10 = -real(ip1_mode_angle_10)./abs(ip1_mode_angle_10);
ip1_mode_damp_freq_angle_10 = [ip1_mode_damp_angle_10(1:end), ip1_mode_freq_angle_10(1:end)];

mu_ip1_mode_damp_freq_angle_10 = mean(ip1_mode_damp_freq_angle_10, 1);        % 1×2 mean vector
Sigma_ip1_mode_damp_freq_angle_10 = cov(ip1_mode_damp_freq_angle_10);         % 2×2 covariance matrix

pdf_vals_ip1_mode_damp_freq_angle_10 = mvnpdf(ip1_mode_damp_freq_angle_10, mu_ip1_mode_damp_freq_angle_10, Sigma_ip1_mode_damp_freq_angle_10);

figure; hold on; axis equal
scatter(ip1_mode_damp_freq_angle_10(:,1), ip1_mode_damp_freq_angle_10(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_ip1_mode_damp_freq_angle_10);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_ip1_mode_damp_freq_angle_10 = (k * V * sqrt(D) * unit_circle)' + mu_ip1_mode_damp_freq_angle_10;
    plot(ellipse_ip1_mode_damp_freq_angle_10(:,1), ellipse_ip1_mode_damp_freq_angle_10(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

ip1_mode_freq_angle_10_lb = min(ellipse_ip1_mode_damp_freq_angle_10(:, 2));
ip1_mode_freq_angle_10_ub = max(ellipse_ip1_mode_damp_freq_angle_10(:, 2));

%% GVT angle 10, freq OOP2
load('gvt_poles.mat');
modal_data_angle_10 = pl_all{1, 5};
figure; plot(modal_data_angle_10, 'k.')

r_min =11.1;
r_max = 11.5;

oop2_mode_angle_10 = modal_data_angle_10(abs(modal_data_angle_10) >= r_min & abs(modal_data_angle_10) <= r_max);
oop2_mode_freq_angle_10 = abs(oop2_mode_angle_10);
oop2_mode_damp_angle_10 = -real(oop2_mode_angle_10)./abs(oop2_mode_angle_10);
oop2_mode_damp_freq_angle_10 = [oop2_mode_damp_angle_10(1:end), oop2_mode_freq_angle_10(1:end)];

mu_oop2_mode_damp_freq_angle_10 = mean(oop2_mode_damp_freq_angle_10, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_10 = cov(oop2_mode_damp_freq_angle_10);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_10 = mvnpdf(oop2_mode_damp_freq_angle_10, mu_oop2_mode_damp_freq_angle_10, Sigma_oop2_mode_damp_freq_angle_10);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_10(:,1), oop2_mode_damp_freq_angle_10(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_10);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_10 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_10;
    plot(ellipse_oop2_mode_damp_freq_angle_10(:,1), ellipse_oop2_mode_damp_freq_angle_10(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_10_lb = min(ellipse_oop2_mode_damp_freq_angle_10(:, 2));
oop2_mode_freq_angle_10_ub = max(ellipse_oop2_mode_damp_freq_angle_10(:, 2));

%% GVT angle 10, freq TOR1
load('gvt_poles.mat');
modal_data_angle_10 = pl_all{1, 5};
% figure; plot(modal_data_angle_10, 'k.')

r_min = 16.5;
r_max = 17;

tor1_mode_angle_10 = modal_data_angle_10(abs(modal_data_angle_10) >= r_min & abs(modal_data_angle_10) <= r_max);
tor1_mode_freq_angle_10 = abs(tor1_mode_angle_10);
tor1_mode_damp_angle_10 = -real(tor1_mode_angle_10)./abs(tor1_mode_angle_10);
tor1_mode_damp_freq_angle_10 = [tor1_mode_damp_angle_10(1:end), tor1_mode_freq_angle_10(1:end)];

mu_tor1_mode_damp_freq_angle_10 = mean(tor1_mode_damp_freq_angle_10, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_10 = cov(tor1_mode_damp_freq_angle_10);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_10 = mvnpdf(tor1_mode_damp_freq_angle_10, mu_tor1_mode_damp_freq_angle_10, Sigma_tor1_mode_damp_freq_angle_10);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_10(:,1), tor1_mode_damp_freq_angle_10(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_10);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_10 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_10;
    plot(ellipse_tor1_mode_damp_freq_angle_10(:,1), ellipse_tor1_mode_damp_freq_angle_10(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_10_lb = min(ellipse_tor1_mode_damp_freq_angle_10(:, 2));
tor1_mode_freq_angle_10_ub = max(ellipse_tor1_mode_damp_freq_angle_10(:, 2));

%% GVT angle 20, freq OOP1
load('gvt_poles.mat');
modal_data_angle_20 = pl_all{1, 4};
figure; plot(modal_data_angle_20, 'k.')

r_min = 1;
r_max = 3;

oop1_mode_angle_20 = modal_data_angle_20(abs(modal_data_angle_20) >= r_min & abs(modal_data_angle_20) <= r_max);
oop1_mode_freq_angle_20 = abs(oop1_mode_angle_20);
oop1_mode_damp_angle_20 = -real(oop1_mode_angle_20)./abs(oop1_mode_angle_20);
oop1_mode_damp_freq_angle_20 = [oop1_mode_damp_angle_20(1:end), oop1_mode_freq_angle_20(1:end)];

mu_oop1_mode_damp_freq_angle_20 = mean(oop1_mode_damp_freq_angle_20, 1);        % 1×2 mean vector
Sigma_oop1_mode_damp_freq_angle_20 = cov(oop1_mode_damp_freq_angle_20);         % 2×2 covariance matrix

pdf_vals_oop1_mode_damp_freq_angle_20 = mvnpdf(oop1_mode_damp_freq_angle_20, mu_oop1_mode_damp_freq_angle_20, Sigma_oop1_mode_damp_freq_angle_20);

figure; hold on; axis equal
scatter(oop1_mode_damp_freq_angle_20(:,1), oop1_mode_damp_freq_angle_20(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop1_mode_damp_freq_angle_20);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop1_mode_damp_freq_angle_20 = (k * V * sqrt(D) * unit_circle)' + mu_oop1_mode_damp_freq_angle_20;
    plot(ellipse_oop1_mode_damp_freq_angle_20(:,1), ellipse_oop1_mode_damp_freq_angle_20(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop1_mode_freq_angle_20_lb = min(ellipse_oop1_mode_damp_freq_angle_20(:, 2));
oop1_mode_freq_angle_20_ub = max(ellipse_oop1_mode_damp_freq_angle_20(:, 2));

%% GVT angle 20, freq IP1
load('gvt_poles.mat');
modal_data_angle_20 = pl_all{1, 4};
% figure; plot(modal_data_angle_20, 'k.')

r_min = 8.5;
r_max = 9;

ip1_mode_angle_20 = modal_data_angle_20(abs(modal_data_angle_20) >= r_min & abs(modal_data_angle_20) <= r_max);
ip1_mode_freq_angle_20 = abs(ip1_mode_angle_20);
ip1_mode_damp_angle_20 = -real(ip1_mode_angle_20)./abs(ip1_mode_angle_20);
ip1_mode_damp_freq_angle_20 = [ip1_mode_damp_angle_20(1:end), ip1_mode_freq_angle_20(1:end)];

mu_ip1_mode_damp_freq_angle_20 = mean(ip1_mode_damp_freq_angle_20, 1);        % 1×2 mean vector
Sigma_ip1_mode_damp_freq_angle_20 = cov(ip1_mode_damp_freq_angle_20);         % 2×2 covariance matrix

pdf_vals_ip1_mode_damp_freq_angle_20 = mvnpdf(ip1_mode_damp_freq_angle_20, mu_ip1_mode_damp_freq_angle_20, Sigma_ip1_mode_damp_freq_angle_20);

figure; hold on; axis equal
scatter(ip1_mode_damp_freq_angle_20(:,1), ip1_mode_damp_freq_angle_20(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_ip1_mode_damp_freq_angle_20);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_ip1_mode_damp_freq_angle_20 = (k * V * sqrt(D) * unit_circle)' + mu_ip1_mode_damp_freq_angle_20;
    plot(ellipse_ip1_mode_damp_freq_angle_20(:,1), ellipse_ip1_mode_damp_freq_angle_20(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

ip1_mode_freq_angle_20_lb = min(ellipse_ip1_mode_damp_freq_angle_20(:, 2));
ip1_mode_freq_angle_20_ub = max(ellipse_ip1_mode_damp_freq_angle_20(:, 2));

%% GVT angle 20, freq OOP2
load('gvt_poles.mat');
modal_data_angle_20 = pl_all{1, 4};
figure; plot(modal_data_angle_20, 'k.')

r_min =11;
r_max = 11.5;

oop2_mode_angle_20 = modal_data_angle_20(abs(modal_data_angle_20) >= r_min & abs(modal_data_angle_20) <= r_max);
oop2_mode_freq_angle_20 = abs(oop2_mode_angle_20);
oop2_mode_damp_angle_20 = -real(oop2_mode_angle_20)./abs(oop2_mode_angle_20);
oop2_mode_damp_freq_angle_20 = [oop2_mode_damp_angle_20(1:end), oop2_mode_freq_angle_20(1:end)];

mu_oop2_mode_damp_freq_angle_20 = mean(oop2_mode_damp_freq_angle_20, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_20 = cov(oop2_mode_damp_freq_angle_20);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_20 = mvnpdf(oop2_mode_damp_freq_angle_20, mu_oop2_mode_damp_freq_angle_20, Sigma_oop2_mode_damp_freq_angle_20);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_20(:,1), oop2_mode_damp_freq_angle_20(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_20);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_20 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_20;
    plot(ellipse_oop2_mode_damp_freq_angle_20(:,1), ellipse_oop2_mode_damp_freq_angle_20(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_20_lb = min(ellipse_oop2_mode_damp_freq_angle_20(:, 2));
oop2_mode_freq_angle_20_ub = max(ellipse_oop2_mode_damp_freq_angle_20(:, 2));

%% GVT angle 20, freq TOR1
load('gvt_poles.mat');
modal_data_angle_20 = pl_all{1, 4};
% figure; plot(modal_data_angle_20, 'k.')

r_min = 16;
r_max = 17;

tor1_mode_angle_20 = modal_data_angle_20(abs(modal_data_angle_20) >= r_min & abs(modal_data_angle_20) <= r_max);
tor1_mode_freq_angle_20 = abs(tor1_mode_angle_20);
tor1_mode_damp_angle_20 = -real(tor1_mode_angle_20)./abs(tor1_mode_angle_20);
tor1_mode_damp_freq_angle_20 = [tor1_mode_damp_angle_20(1:end), tor1_mode_freq_angle_20(1:end)];

mu_tor1_mode_damp_freq_angle_20 = mean(tor1_mode_damp_freq_angle_20, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_20 = cov(tor1_mode_damp_freq_angle_20);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_20 = mvnpdf(tor1_mode_damp_freq_angle_20, mu_tor1_mode_damp_freq_angle_20, Sigma_tor1_mode_damp_freq_angle_20);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_20(:,1), tor1_mode_damp_freq_angle_20(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_20);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_20 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_20;
    plot(ellipse_tor1_mode_damp_freq_angle_20(:,1), ellipse_tor1_mode_damp_freq_angle_20(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_20_lb = min(ellipse_tor1_mode_damp_freq_angle_20(:, 2));
tor1_mode_freq_angle_20_ub = max(ellipse_tor1_mode_damp_freq_angle_20(:, 2));

%% GVT angle 30, freq OOP1
load('gvt_poles.mat');
modal_data_angle_30 = pl_all{1, 3};
figure; plot(modal_data_angle_30, 'k.')

r_min = 1;
r_max = 3;

oop1_mode_angle_30 = modal_data_angle_30(abs(modal_data_angle_30) >= r_min & abs(modal_data_angle_30) <= r_max);
oop1_mode_freq_angle_30 = abs(oop1_mode_angle_30);
oop1_mode_damp_angle_30 = -real(oop1_mode_angle_30)./abs(oop1_mode_angle_30);
oop1_mode_damp_freq_angle_30 = [oop1_mode_damp_angle_30(1:end), oop1_mode_freq_angle_30(1:end)];

mu_oop1_mode_damp_freq_angle_30 = mean(oop1_mode_damp_freq_angle_30, 1);        % 1×2 mean vector
Sigma_oop1_mode_damp_freq_angle_30 = cov(oop1_mode_damp_freq_angle_30);         % 2×2 covariance matrix

pdf_vals_oop1_mode_damp_freq_angle_30 = mvnpdf(oop1_mode_damp_freq_angle_30, mu_oop1_mode_damp_freq_angle_30, Sigma_oop1_mode_damp_freq_angle_30);

figure; hold on; axis equal
scatter(oop1_mode_damp_freq_angle_30(:,1), oop1_mode_damp_freq_angle_30(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop1_mode_damp_freq_angle_30);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop1_mode_damp_freq_angle_30 = (k * V * sqrt(D) * unit_circle)' + mu_oop1_mode_damp_freq_angle_30;
    plot(ellipse_oop1_mode_damp_freq_angle_30(:,1), ellipse_oop1_mode_damp_freq_angle_30(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop1_mode_freq_angle_30_lb = min(ellipse_oop1_mode_damp_freq_angle_30(:, 2));
oop1_mode_freq_angle_30_ub = max(ellipse_oop1_mode_damp_freq_angle_30(:, 2));

%% GVT angle 30, freq IP1
load('gvt_poles.mat');
modal_data_angle_30 = pl_all{1, 3};
% figure; plot(modal_data_angle_30, 'k.')

r_min = 8.5;
r_max = 9;

ip1_mode_angle_30 = modal_data_angle_30(abs(modal_data_angle_30) >= r_min & abs(modal_data_angle_30) <= r_max);
ip1_mode_freq_angle_30 = abs(ip1_mode_angle_30);
ip1_mode_damp_angle_30 = -real(ip1_mode_angle_30)./abs(ip1_mode_angle_30);
ip1_mode_damp_freq_angle_30 = [ip1_mode_damp_angle_30(1:end), ip1_mode_freq_angle_30(1:end)];

mu_ip1_mode_damp_freq_angle_30 = mean(ip1_mode_damp_freq_angle_30, 1);        % 1×2 mean vector
Sigma_ip1_mode_damp_freq_angle_30 = cov(ip1_mode_damp_freq_angle_30);         % 2×2 covariance matrix

pdf_vals_ip1_mode_damp_freq_angle_30 = mvnpdf(ip1_mode_damp_freq_angle_30, mu_ip1_mode_damp_freq_angle_30, Sigma_ip1_mode_damp_freq_angle_30);

figure; hold on; axis equal
scatter(ip1_mode_damp_freq_angle_30(:,1), ip1_mode_damp_freq_angle_30(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_ip1_mode_damp_freq_angle_30);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_ip1_mode_damp_freq_angle_30 = (k * V * sqrt(D) * unit_circle)' + mu_ip1_mode_damp_freq_angle_30;
    plot(ellipse_ip1_mode_damp_freq_angle_30(:,1), ellipse_ip1_mode_damp_freq_angle_30(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

ip1_mode_freq_angle_30_lb = min(ellipse_ip1_mode_damp_freq_angle_30(:, 2));
ip1_mode_freq_angle_30_ub = max(ellipse_ip1_mode_damp_freq_angle_30(:, 2));

%% GVT angle 30, freq OOP2
load('gvt_poles.mat');
modal_data_angle_30 = pl_all{1, 3};
figure; plot(modal_data_angle_30, 'k.')

r_min =11;
r_max = 11.45;

oop2_mode_angle_30 = modal_data_angle_30(abs(modal_data_angle_30) >= r_min & abs(modal_data_angle_30) <= r_max);
oop2_mode_freq_angle_30 = abs(oop2_mode_angle_30);
oop2_mode_damp_angle_30 = -real(oop2_mode_angle_30)./abs(oop2_mode_angle_30);
oop2_mode_damp_freq_angle_30 = [oop2_mode_damp_angle_30(1:end), oop2_mode_freq_angle_30(1:end)];

mu_oop2_mode_damp_freq_angle_30 = mean(oop2_mode_damp_freq_angle_30, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_30 = cov(oop2_mode_damp_freq_angle_30);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_30 = mvnpdf(oop2_mode_damp_freq_angle_30, mu_oop2_mode_damp_freq_angle_30, Sigma_oop2_mode_damp_freq_angle_30);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_30(:,1), oop2_mode_damp_freq_angle_30(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_30);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_30 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_30;
    plot(ellipse_oop2_mode_damp_freq_angle_30(:,1), ellipse_oop2_mode_damp_freq_angle_30(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_30_lb = min(ellipse_oop2_mode_damp_freq_angle_30(:, 2));
oop2_mode_freq_angle_30_ub = max(ellipse_oop2_mode_damp_freq_angle_30(:, 2));

%% GVT angle 30, freq TOR1
load('gvt_poles.mat');
modal_data_angle_30 = pl_all{1, 3};
% figure; plot(modal_data_angle_30, 'k.')

r_min = 16.4;
r_max = 17;

tor1_mode_angle_30 = modal_data_angle_30(abs(modal_data_angle_30) >= r_min & abs(modal_data_angle_30) <= r_max);
tor1_mode_freq_angle_30 = abs(tor1_mode_angle_30);
tor1_mode_damp_angle_30 = -real(tor1_mode_angle_30)./abs(tor1_mode_angle_30);
tor1_mode_damp_freq_angle_30 = [tor1_mode_damp_angle_30(1:end), tor1_mode_freq_angle_30(1:end)];

mu_tor1_mode_damp_freq_angle_30 = mean(tor1_mode_damp_freq_angle_30, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_30 = cov(tor1_mode_damp_freq_angle_30);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_30 = mvnpdf(tor1_mode_damp_freq_angle_30, mu_tor1_mode_damp_freq_angle_30, Sigma_tor1_mode_damp_freq_angle_30);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_30(:,1), tor1_mode_damp_freq_angle_30(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_30);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_30 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_30;
    plot(ellipse_tor1_mode_damp_freq_angle_30(:,1), ellipse_tor1_mode_damp_freq_angle_30(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_30_lb = min(ellipse_tor1_mode_damp_freq_angle_30(:, 2));
tor1_mode_freq_angle_30_ub = max(ellipse_tor1_mode_damp_freq_angle_30(:, 2));

%% GVT angle 60, freq OOP1
load('gvt_poles.mat');
modal_data_angle_60 = pl_all{1, 2};
figure; plot(modal_data_angle_60, 'k.')

r_min = 1;
r_max = 3;

oop1_mode_angle_60 = modal_data_angle_60(abs(modal_data_angle_60) >= r_min & abs(modal_data_angle_60) <= r_max);
oop1_mode_freq_angle_60 = abs(oop1_mode_angle_60);
oop1_mode_damp_angle_60 = -real(oop1_mode_angle_60)./abs(oop1_mode_angle_60);
oop1_mode_damp_freq_angle_60 = [oop1_mode_damp_angle_60(1:end), oop1_mode_freq_angle_60(1:end)];

mu_oop1_mode_damp_freq_angle_60 = mean(oop1_mode_damp_freq_angle_60, 1);        % 1×2 mean vector
Sigma_oop1_mode_damp_freq_angle_60 = cov(oop1_mode_damp_freq_angle_60);         % 2×2 covariance matrix

pdf_vals_oop1_mode_damp_freq_angle_60 = mvnpdf(oop1_mode_damp_freq_angle_60, mu_oop1_mode_damp_freq_angle_60, Sigma_oop1_mode_damp_freq_angle_60);

figure; hold on; axis equal
scatter(oop1_mode_damp_freq_angle_60(:,1), oop1_mode_damp_freq_angle_60(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop1_mode_damp_freq_angle_60);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop1_mode_damp_freq_angle_60 = (k * V * sqrt(D) * unit_circle)' + mu_oop1_mode_damp_freq_angle_60;
    plot(ellipse_oop1_mode_damp_freq_angle_60(:,1), ellipse_oop1_mode_damp_freq_angle_60(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop1_mode_freq_angle_60_lb = min(ellipse_oop1_mode_damp_freq_angle_60(:, 2));
oop1_mode_freq_angle_60_ub = max(ellipse_oop1_mode_damp_freq_angle_60(:, 2));

%% GVT angle 60, freq IP1
load('gvt_poles.mat');
modal_data_angle_60 = pl_all{1, 2};
figure; plot(modal_data_angle_60, 'k.')

r_min = 8.5;
r_max = 9.5;

ip1_mode_angle_60 = modal_data_angle_60(abs(modal_data_angle_60) >= r_min & abs(modal_data_angle_60) <= r_max);
ip1_mode_freq_angle_60 = abs(ip1_mode_angle_60);
ip1_mode_damp_angle_60 = -real(ip1_mode_angle_60)./abs(ip1_mode_angle_60);
ip1_mode_damp_freq_angle_60 = [ip1_mode_damp_angle_60(1:end), ip1_mode_freq_angle_60(1:end)];

ip1_mode_freq_angle_60_lb = min(ip1_mode_freq_angle_60(:, 1));
ip1_mode_freq_angle_60_ub = max(ip1_mode_freq_angle_60(:, 1));

%% GVT angle 60, freq OOP2
load('gvt_poles.mat');
modal_data_angle_60 = pl_all{1, 2};
figure; plot(modal_data_angle_60, 'k.')

r_min =11;
r_max = 11.5;

oop2_mode_angle_60 = modal_data_angle_60(abs(modal_data_angle_60) >= r_min & abs(modal_data_angle_60) <= r_max);
oop2_mode_freq_angle_60 = abs(oop2_mode_angle_60);
oop2_mode_damp_angle_60 = -real(oop2_mode_angle_60)./abs(oop2_mode_angle_60);
oop2_mode_damp_freq_angle_60 = [oop2_mode_damp_angle_60(1:end), oop2_mode_freq_angle_60(1:end)];

mu_oop2_mode_damp_freq_angle_60 = mean(oop2_mode_damp_freq_angle_60, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_60 = cov(oop2_mode_damp_freq_angle_60);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_60 = mvnpdf(oop2_mode_damp_freq_angle_60, mu_oop2_mode_damp_freq_angle_60, Sigma_oop2_mode_damp_freq_angle_60);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_60(:,1), oop2_mode_damp_freq_angle_60(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_60);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_60 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_60;
    plot(ellipse_oop2_mode_damp_freq_angle_60(:,1), ellipse_oop2_mode_damp_freq_angle_60(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_60_lb = min(ellipse_oop2_mode_damp_freq_angle_60(:, 2));
oop2_mode_freq_angle_60_ub = max(ellipse_oop2_mode_damp_freq_angle_60(:, 2));

%% GVT angle 60, freq TOR1
load('gvt_poles.mat');
modal_data_angle_60 = pl_all{1, 2};
figure; plot(modal_data_angle_60, 'k.')

r_min = 15;
r_max = 17;

tor1_mode_angle_60 = modal_data_angle_60(abs(modal_data_angle_60) >= r_min & abs(modal_data_angle_60) <= r_max);
tor1_mode_freq_angle_60 = abs(tor1_mode_angle_60);
tor1_mode_damp_angle_60 = -real(tor1_mode_angle_60)./abs(tor1_mode_angle_60);
tor1_mode_damp_freq_angle_60 = [tor1_mode_damp_angle_60(1:end), tor1_mode_freq_angle_60(1:end)];

mu_tor1_mode_damp_freq_angle_60 = mean(tor1_mode_damp_freq_angle_60, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_60 = cov(tor1_mode_damp_freq_angle_60);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_60 = mvnpdf(tor1_mode_damp_freq_angle_60, mu_tor1_mode_damp_freq_angle_60, Sigma_tor1_mode_damp_freq_angle_60);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_60(:,1), tor1_mode_damp_freq_angle_60(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_60);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_60 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_60;
    plot(ellipse_tor1_mode_damp_freq_angle_60(:,1), ellipse_tor1_mode_damp_freq_angle_60(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_60_lb = min(ellipse_tor1_mode_damp_freq_angle_60(:, 2));
tor1_mode_freq_angle_60_ub = max(ellipse_tor1_mode_damp_freq_angle_60(:, 2));

%% GVT angle 90, freq OOP1
load('gvt_poles.mat');
modal_data_angle_90 = pl_all{1, 1};
figure; plot(modal_data_angle_90, 'k.')

r_min = 1;
r_max = 3;

oop1_mode_angle_90 = modal_data_angle_90(abs(modal_data_angle_90) >= r_min & abs(modal_data_angle_90) <= r_max);
oop1_mode_freq_angle_90 = abs(oop1_mode_angle_90);
oop1_mode_damp_angle_90 = -real(oop1_mode_angle_90)./abs(oop1_mode_angle_90);
oop1_mode_damp_freq_angle_90 = [oop1_mode_damp_angle_90(1:end), oop1_mode_freq_angle_90(1:end)];

mu_oop1_mode_damp_freq_angle_90 = mean(oop1_mode_damp_freq_angle_90, 1);        % 1×2 mean vector
Sigma_oop1_mode_damp_freq_angle_90 = cov(oop1_mode_damp_freq_angle_90);         % 2×2 covariance matrix

pdf_vals_oop1_mode_damp_freq_angle_90 = mvnpdf(oop1_mode_damp_freq_angle_90, mu_oop1_mode_damp_freq_angle_90, Sigma_oop1_mode_damp_freq_angle_90);

figure; hold on; axis equal
scatter(oop1_mode_damp_freq_angle_90(:,1), oop1_mode_damp_freq_angle_90(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop1_mode_damp_freq_angle_90);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop1_mode_damp_freq_angle_90 = (k * V * sqrt(D) * unit_circle)' + mu_oop1_mode_damp_freq_angle_90;
    plot(ellipse_oop1_mode_damp_freq_angle_90(:,1), ellipse_oop1_mode_damp_freq_angle_90(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop1_mode_freq_angle_90_lb = min(ellipse_oop1_mode_damp_freq_angle_90(:, 2));
oop1_mode_freq_angle_90_ub = max(ellipse_oop1_mode_damp_freq_angle_90(:, 2));

%% GVT angle 90, freq IP1
load('gvt_poles.mat');
modal_data_angle_90 = pl_all{1, 1};
figure; plot(modal_data_angle_90, 'k.')

r_min = 9;
r_max = 10;

ip1_mode_angle_90 = modal_data_angle_90(abs(modal_data_angle_90) >= r_min & abs(modal_data_angle_90) <= r_max);
ip1_mode_freq_angle_90 = abs(ip1_mode_angle_90);
ip1_mode_damp_angle_90 = -real(ip1_mode_angle_90)./abs(ip1_mode_angle_90);
ip1_mode_damp_freq_angle_90 = [ip1_mode_damp_angle_90(1:end), ip1_mode_freq_angle_90(1:end)];

mu_ip1_mode_damp_freq_angle_90 = mean(ip1_mode_damp_freq_angle_90, 1);        % 1×2 mean vector
Sigma_ip1_mode_damp_freq_angle_90 = cov(ip1_mode_damp_freq_angle_90);         % 2×2 covariance matrix

pdf_vals_ip1_mode_damp_freq_angle_90 = mvnpdf(ip1_mode_damp_freq_angle_90, mu_ip1_mode_damp_freq_angle_90, Sigma_ip1_mode_damp_freq_angle_90);

figure; hold on; axis equal
scatter(ip1_mode_damp_freq_angle_90(:,1), ip1_mode_damp_freq_angle_90(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_ip1_mode_damp_freq_angle_90);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_ip1_mode_damp_freq_angle_90 = (k * V * sqrt(D) * unit_circle)' + mu_ip1_mode_damp_freq_angle_90;
    plot(ellipse_ip1_mode_damp_freq_angle_90(:,1), ellipse_ip1_mode_damp_freq_angle_90(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

ip1_mode_freq_angle_90_lb = min(ellipse_ip1_mode_damp_freq_angle_90(:, 2));
ip1_mode_freq_angle_90_ub = max(ellipse_ip1_mode_damp_freq_angle_90(:, 2));

%% GVT angle 90, freq OOP2
load('gvt_poles.mat');
modal_data_angle_90 = pl_all{1, 1};
% figure; plot(modal_data_angle_90, 'k.')

r_min =11;
r_max = 12;

oop2_mode_angle_90 = modal_data_angle_90(abs(modal_data_angle_90) >= r_min & abs(modal_data_angle_90) <= r_max);
oop2_mode_freq_angle_90 = abs(oop2_mode_angle_90);
oop2_mode_damp_angle_90 = -real(oop2_mode_angle_90)./abs(oop2_mode_angle_90);
oop2_mode_damp_freq_angle_90 = [oop2_mode_damp_angle_90(1:end), oop2_mode_freq_angle_90(1:end)];

mu_oop2_mode_damp_freq_angle_90 = mean(oop2_mode_damp_freq_angle_90, 1);        % 1×2 mean vector
Sigma_oop2_mode_damp_freq_angle_90 = cov(oop2_mode_damp_freq_angle_90);         % 2×2 covariance matrix

pdf_vals_oop2_mode_damp_freq_angle_90 = mvnpdf(oop2_mode_damp_freq_angle_90, mu_oop2_mode_damp_freq_angle_90, Sigma_oop2_mode_damp_freq_angle_90);

figure; hold on; axis equal
scatter(oop2_mode_damp_freq_angle_90(:,1), oop2_mode_damp_freq_angle_90(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop2_mode_damp_freq_angle_90);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_oop2_mode_damp_freq_angle_90 = (k * V * sqrt(D) * unit_circle)' + mu_oop2_mode_damp_freq_angle_90;
    plot(ellipse_oop2_mode_damp_freq_angle_90(:,1), ellipse_oop2_mode_damp_freq_angle_90(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

oop2_mode_freq_angle_90_lb = min(ellipse_oop2_mode_damp_freq_angle_90(:, 2));
oop2_mode_freq_angle_90_ub = max(ellipse_oop2_mode_damp_freq_angle_90(:, 2));

%% GVT angle 90, freq TOR1
load('gvt_poles.mat');
modal_data_angle_90 = pl_all{1, 1};
figure; plot(modal_data_angle_90, 'k.')

r_min = 15;
r_max = 17;

tor1_mode_angle_90 = modal_data_angle_90(abs(modal_data_angle_90) >= r_min & abs(modal_data_angle_90) <= r_max);
tor1_mode_freq_angle_90 = abs(tor1_mode_angle_90);
tor1_mode_damp_angle_90 = -real(tor1_mode_angle_90)./abs(tor1_mode_angle_90);
tor1_mode_damp_freq_angle_90 = [tor1_mode_damp_angle_90(1:end), tor1_mode_freq_angle_90(1:end)];

mu_tor1_mode_damp_freq_angle_90 = mean(tor1_mode_damp_freq_angle_90, 1);        % 1×2 mean vector
Sigma_tor1_mode_damp_freq_angle_90 = cov(tor1_mode_damp_freq_angle_90);         % 2×2 covariance matrix

pdf_vals_tor1_mode_damp_freq_angle_90 = mvnpdf(tor1_mode_damp_freq_angle_90, mu_tor1_mode_damp_freq_angle_90, Sigma_tor1_mode_damp_freq_angle_90);

figure; hold on; axis equal
scatter(tor1_mode_damp_freq_angle_90(:,1), tor1_mode_damp_freq_angle_90(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_tor1_mode_damp_freq_angle_90);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse_tor1_mode_damp_freq_angle_90 = (k * V * sqrt(D) * unit_circle)' + mu_tor1_mode_damp_freq_angle_90;
    plot(ellipse_tor1_mode_damp_freq_angle_90(:,1), ellipse_tor1_mode_damp_freq_angle_90(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')

tor1_mode_freq_angle_90_lb = min(ellipse_tor1_mode_damp_freq_angle_90(:, 2));
tor1_mode_freq_angle_90_ub = max(ellipse_tor1_mode_damp_freq_angle_90(:, 2));

%%
oop1_mode_freq_lb = [oop1_mode_freq_angle_0_lb, oop1_mode_freq_angle_10_lb, oop1_mode_freq_angle_20_lb, oop1_mode_freq_angle_30_lb, oop1_mode_freq_angle_60_lb, oop1_mode_freq_angle_90_lb];
oop1_mode_freq_ub = [oop1_mode_freq_angle_0_ub, oop1_mode_freq_angle_10_ub, oop1_mode_freq_angle_20_ub, oop1_mode_freq_angle_30_ub, oop1_mode_freq_angle_60_ub, oop1_mode_freq_angle_90_ub];

oop2_mode_freq_lb = [oop2_mode_freq_angle_0_lb, oop2_mode_freq_angle_10_lb, oop2_mode_freq_angle_20_lb, oop2_mode_freq_angle_30_lb, oop2_mode_freq_angle_60_lb, oop2_mode_freq_angle_90_lb];
oop2_mode_freq_ub = [oop2_mode_freq_angle_0_ub, oop2_mode_freq_angle_10_ub, oop2_mode_freq_angle_20_ub, oop2_mode_freq_angle_30_ub, oop2_mode_freq_angle_60_ub, oop2_mode_freq_angle_90_ub];

ip1_mode_freq_lb = [ip1_mode_freq_angle_0_lb, ip1_mode_freq_angle_10_lb, ip1_mode_freq_angle_20_lb, ip1_mode_freq_angle_30_lb, ip1_mode_freq_angle_60_lb, ip1_mode_freq_angle_90_lb];
ip1_mode_freq_ub = [ip1_mode_freq_angle_0_ub, ip1_mode_freq_angle_10_ub, ip1_mode_freq_angle_20_ub, ip1_mode_freq_angle_30_ub, ip1_mode_freq_angle_60_ub, ip1_mode_freq_angle_90_ub];

tor1_mode_freq_lb = [tor1_mode_freq_angle_0_lb, tor1_mode_freq_angle_10_lb, tor1_mode_freq_angle_20_lb, tor1_mode_freq_angle_30_lb, tor1_mode_freq_angle_60_lb, tor1_mode_freq_angle_90_lb];
tor1_mode_freq_ub = [tor1_mode_freq_angle_0_ub, tor1_mode_freq_angle_10_ub, tor1_mode_freq_angle_20_ub, tor1_mode_freq_angle_30_ub, tor1_mode_freq_angle_60_ub, tor1_mode_freq_angle_90_ub];

%%
save('GVT_freq_lower_and_upper_bounds.mat', 'oop1_mode_freq_lb', 'oop1_mode_freq_ub', 'oop2_mode_freq_lb', 'oop2_mode_freq_ub', 'ip1_mode_freq_lb', 'ip1_mode_freq_ub', 'tor1_mode_freq_lb', 'tor1_mode_freq_ub')

%%
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
beta_y_g1_1_experimental_data_set = exprData{1, 1}.beta_y(2:end);  
beta_y_g1_1_lb = 0.99*beta_y_g1_1_experimental_data_set-0.01*;


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