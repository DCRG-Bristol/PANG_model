%%
load('WTTests\WT_testData\SJD_wing_proc_modalData.mat');
modal_data_idx_angle_5_idx_speed_6 = modalData.stabPls{5, 6};
figure; plot(modal_data_idx_angle_5_idx_speed_6, 'k.')

r_min = 1;
r_max = 3;

oop1_mode = modal_data_idx_angle_5_idx_speed_6(abs(modal_data_idx_angle_5_idx_speed_6) >= r_min & abs(modal_data_idx_angle_5_idx_speed_6) <= r_max);
oop1_mode_freq = abs(oop1_mode);
oop1_mode_damp = -real(oop1_mode)./abs(oop1_mode);
oop1_mode_damp_freq = [oop1_mode_damp, oop1_mode_freq];

mu_oop1_mode_damp_freq = mean(oop1_mode_damp_freq, 1);        % 1×2 mean vector
Sigma_oop1_mode_damp_freq = cov(oop1_mode_damp_freq);         % 2×2 covariance matrix

pdf_vals_oop1_mode_damp_freq = mvnpdf(oop1_mode_damp_freq, mu_oop1_mode_damp_freq, Sigma_oop1_mode_damp_freq);

figure; hold on; axis equal
scatter(oop1_mode_damp_freq(:,1), oop1_mode_damp_freq(:,2), 10, 'filled')

theta = linspace(0,2*pi,200);
unit_circle = [cos(theta); sin(theta)];

[V,D] = eig(Sigma_oop1_mode_damp_freq);

for k = 1:3   % 1σ, 2σ, 3σ
    ellipse = (k * V * sqrt(D) * unit_circle)' + mu_oop1_mode_damp_freq;
    plot(ellipse(:,1), ellipse(:,2), 'LineWidth', 2)
end

legend('Data','1σ','2σ','3σ')


K = 2;   % number of Gaussians
gm = fitgmdist(oop1_mode_damp_freq, K);
idx = cluster(gm, oop1_mode_damp_freq);


p = pdf(gm, oop1_mode_damp_freq);            % total mixture density
p_k = posterior(gm, oop1_mode_damp_freq);    % N×K soft assignments


figure; hold on; axis equal
scatter(oop1_mode_damp_freq(:,1), oop1_mode_damp_freq(:,2), 8, 'k', 'filled')

theta = linspace(0,2*pi,200);
circle = [cos(theta); sin(theta)];

for k = 1:K
    [V,D] = eig(gm.Sigma(:,:,k));
    ellipse = (V * sqrt(D) * circle)' + gm.mu(k,:);
    plot(ellipse(:,1), ellipse(:,2), 'LineWidth', 2)
end

