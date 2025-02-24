function [dResp1, dResp2] =  M2_NDNF_VIP2R_Arch_calc_dResponses

% input: index of significantly responsive units, either exc or inh

pre_time = 5;
post_time = 5;
bin_time = 0.001;

[psth_spx_1, psth_spx_2, idx] =  M2_NDNF_VIP2R_Arch_find_shocks_v01(pre_time, post_time, bin_time);

baselines1 = sum(psth_spx_1(:, pre_time/bin_time-500:pre_time/bin_time),2);
responses1 = sum(psth_spx_1(:, pre_time/bin_time+1:pre_time/bin_time+500),2);
delta1 = responses1(idx) - baselines1(idx);

baselines2 = sum(psth_spx_2(:, pre_time/bin_time-500:pre_time/bin_time),2);
responses2 = sum(psth_spx_2(:, pre_time/bin_time+1:pre_time/bin_time+500),2);
delta2 = responses2(idx) - baselines2(idx);


figure;
hold on;
for i = 1:numel(idx)
    % Plot individual data points
    plot([1, 2], [delta1(i), delta2(i)], 'k-'); % Line connecting paired points
    plot(1, delta1(i), 'ro'); % Data point from vector1
    plot(2, delta2(i), 'bo'); % Data point from vector2
end
hold off;

% Combine the data into one matrix
data = [delta1, delta2];

% Create boxplots
figure;
boxplot(data, 'Labels', {'Vector 1', 'Vector 2'});
title('Boxplots of Two Vectors');
ylabel('Values');

ylim([-2 10])