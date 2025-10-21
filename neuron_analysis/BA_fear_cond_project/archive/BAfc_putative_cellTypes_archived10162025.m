function [cell_metrics] = BAfc_putative_cellTypes(varargin)
% BAfc_simple_clustering.m
% Simple k-means clustering to separate high firing rate, narrow waveform neurons
% Uses firingRate, troughToPeak, and waveform half-width
%clear all, close all
prs = inputParser;
addParameter(prs,'cell_metrics',[],@isstruct)
addParameter(prs,'plot',false,@islogical) % plot results
addParameter(prs,'width_critical',0.4,@isnumeric) % border between narrow and wide waveforms
addParameter(prs,'fr_critical',10,@isnumerci) % border between fast and regular firing rate
% Bienvenu cikkben a 10Hz jo tampont a PV IN-ekre
addParameter(prs,'ttp_d_crit',0.16,@isnumerci) % border between narrow and wide
addParameter(prs,'ab_ratio',0.1,@isnumeric)
parse(prs,varargin{:})
g = prs.Results;

cell_metrics = g.cell_metrics;

num_neurons_total = numel(cell_metrics.cellID);

% Extract ALL data into separate variables (before any filtering)
% This keeps g.cell_metrics intact as reference
all_cellID = cell_metrics.cellID;
all_sessionName = cell_metrics.sessionName;
all_firingRate = cell_metrics.firingRate';
all_troughToPeak = cell_metrics.troughToPeak';
all_ab_ratio = cell_metrics.ab_ratio';
all_acg_tau_decay = cell_metrics.acg_tau_decay';
all_acg_tau_rise = cell_metrics.acg_tau_rise';
all_polarity = cell_metrics.polarity';
all_waveforms = cell2mat(cell_metrics.waveforms.filt');
all_acg_wide = cell_metrics.acg.wide;  % 1001 x num_neurons
all_isi_log10 = cell_metrics.isi.log10;  % 100 x num_neurons

% Initialize index tracking
% original_idx will always map current neurons back to the original dataset
original_idx = (1:num_neurons_total)';  % Start with all neurons

% Remove inverted polarity neurons
is_inverted = all_polarity > 0;
num_inverted = sum(is_inverted);

% Keep only standard polarity neurons
keep_idx = ~is_inverted;
original_idx = original_idx(keep_idx);  % Update index tracking

% Filter all variables
cellID = all_cellID(keep_idx);
sessionName = all_sessionName(keep_idx);
firingRate = all_firingRate(keep_idx);
troughToPeak = all_troughToPeak(keep_idx);
ab_ratio = all_ab_ratio(keep_idx);
acg_tau_decay = all_acg_tau_decay(keep_idx);
acg_tau_rise = all_acg_tau_rise(keep_idx);
polarity = all_polarity(keep_idx);
waveforms = all_waveforms(keep_idx, :);
acg_wide = all_acg_wide(:, keep_idx);
isi_log10 = all_isi_log10(:, keep_idx);

num_neurons = sum(keep_idx);

%% ===== PREPARE FEATURES =====
% Perform PCA on ACG data
[coeff_acg, score_acg, ~, ~, explained_acg] = pca(acg_wide');  % Transpose so neurons are rows
acg_pc1 = score_acg(:, 1);
acg_pc2 = score_acg(:, 2);

fprintf('  PC1 explains %.1f%% of variance\n', explained_acg(1));
fprintf('  PC2 explains %.1f%% of variance\n', explained_acg(2));
fprintf('  PC1+PC2 explain %.1f%% of variance\n\n', sum(explained_acg(1:2)));

% Perform PCA on waveforms
[coeff_wf, score_wf, ~, ~, explained_wf] = pca(waveforms);  % Waveforms already has neurons as rows
waveform_pc1 = score_wf(:, 1);
waveform_pc2 = score_wf(:, 2);

% Perform PCA on ISI log10
[coeff_isi, score_isi, ~, ~, explained_isi] = pca(isi_log10');  % Transpose so neurons are rows
isi_pc1 = score_isi(:, 1);
isi_pc2 = score_isi(:, 2);


%% ===== PREPARE FOR CLUSTERING =====

% INDEX TRACKING NOTE:
% After removing inverted neurons and NaN values, the vectors shrink.
% The variable 'original_idx_valid' maintains the mapping:
%   original_idx_valid(i) = original neuron number for valid neuron i
% This allows tracing any result back to the original loaded dataset.

% Remove neurons with missing data
valid_idx = ~isnan(firingRate) & ~isnan(troughToPeak) & ~isnan(ab_ratio) & ...
            ~isnan(acg_tau_decay) & ~isnan(acg_tau_rise) & ...
            ~isnan(acg_pc1) & ~isnan(acg_pc2) & ...
            ~isnan(waveform_pc1) & ~isnan(waveform_pc2) & ...
            ~isnan(isi_pc1) & ~isnan(isi_pc2);
num_valid = sum(valid_idx);

% Update index tracking for valid neurons
original_idx_valid = original_idx(valid_idx);

% Extract valid data for clustering
valid_cellID = cellID(valid_idx);
valid_sessionName = sessionName(valid_idx);
valid_fr = firingRate(valid_idx);
valid_ttp = troughToPeak(valid_idx);
valid_ab = ab_ratio(valid_idx);
valid_tau_decay = acg_tau_decay(valid_idx);
valid_tau_rise = acg_tau_rise(valid_idx);
valid_acg_pc1 = acg_pc1(valid_idx);
valid_acg_pc2 = acg_pc2(valid_idx);
valid_waveform_pc1 = waveform_pc1(valid_idx);
valid_waveform_pc2 = waveform_pc2(valid_idx);
valid_isi_pc1 = isi_pc1(valid_idx);
valid_isi_pc2 = isi_pc2(valid_idx);

% features = [zscore(valid_fr), ...
%            zscore(valid_ttp), ...
%            zscore(valid_ab), ...
%            zscore(valid_acg_pc1), ...
%            zscore(valid_acg_pc2), ...
%            zscore(valid_waveform_pc1), ...
%            zscore(valid_waveform_pc2)];
features = [zscore(valid_fr), ...
           zscore(valid_ttp), ...
           zscore(valid_ab), ...
           ];

%% ===== K-MEANS CLUSTERING =====
[cluster_idx, centroids] = kmeans(features, 2, 'Replicates', 50, 'MaxIter', 1000);

% Determine which cluster has higher firing rate (IN)
cluster1_fr = mean(valid_fr(cluster_idx == 1));
cluster2_fr = mean(valid_fr(cluster_idx == 2));

if cluster1_fr > cluster2_fr
    IN_cluster = 1;
    PN_cluster = 2;
else
    IN_cluster = 2;
    PN_cluster = 1;
end

% Create initial two-way classification
is_IN = (cluster_idx == IN_cluster);
is_PN = (cluster_idx == PN_cluster);
is_unknown = false(num_valid, 1);  % Initialize empty unknown array

% % Apply manual criteria: INs with negative ab_ratio cannot be true INs
% % These are likely false positives and should be marked as unknown
% manual_exclusion = is_IN & (valid_ttp > 0.5);
% num_manual_excluded = sum(manual_exclusion);
% 
% if num_manual_excluded > 0
%     % Move manually excluded neurons from IN to unknown
%     is_IN(manual_exclusion) = false;
%     is_unknown(manual_exclusion) = true;
% end


%% ===== VISUALIZATION =====

fprintf('Creating visualization...\n\n');

colors = BAfc_colors();
color_IN = colors.IN_primary;
color_PN = colors.PN_primary;
color_unknown = [0.6, 0.6, 0.6];  % Grey for unknown

fig = figure('Position', [50, 50, 1800, 1000]);

% (1,1) Firing Rate vs AB Ratio
subplot(3, 3, 1);
hold on;
scatter(valid_ab(is_PN), valid_fr(is_PN), 50, color_PN, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
scatter(valid_ab(is_IN), valid_fr(is_IN), 50, color_IN, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
xlabel('AB Ratio (Waveform Asymmetry)');
ylabel('Firing Rate (Hz)');
set(gca, 'YScale', 'log');
title('Firing Rate vs AB Ratio');
legend({'PN', 'IN'}, 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

% (1,2) Firing Rate vs Trough-to-Peak
subplot(3, 3, 2);
hold on;
scatter(valid_ttp(is_PN), valid_fr(is_PN), 50, color_PN, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
scatter(valid_ttp(is_IN), valid_fr(is_IN), 50, color_IN, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
xlabel('Trough-to-Peak (ms)');
ylabel('Firing Rate (Hz)');
set(gca, 'YScale', 'log');
title('Firing Rate vs Trough-to-Peak');
legend({'PN', 'IN'}, 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

% (1,3) AB Ratio vs Trough-to-Peak
subplot(3, 3, 3);
hold on;
scatter(valid_ttp(is_PN), valid_ab(is_PN), 50, color_PN, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
scatter(valid_ttp(is_IN), valid_ab(is_IN), 50, color_IN, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
xlabel('Trough-to-Peak (ms)');
ylabel('AB Ratio');
title('AB Ratio vs Trough-to-Peak');
legend({'PN', 'IN'}, 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

% (2,1) Firing Rate Cumulative Distribution
subplot(3, 3, 4);
hold on;
% Sort firing rates for cumulative plot
[sorted_fr_pn, ~] = sort(valid_fr(is_PN));
cdf_pn = (1:length(sorted_fr_pn))' / length(sorted_fr_pn);
plot(sorted_fr_pn, cdf_pn, 'LineWidth', 2.5, 'Color', color_PN);
[sorted_fr_in, ~] = sort(valid_fr(is_IN));
cdf_in = (1:length(sorted_fr_in))' / length(sorted_fr_in);
plot(sorted_fr_in, cdf_in, 'LineWidth', 2.5, 'Color', color_IN);
xlabel('Firing Rate (Hz)');
ylabel('Cumulative Probability');
title('Firing Rate CDF');
legend({'PN', 'IN'}, 'Location', 'southeast');
set(gca, 'XScale', 'log');
grid on;
set(gca, 'FontSize', 11);
ylim([0 1]);

% (2,2) AB Ratio Distribution
subplot(3, 3, 5);
hold on;
histogram(valid_ab(is_PN), 30, 'FaceColor', color_PN, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
histogram(valid_ab(is_IN), 30, 'FaceColor', color_IN, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xlabel('AB Ratio');
ylabel('Count');
title('AB Ratio Distribution');
legend({'PN', 'IN'}, 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

% (2,3) Trough-to-Peak Distribution
subplot(3, 3, 6);
hold on;
histogram(valid_ttp(is_PN), 30, 'FaceColor', color_PN, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
histogram(valid_ttp(is_IN), 30, 'FaceColor', color_IN, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xlabel('Trough-to-Peak (ms)');
ylabel('Count');
title('Trough-to-Peak Distribution');
legend({'PN', 'IN'}, 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

% (3,1-3) Statistics Summary spanning 3 columns
subplot(3, 3, [7, 8, 9]);
axis off;

% Calculate total unknown (manual exclusions + inverted)
total_unknown = sum(is_unknown) + num_inverted;

text(0.1, 0.95, 'Clustering Summary', 'FontSize', 14, 'FontWeight', 'bold');
text(0.1, 0.85, 'Method: k-means (k=2) + manual criteria', 'FontSize', 10);
text(0.1, 0.78, sprintf('Total neurons: %d', num_neurons_total), 'FontSize', 10);

text(0.1, 0.60, sprintf('Interneurons: %d (%.1f%%)', sum(is_IN), 100*sum(is_IN)/num_neurons_total), ...
    'FontSize', 11, 'Color', color_IN, 'FontWeight', 'bold');
text(0.1, 0.53, sprintf('  FR: %.1f ± %.1f Hz', mean(valid_fr(is_IN)), std(valid_fr(is_IN))), ...
    'FontSize', 8, 'Color', color_IN);
text(0.1, 0.48, sprintf('  TTP: %.3f ± %.3f ms', mean(valid_ttp(is_IN)), std(valid_ttp(is_IN))), ...
    'FontSize', 8, 'Color', color_IN);
text(0.1, 0.43, sprintf('  AB: %.3f ± %.3f', mean(valid_ab(is_IN)), std(valid_ab(is_IN))), ...
    'FontSize', 8, 'Color', color_IN);

text(0.1, 0.32, sprintf('Principal Neurons: %d (%.1f%%)', sum(is_PN), 100*sum(is_PN)/num_neurons_total), ...
    'FontSize', 11, 'Color', color_PN, 'FontWeight', 'bold');
text(0.1, 0.25, sprintf('  FR: %.1f ± %.1f Hz', mean(valid_fr(is_PN)), std(valid_fr(is_PN))), ...
    'FontSize', 8, 'Color', color_PN);
text(0.1, 0.20, sprintf('  TTP: %.3f ± %.3f ms', mean(valid_ttp(is_PN)), std(valid_ttp(is_PN))), ...
    'FontSize', 8, 'Color', color_PN);
text(0.1, 0.15, sprintf('  AB: %.3f ± %.3f', mean(valid_ab(is_PN)), std(valid_ab(is_PN))), ...
    'FontSize', 8, 'Color', color_PN);

text(0.1, 0.04, sprintf('Unknown: %d (%.1f%%)', total_unknown, 100*total_unknown/num_neurons_total), ...
    'FontSize', 11, 'Color', color_unknown, 'FontWeight', 'bold');
text(0.1, -0.03, sprintf('  Inverted polarity: %d', num_inverted), ...
    'FontSize', 8, 'Color', color_unknown);
if sum(is_unknown) > 0
    text(0.1, -0.08, sprintf('  Manual exclusions: %d (IN with TTP>0.5)', sum(is_unknown)), ...
        'FontSize', 8, 'Color', color_unknown);
end

sgtitle('Neuron Classification: k-means + Manual Criteria', 'FontSize', 15, 'FontWeight', 'bold');

%% ===== SAVE RESULTS =====
% Create cell type assignments in the original g.cell_metrics structure
cell_metrics.putativeCellType(1:end) = {'unknown'};
cell_metrics.putativeCellType(original_idx_valid(is_IN)) = {'IN'};
cell_metrics.putativeCellType(original_idx_valid(is_PN)) = {'PN'};
% Outliers remain as 'unknown'

% Create full-length cell type array (for all neurons, not just valid)
all_cellType = cell(num_neurons_total, 1);
all_cellType(:) = {'unknown'};
all_cellType(original_idx_valid(is_IN)) = {'IN'};
all_cellType(original_idx_valid(is_PN)) = {'PN'};
% Outliers remain as 'unknown'
