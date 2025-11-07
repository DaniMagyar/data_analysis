% BAfc_figure_4_supp_latency.m
% Supplementary figure: Onset latency comparison between selective and multisensory neurons
% For each brain region (LA and Astria) and stimulus type (CS and US)

clear all; close all

%% Setup
recordings = {...
    'MD292_002_kilosort','MD293_001_kilosort','MD294_001_kilosort','MD295_001_kilosort',...
    'MD296_001_kilosort','MD297_001_kilosort','MD298_001_kilosort','MD299_001_kilosort',...
    'MD300_001_kilosort','MD304_001_kilosort','MD305_001_kilosort','MD307_001_kilosort',...
    'MD309_001_kilosort','MD310_001_kilosort','MD311_002_kilosort','MD312_001_kilosort',...
    'MD313_001_kilosort','MD314_001_kilosort','MD315_001_kilosort','MD316_002_kilosort',...
    'MD317_001_kilosort','MD318_001_kilosort','MD318_002_kilosort','MD319_003_kilosort'};

ttl = {'triptest_sound_only','triptest_shocks_only','triptest_both'};
brain_regions = {'LA', 'Astria'};

cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
cell_metrics = BAfc_match_celltypes('cell_metrics', cell_metrics);

g.colors = BAfc_colors;
g.fontSize1 = 14;
g.fontSize2 = 12;
g.bin_time = 0.001;
g.pre_time = 5;
g.post_time = 0.5;
g.monosyn_window = 0.025;  % 0-25ms
g.smoothvalue = 7;

g.onset_threshold = 3;
g.min_consec_bins = max(1, round(0.001 / g.bin_time));
g.alpha = 0.0;

% Responsiveness detection method
g.use_two_rule = true;  % true: two-rule (Rule 1 OR Rule 2), false: one-rule (z-score only)

% Two-rule responsiveness parameters (used if g.use_two_rule = true)
g.zscore_threshold_rule1 = 3;   % Rule 1: z-score threshold
g.prob_threshold_rule1 = 0.25;  % Rule 1: probability threshold
g.zscore_threshold_rule2 = 10;  % Rule 2: z-score threshold
g.prob_threshold_rule2 = 0.1;   % Rule 2: probability threshold

% One-rule responsiveness parameter (used if g.use_two_rule = false)
g.zscore_threshold_one_rule = 5;  % One-rule: z-score threshold only

%% Calculate PSTHs once
fprintf('Calculating PSTHs...\n');
psthZ_full = cell(1,3);
baseline_bins = round(g.pre_time / g.bin_time);
for hmp = 1:3
    psth_spx_og = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{hmp}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Z-score using baseline period only
    baseline_mean = mean(psth_spx_og(:, 1:baseline_bins), 2);
    baseline_std = std(psth_spx_og(:, 1:baseline_bins), 0, 2);
    baseline_std(baseline_std == 0) = 1;  % Avoid division by zero
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;

    psth_spx = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);
    psthZ_full{hmp} = psth_spx;
end

%% Monosynaptic detection for each brain region
results_all = cell(1,2);
monosyn_window_bins = round((g.pre_time)/g.bin_time+1 : (g.pre_time+g.monosyn_window)/g.bin_time);

for br = 1:2
    fprintf('\nProcessing %s...\n', brain_regions{br});
    idx_neurons = strcmp(cell_metrics.brainRegion, brain_regions{br});

    % For LA, use PNs only; for Astria, use all
    if strcmp(brain_regions{br}, 'LA')
        idx_PN = idx_neurons & strcmp(cell_metrics.putativeCellType, 'PN');
    else
        idx_PN = idx_neurons;
    end

    n_PN = sum(idx_PN);

    if n_PN == 0
        continue;
    end

    % Get number of trials
    num_trials_CS = size(cell_metrics.general.(ttl{1}){1}, 1);
    num_trials_US = size(cell_metrics.general.(ttl{2}){1}, 1);

    % Extract monosynaptic window responses for CS, US
    psth_CS_PN = psthZ_full{1}(idx_PN, :);
    psth_US_PN = psthZ_full{2}(idx_PN, :);

    % Calculate spike probabilities trial-by-trial
    [~, ~, postAP_norm_CS] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{1}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    [~, ~, postAP_norm_US] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{2}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Process PNs
    fprintf('  Processing neurons (%d)...\n', n_PN);
    CS_peak_PN = max(psth_CS_PN(:, monosyn_window_bins), [], 2);
    US_peak_PN = max(psth_US_PN(:, monosyn_window_bins), [], 2);

    neuron_indices_PN = find(idx_PN);
    CS_prob_PN = zeros(n_PN, 1);
    US_prob_PN = zeros(n_PN, 1);

    for n = 1:n_PN
        global_idx = neuron_indices_PN(n);

        % CS probability - count trials with at least 1 spike in 0-25ms window
        if ~isempty(postAP_norm_CS{global_idx})
            responsive_trials_CS = 0;
            for trial = 1:length(postAP_norm_CS{global_idx})
                trial_spikes = postAP_norm_CS{global_idx}{trial};
                if any(trial_spikes > 0 & trial_spikes <= g.monosyn_window)
                    responsive_trials_CS = responsive_trials_CS + 1;
                end
            end
            CS_prob_PN(n) = responsive_trials_CS / num_trials_CS;
        end

        % US probability
        if ~isempty(postAP_norm_US{global_idx})
            responsive_trials_US = 0;
            for trial = 1:length(postAP_norm_US{global_idx})
                trial_spikes = postAP_norm_US{global_idx}{trial};
                if any(trial_spikes > 0 & trial_spikes <= g.monosyn_window)
                    responsive_trials_US = responsive_trials_US + 1;
                end
            end
            US_prob_PN(n) = responsive_trials_US / num_trials_US;
        end
    end

    % Responsiveness for PNs
    if g.use_two_rule
        % Two-rule: (Rule 1 OR Rule 2)
        CS_responsive_PN = (CS_peak_PN >= g.zscore_threshold_rule1 & CS_prob_PN >= g.prob_threshold_rule1) | ...
                           (CS_peak_PN >= g.zscore_threshold_rule2 & CS_prob_PN >= g.prob_threshold_rule2);
        US_responsive_PN = (US_peak_PN >= g.zscore_threshold_rule1 & US_prob_PN >= g.prob_threshold_rule1) | ...
                           (US_peak_PN >= g.zscore_threshold_rule2 & US_prob_PN >= g.prob_threshold_rule2);
    else
        % One-rule: z-score only
        CS_responsive_PN = CS_peak_PN >= g.zscore_threshold_one_rule;
        US_responsive_PN = US_peak_PN >= g.zscore_threshold_one_rule;
    end

    fprintf('  CS responsive: %d, US responsive: %d\n', sum(CS_responsive_PN), sum(US_responsive_PN));

    % Categorize neurons
    Clusters_PN = zeros(n_PN, 1);
    Clusters_PN(CS_responsive_PN & ~US_responsive_PN) = 1;  % CS-selective
    Clusters_PN(US_responsive_PN & ~CS_responsive_PN) = 2;  % US-selective
    Clusters_PN(CS_responsive_PN & US_responsive_PN) = 3;   % Multisensory

    % Calculate onset latencies
    CS_onset_lat_PN = nan(n_PN, 1);
    US_onset_lat_PN = nan(n_PN, 1);

    for n = 1:n_PN
        [CS_onset_lat_PN(n), ~] = compute_onset_offset_latency(psth_CS_PN(n, :), monosyn_window_bins, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [US_onset_lat_PN(n), ~] = compute_onset_offset_latency(psth_US_PN(n, :), monosyn_window_bins, g.onset_threshold, g.min_consec_bins, g.bin_time);
    end

    fprintf('  CS onset latencies (ms): valid=%d\n', sum(~isnan(CS_onset_lat_PN)));
    fprintf('  US onset latencies (ms): valid=%d\n', sum(~isnan(US_onset_lat_PN)));

    % Store results
    results_all{br}.Clusters_PN = Clusters_PN;
    results_all{br}.CS_onset_lat_PN = CS_onset_lat_PN;
    results_all{br}.US_onset_lat_PN = US_onset_lat_PN;

    fprintf('  CS-sel: %d | US-sel: %d | Multi: %d\n', ...
        sum(Clusters_PN==1), sum(Clusters_PN==2), sum(Clusters_PN==3));
end

%% Create supplementary figure - selective vs multisensory comparison
fig = figure('Units', 'pixels', 'Position', [100, 100, 800, 400], 'Visible', 'on');

% Collect latency data for each group
% LA CS-selective and Multisensory (CS latency)
LA_CS_sel_lat = results_all{1}.CS_onset_lat_PN(results_all{1}.Clusters_PN == 1) * 1000;  % Convert to ms
LA_Multi_CS_lat = results_all{1}.CS_onset_lat_PN(results_all{1}.Clusters_PN == 3) * 1000;

% Astria CS-selective and Multisensory (CS latency)
Astria_CS_sel_lat = results_all{2}.CS_onset_lat_PN(results_all{2}.Clusters_PN == 1) * 1000;
Astria_Multi_CS_lat = results_all{2}.CS_onset_lat_PN(results_all{2}.Clusters_PN == 3) * 1000;

% LA US-selective and Multisensory (US latency)
LA_US_sel_lat = results_all{1}.US_onset_lat_PN(results_all{1}.Clusters_PN == 2) * 1000;
LA_Multi_US_lat = results_all{1}.US_onset_lat_PN(results_all{1}.Clusters_PN == 3) * 1000;

% Astria US-selective and Multisensory (US latency)
Astria_US_sel_lat = results_all{2}.US_onset_lat_PN(results_all{2}.Clusters_PN == 2) * 1000;
Astria_Multi_US_lat = results_all{2}.US_onset_lat_PN(results_all{2}.Clusters_PN == 3) * 1000;

% Remove NaN values
LA_CS_sel_lat = LA_CS_sel_lat(~isnan(LA_CS_sel_lat));
LA_Multi_CS_lat = LA_Multi_CS_lat(~isnan(LA_Multi_CS_lat));
Astria_CS_sel_lat = Astria_CS_sel_lat(~isnan(Astria_CS_sel_lat));
Astria_Multi_CS_lat = Astria_Multi_CS_lat(~isnan(Astria_Multi_CS_lat));
LA_US_sel_lat = LA_US_sel_lat(~isnan(LA_US_sel_lat));
LA_Multi_US_lat = LA_Multi_US_lat(~isnan(LA_Multi_US_lat));
Astria_US_sel_lat = Astria_US_sel_lat(~isnan(Astria_US_sel_lat));
Astria_Multi_US_lat = Astria_Multi_US_lat(~isnan(Astria_Multi_US_lat));

hold on;

% Positions for 4 pairs of boxplots
positions = [1 2, 4 5, 7 8, 10 11];
all_data = {LA_CS_sel_lat, LA_Multi_CS_lat, Astria_CS_sel_lat, Astria_Multi_CS_lat, ...
            LA_US_sel_lat, LA_Multi_US_lat, Astria_US_sel_lat, Astria_Multi_US_lat};

% Plot boxplots
boxplot([all_data{1}; all_data{2}; all_data{3}; all_data{4}; ...
         all_data{5}; all_data{6}; all_data{7}; all_data{8}], ...
        [ones(length(all_data{1}),1)*1; ones(length(all_data{2}),1)*2; ...
         ones(length(all_data{3}),1)*4; ones(length(all_data{4}),1)*5; ...
         ones(length(all_data{5}),1)*7; ones(length(all_data{6}),1)*8; ...
         ones(length(all_data{7}),1)*10; ones(length(all_data{8}),1)*11], ...
        'Positions', positions, 'Widths', 0.6, 'Colors', 'k');

% Overlay scatter points with jitter
jitter_amount = 0.15;
for i = 1:8
    if ~isempty(all_data{i})
        x_jitter = positions(i) + (rand(length(all_data{i}), 1) - 0.5) * jitter_amount;

        % Color based on group type (selective vs multisensory)
        if mod(i, 2) == 1  % Selective groups (odd indices)
            color = [0.5 0.5 0.5];  % Gray
        else  % Multisensory groups (even indices)
            color = [0.6 0.2 0.6];  % Purple
        end

        scatter(x_jitter, all_data{i}, 20, color, 'filled', 'MarkerFaceAlpha', 0.4);
    end
end

hold off;

% Formatting
xlim([0 12]);
ylim([0 50]);
xticks([1.5 4.5 7.5 10.5]);
xticklabels({'LA CS', 'Astria CS', 'LA US', 'Astria US'});
ylabel('Onset Latency (ms)', 'FontSize', g.fontSize2);
title('Onset Latencies: Selective vs Multisensory', 'FontSize', g.fontSize1, 'FontWeight', 'bold');
set(gca, 'FontSize', g.fontSize2);
box off;

% Add statistics - ranksum test for each pair (selective vs multisensory)
y_pos_sig = 45;  % Position for significance markers
pairs_to_test = {
    {LA_CS_sel_lat, LA_Multi_CS_lat, 1.5};      % LA CS pair, x-position for text
    {Astria_CS_sel_lat, Astria_Multi_CS_lat, 4.5};  % Astria CS pair
    {LA_US_sel_lat, LA_Multi_US_lat, 7.5};      % LA US pair
    {Astria_US_sel_lat, Astria_Multi_US_lat, 10.5};  % Astria US pair
};

hold on;
for p = 1:4
    data1 = pairs_to_test{p}{1};
    data2 = pairs_to_test{p}{2};
    x_center = pairs_to_test{p}{3};

    if ~isempty(data1) && ~isempty(data2) && length(data1) > 0 && length(data2) > 0
        [p_val, ~] = ranksum(data1, data2);

        % Draw line connecting the two boxplots
        x_left = x_center - 0.5;
        x_right = x_center + 0.5;
        plot([x_left x_right], [y_pos_sig y_pos_sig], 'k-', 'LineWidth', 1.5);

        if p_val < 0.05
            sig_text = get_sig_stars(p_val);
        else
            sig_text = 'n.s.';
        end

        text(x_center, y_pos_sig, sig_text, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
    end
end
hold off;

% Add legend
legend_ax = axes('Position', [0.75 0.8 0.15 0.1], 'Visible', 'off');
hold(legend_ax, 'on');
scatter(legend_ax, 1, 1, 50, [0.5 0.5 0.5], 'filled');
scatter(legend_ax, 1, 2, 50, [0.6 0.2 0.6], 'filled');
legend(legend_ax, {'Selective', 'Multisensory'}, 'Location', 'best', 'FontSize', g.fontSize2);
hold(legend_ax, 'off');

fprintf('\nDone.\n');

%% Helper functions
function sig_text = get_sig_stars(p_value)
    if p_value < 0.001
        sig_text = '***';
    elseif p_value < 0.01
        sig_text = '**';
    elseif p_value < 0.05
        sig_text = '*';
    else
        sig_text = '';
    end
end

function [onset_lat, offset_lat] = compute_onset_offset_latency(z_trace, event_inds, threshold, min_consec, bin_time)
    seg = z_trace(event_inds);
    if any(isnan(seg))
        onset_lat = NaN;
        offset_lat = NaN;
        return
    end

    % Exclude artifact period (0-12ms)
    artifact_period = 0.012;  % 12ms
    artifact_bins = round(artifact_period / bin_time);

    % Only search for onset after artifact period
    if length(seg) <= artifact_bins
        onset_lat = NaN;
        offset_lat = NaN;
        return
    end

    seg_post_artifact = seg(artifact_bins+1:end);
    isAbove = seg_post_artifact >= threshold;
    winsum = conv(double(isAbove), ones(1, min_consec), 'same');
    onset_idx_relative = find(winsum >= min_consec, 1, 'first');

    if isempty(onset_idx_relative)
        onset_lat = NaN;
        offset_lat = NaN;
    else
        % Convert onset_idx back to absolute time (including artifact period)
        onset_idx = artifact_bins + onset_idx_relative;
        onset_lat = (onset_idx - 1) * bin_time;

        seg_after_onset = seg(onset_idx:end);
        isBelow = seg_after_onset < threshold;
        winsum_below = conv(double(isBelow), ones(1, min_consec), 'same');
        offset_idx_relative = find(winsum_below >= min_consec, 1, 'first');

        if isempty(offset_idx_relative)
            % No offset detected - use end of response window
            offset_lat = (length(event_inds) - 1) * bin_time;
        else
            offset_idx = onset_idx + offset_idx_relative - 1;
            offset_lat = (offset_idx - 1) * bin_time;
        end
    end
end
