% BAfc_figure_5_v2.m
% Optogenetic manipulation of monosynaptic responses for LA and Astria
% Compares CS+US vs CS+US+light in the monosynaptic window (0-25ms)
% Detects neurons responsive in EITHER condition (no-light OR light)
% To run gui: BAfc_monosyn_raster_ui(g_ui, monosyn_results, ttl)

clear all; close all

%% Setup - recordings with optogenetic manipulation
recordings = {...
  'MD298_001_kilosort',...
  'MD299_001_kilosort',...
  'MD300_001_kilosort',...
  'MD304_001_kilosort',...
  'MD307_001_kilosort',...
  'MD309_001_kilosort',...
  'MD315_001_kilosort',...
  'MD316_002_kilosort',...
  'MD317_001_kilosort',...
  'MD318_001_kilosort',...
  'MD318_002_kilosort',...
  'MD319_003_kilosort'};

% TTL types: CS+US only vs CS+US+light
% ttl = {'triptest_both', 'triptest_both_light'};
% ttl = {'triptest_sound_only', 'triptest_sound_only_light'};
ttl = {'triptest_shocks_only', 'triptest_shocks_only_light'};
brain_regions = {'LA', 'Astria'};

cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
cell_metrics = BAfc_match_celltypes('cell_metrics', cell_metrics);

g.colors = BAfc_colors;
g.fontSize1 = 14;
g.fontSize2 = 12;
g.bin_time = 0.001;
g.pre_time = 5;
g.post_time = 0.5;
g.monosyn_window = 0.05;  % 0-25ms
g.smoothvalue = 7;
g.plotwin = [0.05 0.05];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.xlinewidth = 2;
g.clim_percentile = 95;

g.onset_threshold = 5;  % Lower threshold for monosynaptic detection
g.min_consec_bins = max(1, round(0.001 / g.bin_time));  % 1ms for brief monosynaptic responses
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

%% Calculate PSTHs for both conditions
fprintf('Calculating PSTHs...\n');
psthZ_full = cell(1,2);
psthHz_full = cell(1,2);
baseline_bins = round(g.pre_time / g.bin_time);

for hmp = 1:2
    psth_spx_og = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{hmp}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Convert to Hz: divide by number of trials and bin time
    num_trials = size(cell_metrics.general.(ttl{hmp}){1}, 1);
    psth_hz = psth_spx_og / (num_trials * g.bin_time);
    psth_hz_smooth = smoothdata(psth_hz, 2, 'sgolay', g.smoothvalue);
    psthHz_full{hmp} = psth_hz_smooth;

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

    % For LA, separate PNs and INs
    if strcmp(brain_regions{br}, 'LA')
        idx_PN = idx_neurons & strcmp(cell_metrics.putativeCellType, 'PN');
        idx_IN = idx_neurons & strcmp(cell_metrics.putativeCellType, 'IN');
    else
        idx_PN = idx_neurons;  % For Astria, treat all as one group
        idx_IN = false(size(idx_neurons));
    end

    n_neurons = sum(idx_neurons);

    if n_neurons == 0
        continue;
    end

    % Get number of trials for each condition
    num_trials_nolight = size(cell_metrics.general.(ttl{1}){1}, 1);
    num_trials_light = size(cell_metrics.general.(ttl{2}){1}, 1);

    % Extract monosynaptic window responses for both conditions
    % Separate PNs and INs
    psth_nolight_PN = psthZ_full{1}(idx_PN, :);
    psth_light_PN = psthZ_full{2}(idx_PN, :);

    psth_nolight_Hz_PN = psthHz_full{1}(idx_PN, :);
    psth_light_Hz_PN = psthHz_full{2}(idx_PN, :);

    psth_nolight_IN = psthZ_full{1}(idx_IN, :);
    psth_light_IN = psthZ_full{2}(idx_IN, :);

    psth_nolight_Hz_IN = psthHz_full{1}(idx_IN, :);
    psth_light_Hz_IN = psthHz_full{2}(idx_IN, :);

    n_PN = sum(idx_PN);
    n_IN = sum(idx_IN);

    % Calculate spike probabilities trial-by-trial for both conditions
    [~, ~, postAP_norm_nolight] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{1}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    [~, ~, postAP_norm_light] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{2}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Process PNs
    fprintf('  Processing PNs (%d neurons)...\n', n_PN);
    peak_nolight_PN = max(psth_nolight_PN(:, monosyn_window_bins), [], 2);
    peak_light_PN = max(psth_light_PN(:, monosyn_window_bins), [], 2);

    neuron_indices_PN = find(idx_PN);
    prob_nolight_PN = zeros(n_PN, 1);
    prob_light_PN = zeros(n_PN, 1);

    for n = 1:n_PN
        global_idx = neuron_indices_PN(n);

        % No-light probability
        if ~isempty(postAP_norm_nolight{global_idx})
            responsive_trials = 0;
            for trial = 1:length(postAP_norm_nolight{global_idx})
                trial_spikes = postAP_norm_nolight{global_idx}{trial};
                if any(trial_spikes > 0 & trial_spikes <= g.monosyn_window)
                    responsive_trials = responsive_trials + 1;
                end
            end
            prob_nolight_PN(n) = responsive_trials / num_trials_nolight;
        end

        % Light probability
        if ~isempty(postAP_norm_light{global_idx})
            responsive_trials = 0;
            for trial = 1:length(postAP_norm_light{global_idx})
                trial_spikes = postAP_norm_light{global_idx}{trial};
                if any(trial_spikes > 0 & trial_spikes <= g.monosyn_window)
                    responsive_trials = responsive_trials + 1;
                end
            end
            prob_light_PN(n) = responsive_trials / num_trials_light;
        end
    end

    % Responsiveness for PNs in EITHER condition
    if g.use_two_rule
        % Two-rule: (Rule 1 OR Rule 2)
        responsive_nolight_PN = (peak_nolight_PN >= g.zscore_threshold_rule1 & prob_nolight_PN >= g.prob_threshold_rule1) | ...
                                (peak_nolight_PN >= g.zscore_threshold_rule2 & prob_nolight_PN >= g.prob_threshold_rule2);
        responsive_light_PN = (peak_light_PN >= g.zscore_threshold_rule1 & prob_light_PN >= g.prob_threshold_rule1) | ...
                              (peak_light_PN >= g.zscore_threshold_rule2 & prob_light_PN >= g.prob_threshold_rule2);
    else
        % One-rule: z-score only
        responsive_nolight_PN = peak_nolight_PN >= g.zscore_threshold_one_rule;
        responsive_light_PN = peak_light_PN >= g.zscore_threshold_one_rule;
    end

    % Combine: responsive in EITHER condition
    responsive_either_PN = responsive_nolight_PN | responsive_light_PN;

    fprintf('  PNs - No-light responsive: %d, Light responsive: %d, Either: %d\n', ...
        sum(responsive_nolight_PN), sum(responsive_light_PN), sum(responsive_either_PN));

    % Process INs
    fprintf('  Processing INs (%d neurons)...\n', n_IN);
    peak_nolight_IN = max(psth_nolight_IN(:, monosyn_window_bins), [], 2);
    peak_light_IN = max(psth_light_IN(:, monosyn_window_bins), [], 2);

    neuron_indices_IN = find(idx_IN);
    prob_nolight_IN = zeros(n_IN, 1);
    prob_light_IN = zeros(n_IN, 1);

    for n = 1:n_IN
        global_idx = neuron_indices_IN(n);

        % No-light probability
        if ~isempty(postAP_norm_nolight{global_idx})
            responsive_trials = 0;
            for trial = 1:length(postAP_norm_nolight{global_idx})
                trial_spikes = postAP_norm_nolight{global_idx}{trial};
                if any(trial_spikes > 0 & trial_spikes <= g.monosyn_window)
                    responsive_trials = responsive_trials + 1;
                end
            end
            prob_nolight_IN(n) = responsive_trials / num_trials_nolight;
        end

        % Light probability
        if ~isempty(postAP_norm_light{global_idx})
            responsive_trials = 0;
            for trial = 1:length(postAP_norm_light{global_idx})
                trial_spikes = postAP_norm_light{global_idx}{trial};
                if any(trial_spikes > 0 & trial_spikes <= g.monosyn_window)
                    responsive_trials = responsive_trials + 1;
                end
            end
            prob_light_IN(n) = responsive_trials / num_trials_light;
        end
    end

    % Responsiveness for INs in EITHER condition
    if g.use_two_rule
        % Two-rule: (Rule 1 OR Rule 2)
        responsive_nolight_IN = (peak_nolight_IN >= g.zscore_threshold_rule1 & prob_nolight_IN >= g.prob_threshold_rule1) | ...
                                (peak_nolight_IN >= g.zscore_threshold_rule2 & prob_nolight_IN >= g.prob_threshold_rule2);
        responsive_light_IN = (peak_light_IN >= g.zscore_threshold_rule1 & prob_light_IN >= g.prob_threshold_rule1) | ...
                              (peak_light_IN >= g.zscore_threshold_rule2 & prob_light_IN >= g.prob_threshold_rule2);
    else
        % One-rule: z-score only
        responsive_nolight_IN = peak_nolight_IN >= g.zscore_threshold_one_rule;
        responsive_light_IN = peak_light_IN >= g.zscore_threshold_one_rule;
    end

    % Combine: responsive in EITHER condition
    responsive_either_IN = responsive_nolight_IN | responsive_light_IN;

    fprintf('  INs - No-light responsive: %d, Light responsive: %d, Either: %d\n', ...
        sum(responsive_nolight_IN), sum(responsive_light_IN), sum(responsive_either_IN));

    % Store results
    results_all{br}.responsive_nolight_PN = responsive_nolight_PN;
    results_all{br}.responsive_light_PN = responsive_light_PN;
    results_all{br}.responsive_either_PN = responsive_either_PN;
    results_all{br}.responsive_nolight_IN = responsive_nolight_IN;
    results_all{br}.responsive_light_IN = responsive_light_IN;
    results_all{br}.responsive_either_IN = responsive_either_IN;

    results_all{br}.psth_nolight_PN = psth_nolight_PN;
    results_all{br}.psth_light_PN = psth_light_PN;
    results_all{br}.psth_nolight_Hz_PN = psth_nolight_Hz_PN;
    results_all{br}.psth_light_Hz_PN = psth_light_Hz_PN;

    results_all{br}.psth_nolight_IN = psth_nolight_IN;
    results_all{br}.psth_light_IN = psth_light_IN;
    results_all{br}.psth_nolight_Hz_IN = psth_nolight_Hz_IN;
    results_all{br}.psth_light_Hz_IN = psth_light_Hz_IN;

    results_all{br}.peak_nolight_PN = peak_nolight_PN;
    results_all{br}.peak_light_PN = peak_light_PN;
    results_all{br}.peak_nolight_IN = peak_nolight_IN;
    results_all{br}.peak_light_IN = peak_light_IN;

    results_all{br}.prob_nolight_PN = prob_nolight_PN;
    results_all{br}.prob_light_PN = prob_light_PN;
    results_all{br}.prob_nolight_IN = prob_nolight_IN;
    results_all{br}.prob_light_IN = prob_light_IN;

    results_all{br}.n_PN = n_PN;
    results_all{br}.n_IN = n_IN;
    results_all{br}.n_neurons = n_neurons;
end

fprintf('\nMonosynaptic detection complete.\n');
fprintf('\nLA: %d PNs, %d INs\n', results_all{1}.n_PN, results_all{1}.n_IN);
fprintf('Astria: %d neurons\n', results_all{2}.n_PN);

%% Prepare data for BAfc_monosyn_raster_ui
fprintf('\nPreparing data for raster UI...\n');

% Create g structure for UI
g_ui.cell_metrics = cell_metrics;
g_ui.params.monosyn_window = g.monosyn_window;
g_ui.params.pre_time_short = 0.1;  % 100ms for display
g_ui.params.post_time_short = 0.1;  % 100ms for display

% Two-rule parameters
if g.use_two_rule
    g_ui.params.zscore_threshold = g.zscore_threshold_rule1;
    g_ui.params.prob_threshold = g.prob_threshold_rule1;
    g_ui.params.zscore_threshold_strict = g.zscore_threshold_rule2;
    g_ui.params.prob_threshold_lenient = g.prob_threshold_rule2;
else
    % One-rule: use same threshold for both
    g_ui.params.zscore_threshold = g.zscore_threshold_one_rule;
    g_ui.params.prob_threshold = 0;  % No probability threshold in one-rule mode
    g_ui.params.zscore_threshold_strict = g.zscore_threshold_one_rule;
    g_ui.params.prob_threshold_lenient = 0;
end

% Create monosyn_results structure for each TTL
monosyn_results = struct();

for tt = 1:length(ttl)
    ttl_fn = matlab.lang.makeValidName(ttl{tt});

    % Collect data from all brain regions
    all_neuron_idx = [];
    all_peak_zscore = [];
    all_response_probability = [];
    all_mean_latency = [];
    responsive_neuron_idx = [];

    for br = 1:2
        if isempty(results_all{br})
            continue;
        end

        res = results_all{br};

        % Get neuron indices for this region
        idx_neurons = strcmp(cell_metrics.brainRegion, brain_regions{br});
        if strcmp(brain_regions{br}, 'LA')
            idx_PN = idx_neurons & strcmp(cell_metrics.putativeCellType, 'PN');
        else
            idx_PN = idx_neurons;
        end
        neuron_indices = find(idx_PN);

        % Get peak z-scores and probabilities based on condition
        if tt == 1  % No-light
            peak_zscore = res.peak_nolight_PN;
            prob = res.prob_nolight_PN;
            responsive = res.responsive_either_PN;  % Use EITHER for UI
        else  % Light
            peak_zscore = res.peak_light_PN;
            prob = res.prob_light_PN;
            responsive = res.responsive_either_PN;  % Use EITHER for UI
        end

        % Mean latencies (placeholder)
        mean_latency = nan(length(neuron_indices), 1);

        % Append to arrays (ensure column vectors)
        all_neuron_idx = [all_neuron_idx; neuron_indices(:)];
        all_peak_zscore = [all_peak_zscore; peak_zscore(:)];
        all_response_probability = [all_response_probability; prob(:)];
        all_mean_latency = [all_mean_latency; mean_latency(:)];
        responsive_idx = neuron_indices(responsive);
        responsive_neuron_idx = [responsive_neuron_idx; responsive_idx(:)];
    end

    % Store in monosyn_results structure
    monosyn_results.(ttl_fn).all_neuron_idx = all_neuron_idx;
    monosyn_results.(ttl_fn).all_peak_zscore = all_peak_zscore;
    monosyn_results.(ttl_fn).all_response_probability = all_response_probability;
    monosyn_results.(ttl_fn).all_mean_latency = all_mean_latency;
    monosyn_results.(ttl_fn).neuron_idx = responsive_neuron_idx;

    fprintf('  %s: %d neurons, %d responsive\n', ttl{tt}, length(all_neuron_idx), length(responsive_neuron_idx));
end

fprintf('\nData prepared. Launch UI with:\n');
fprintf('  BAfc_monosyn_raster_ui(g_ui, monosyn_results, ttl)\n');

%% Statistical comparison: No-light vs Light conditions
fprintf('\n=== Statistical Comparison: Light vs No-Light (CS+US) ===\n');

% Response window for comparison (12-25ms, excluding artifact)
response_start_bin = round((g.pre_time + 0.012) / g.bin_time);
response_end_bin = round((g.pre_time + g.monosyn_window) / g.bin_time);
response_window_sec = [0.012, g.monosyn_window];

comparison_results = struct();

for br = 1:2
    region = brain_regions{br};

    if isempty(results_all{br})
        comparison_results.(region).n_increased = 0;
        comparison_results.(region).n_decreased = 0;
        comparison_results.(region).n_unchanged = 0;
        comparison_results.(region).n_total = 0;
        continue;
    end

    res = results_all{br};

    % Use neurons responsive in EITHER condition
    responsive = res.responsive_either_PN;
    n_responsive = sum(responsive);

    if n_responsive == 0
        comparison_results.(region).n_increased = 0;
        comparison_results.(region).n_decreased = 0;
        comparison_results.(region).n_unchanged = 0;
        comparison_results.(region).n_total = 0;
        continue;
    end

    % Get neuron indices
    idx_neurons = strcmp(cell_metrics.brainRegion, region);
    if strcmp(region, 'LA')
        idx_PN = idx_neurons & strcmp(cell_metrics.putativeCellType, 'PN');
    else
        idx_PN = idx_neurons;
    end
    all_neuron_indices = find(idx_PN);
    responsive_neuron_indices = all_neuron_indices(responsive);

    % Get trial-by-trial spike data
    [~, ~, postAP_norm_nolight] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{1}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    [~, ~, postAP_norm_light] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{2}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Bootstrap test for each neuron
    n_increased = 0;
    n_decreased = 0;
    n_unchanged = 0;

    increased_idx = [];
    decreased_idx = [];
    unchanged_idx = [];

    % Store p-values for all windows for each neuron
    neuron_pvalues = cell(n_responsive, 1);  % Cell array to store p-values per neuron

    % Define time windows: 12-13ms, 12-14ms, 12-15ms, ..., 12-40ms (in 1ms steps)
    artifact_end = 0.012;  % 12ms artifact exclusion
    window_ends = 0.013:0.001:0.040;  % 13ms to 40ms in 1ms steps
    n_windows = length(window_ends);

    for n = 1:n_responsive
        global_idx = responsive_neuron_indices(n);

        % Test each window
        is_increased_any_window = false;
        is_decreased_any_window = false;

        % Store p-values for this neuron
        p_vals = nan(1, n_windows);

        for w = 1:n_windows
            test_window = [artifact_end, window_ends(w)];

            % Get spike counts per trial in this test window
            spikes_nolight_trials = [];
            if ~isempty(postAP_norm_nolight{global_idx})
                for trial = 1:length(postAP_norm_nolight{global_idx})
                    trial_spikes = postAP_norm_nolight{global_idx}{trial};
                    spike_count = sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                    spikes_nolight_trials = [spikes_nolight_trials; spike_count];
                end
            end

            spikes_light_trials = [];
            if ~isempty(postAP_norm_light{global_idx})
                for trial = 1:length(postAP_norm_light{global_idx})
                    trial_spikes = postAP_norm_light{global_idx}{trial};
                    spike_count = sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                    spikes_light_trials = [spikes_light_trials; spike_count];
                end
            end

            % Wilcoxon signed-rank test for this window
            try
                [p, ~, ~] = signrank(spikes_nolight_trials, spikes_light_trials);
                p_vals(w) = p;
                obs_diff = mean(spikes_light_trials) - mean(spikes_nolight_trials);

                if p < 0.05
                    if obs_diff > 0
                        is_increased_any_window = true;
                    else
                        is_decreased_any_window = true;
                    end
                end
            catch
                % Test failed, skip this window
                p_vals(w) = NaN;
            end
        end

        % Store p-values for this neuron
        neuron_pvalues{n} = p_vals;

        % Classify neuron based on ANY window showing significance
        if is_increased_any_window
            n_increased = n_increased + 1;
            increased_idx = [increased_idx; global_idx];
        elseif is_decreased_any_window
            n_decreased = n_decreased + 1;
            decreased_idx = [decreased_idx; global_idx];
        else
            n_unchanged = n_unchanged + 1;
            unchanged_idx = [unchanged_idx; global_idx];
        end
    end

    % Store results
    comparison_results.(region).n_increased = n_increased;
    comparison_results.(region).n_decreased = n_decreased;
    comparison_results.(region).n_unchanged = n_unchanged;
    comparison_results.(region).n_total = n_responsive;
    comparison_results.(region).increased_idx = increased_idx;
    comparison_results.(region).decreased_idx = decreased_idx;
    comparison_results.(region).unchanged_idx = unchanged_idx;
    comparison_results.(region).all_neuron_indices = responsive_neuron_indices;
    comparison_results.(region).neuron_pvalues = neuron_pvalues;

    fprintf('\n%s: %d responsive neurons (in either condition)\n', region, n_responsive);
    fprintf('  Increased (p<0.05): %d (%.1f%%)\n', n_increased, n_increased/n_responsive*100);
    fprintf('  Decreased (p<0.05): %d (%.1f%%)\n', n_decreased, n_decreased/n_responsive*100);
    fprintf('  Unchanged (p>=0.05): %d (%.1f%%)\n', n_unchanged, n_unchanged/n_responsive*100);
end

%% Visualization: Pie charts showing percentage changes
fig_comparison = figure('Units', 'pixels', 'Position', [100, 100, 800, 400], 'Visible', 'on');
t = tiledlayout(fig_comparison, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

region_names = {'LA', 'Astria'};

for r = 1:2
    region = region_names{r};

    % Get data
    n_increased = comparison_results.(region).n_increased;
    n_decreased = comparison_results.(region).n_decreased;
    n_unchanged = comparison_results.(region).n_unchanged;
    n_total = comparison_results.(region).n_total;

    % Pie chart
    ax_pie = nexttile(t, (r-1)*2 + 1);
    pie_data = [n_increased, n_decreased, n_unchanged];
    colors = [0.8 0.2 0.2; 0.2 0.4 0.8; 0.7 0.7 0.7];

    if sum(pie_data) > 0
        pie(ax_pie, pie_data);
        colormap(ax_pie, colors);
    end

    title(sprintf('%s CS+US (n=%d)', region, n_total), 'FontSize', g.fontSize1, 'FontWeight', 'bold');
    legend({'Increased', 'Decreased', 'Unchanged'}, 'Location', 'southoutside', 'FontSize', g.fontSize2);

    % Bar chart with percentages
    ax_bar = nexttile(t, (r-1)*2 + 2);
    pct_data = [n_increased, n_decreased, n_unchanged] / max(n_total, 1) * 100;
    bar_h = bar(ax_bar, 1:3, pct_data, 'FaceColor', 'flat');
    bar_h.CData = colors(1:3, :);

    xticks(ax_bar, 1:3);
    xticklabels(ax_bar, {'Inc', 'Dec', 'Unch'});
    ylabel(ax_bar, 'Percentage (%)', 'FontSize', g.fontSize2);
    ylim(ax_bar, [0 100]);
    title(sprintf('%.1f%% / %.1f%% / %.1f%%', pct_data), 'FontSize', g.fontSize2);
    set(ax_bar, 'FontSize', g.fontSize2);
    box(ax_bar, 'off');

    % Add percentage labels on bars
    hold(ax_bar, 'on');
    for i = 1:2
        if pct_data(i) > 5
            text(ax_bar, i, pct_data(i)/2, sprintf('%.0f%%', pct_data(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Color', 'w');
        end
    end
    hold(ax_bar, 'off');
end

fprintf('\nVisualization complete.\n');

%% Print detailed summary with cellIDs
fprintf('\n========================================\n');
fprintf('DETAILED SUMMARY BY CATEGORY (CS+US)\n');
fprintf('========================================\n');

for r = 1:2
    region = region_names{r};
    fprintf('\n--- %s ---\n', region);

    if comparison_results.(region).n_total == 0
        fprintf('  No responsive neurons\n');
        continue;
    end

    % Helper function to find p-values for a neuron
    all_neuron_indices = comparison_results.(region).all_neuron_indices;
    neuron_pvalues = comparison_results.(region).neuron_pvalues;

    % Print each category with p-values
    fprintf('\nINCREASED (p<0.05, n=%d):\n', comparison_results.(region).n_increased);
    fprintf('  Window labels: 12-13ms, 12-14ms, 12-15ms, ..., 12-40ms (28 windows in 1ms steps)\n');
    if ~isempty(comparison_results.(region).increased_idx)
        for i = 1:length(comparison_results.(region).increased_idx)
            idx = comparison_results.(region).increased_idx(i);
            cellID = cell_metrics.cellID(idx);
            animal = cell_metrics.animal{idx};

            % Find this neuron's p-values
            neuron_pos = find(all_neuron_indices == idx);
            if ~isempty(neuron_pos)
                p_vals = neuron_pvalues{neuron_pos};
                % Find minimum p-value and corresponding window
                [min_p, min_idx] = min(p_vals);
                window_end_ms = (0.013 + (min_idx-1)*0.001) * 1000;  % Convert to ms
                fprintf('  neuron_idx=%d, cellID=%d, animal=%s, min_p=%.4f (12-%.0fms)\n', ...
                    idx, cellID, animal, min_p, window_end_ms);
            else
                fprintf('  neuron_idx=%d, cellID=%d, animal=%s\n', idx, cellID, animal);
            end
        end
    end

    fprintf('\nDECREASED (p<0.05, n=%d):\n', comparison_results.(region).n_decreased);
    fprintf('  Window labels: 12-13ms, 12-14ms, 12-15ms, ..., 12-40ms (28 windows in 1ms steps)\n');
    if ~isempty(comparison_results.(region).decreased_idx)
        for i = 1:length(comparison_results.(region).decreased_idx)
            idx = comparison_results.(region).decreased_idx(i);
            cellID = cell_metrics.cellID(idx);
            animal = cell_metrics.animal{idx};

            % Find this neuron's p-values
            neuron_pos = find(all_neuron_indices == idx);
            if ~isempty(neuron_pos)
                p_vals = neuron_pvalues{neuron_pos};
                % Find minimum p-value and corresponding window
                [min_p, min_idx] = min(p_vals);
                window_end_ms = (0.013 + (min_idx-1)*0.001) * 1000;  % Convert to ms
                fprintf('  neuron_idx=%d, cellID=%d, animal=%s, min_p=%.4f (12-%.0fms)\n', ...
                    idx, cellID, animal, min_p, window_end_ms);
            else
                fprintf('  neuron_idx=%d, cellID=%d, animal=%s\n', idx, cellID, animal);
            end
        end
    end

    fprintf('\nUNCHANGED (p>=0.05, n=%d):\n', comparison_results.(region).n_unchanged);
    fprintf('  Window labels: 12-13ms, 12-14ms, 12-15ms, ..., 12-40ms (28 windows in 1ms steps)\n');
    if ~isempty(comparison_results.(region).unchanged_idx)
        for i = 1:length(comparison_results.(region).unchanged_idx)
            idx = comparison_results.(region).unchanged_idx(i);
            cellID = cell_metrics.cellID(idx);
            animal = cell_metrics.animal{idx};

            % Find this neuron's p-values
            neuron_pos = find(all_neuron_indices == idx);
            if ~isempty(neuron_pos)
                p_vals = neuron_pvalues{neuron_pos};
                % Find minimum p-value and corresponding window
                [min_p, min_idx] = min(p_vals);
                window_end_ms = (0.013 + (min_idx-1)*0.001) * 1000;  % Convert to ms
                fprintf('  neuron_idx=%d, cellID=%d, animal=%s, min_p=%.4f (12-%.0fms)\n', ...
                    idx, cellID, animal, min_p, window_end_ms);
            else
                fprintf('  neuron_idx=%d, cellID=%d, animal=%s\n', idx, cellID, animal);
            end
        end
    end
end

fprintf('\n========================================\n');
fprintf('To check specific neurons in the UI:\n');
fprintf('  BAfc_monosyn_raster_ui(g_ui, monosyn_results, ttl)\n');
fprintf('  Then use "Go to Neuron_idx" field to jump to specific neurons\n');
fprintf('========================================\n');
