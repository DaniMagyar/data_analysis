% BAfc_figure_5.m
% Optogenetic manipulation of monosynaptic responses for LA and Astria
% Processes CS and US separately with their respective light conditions
% Compares: CS vs CS+light, US vs US+light
% To run gui: BAfc_monosyn_raster_ui(g_ui, monosyn_results, ttl)

clear all; close all

%% TESTING METHOD SELECTION
% Choose testing method for optogenetic enhancement detection:
% 'multi_window' - Test multiple 1ms windows from artifact_end to monosyn_window
%                  Neuron enhanced if ANY window p<0.05 with positive change
%                  WARNING: High false positive rate (~86% for 38 windows)
% 'single_window' - Test single pre-defined window
%                   No multiple comparisons issue, statistically sound
g.testing_method = 'single_window';  % Options: 'multi_window' or 'single_window'

% Single window parameters (only used if testing_method = 'single_window')
g.test_window_start = 0.012;  % Start of test window (s) - after artifact
g.test_window_end = 0.050;    % End of test window (s) - monosynaptic window
% Common options:
%   Full monosynaptic: [0.012, 0.050]
%   Early responses:   [0.012, 0.020]
%   Sustained:         [0.020, 0.050]

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

% TTL types: CS, CS+light, US, US+light
ttl = {'triptest_sound_only', 'triptest_sound_only_light', ...
       'triptest_shocks_only', 'triptest_shocks_only_light'};
brain_regions = {'LA', 'Astria'};

cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
cell_metrics = BAfc_match_celltypes('cell_metrics', cell_metrics);

g.colors = BAfc_colors;
g.fontSize1 = 10;
g.fontSize2 = 10;
g.bin_time = 0.001;
g.pre_time = 5;
g.post_time = 0.5;
g.monosyn_window = 0.05;  % 0-50ms
g.smoothvalue = 5;
g.plotwin = [0.05 0.05];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.xlinewidth = 2;
g.clim_percentile = 95;

g.onset_threshold = 5;
g.min_consec_bins = max(1, round(0.001 / g.bin_time));
g.alpha = 0.0;

% Responsiveness detection method
g.use_two_rule = true;

% Two-rule responsiveness parameters
g.zscore_threshold_rule1 = 3;
g.prob_threshold_rule1 = 0.25;
g.zscore_threshold_rule2 = 10;
g.prob_threshold_rule2 = 0.1;

% One-rule responsiveness parameter
g.zscore_threshold_one_rule = 5;

%% Calculate PSTHs for all 4 conditions
fprintf('Calculating PSTHs...\n');
psthZ_full = cell(1,4);
psthHz_full = cell(1,4);
baseline_bins = round(g.pre_time / g.bin_time);

% Calculate Savitzky-Golay filter delay correction
filter_delay = floor(g.smoothvalue / 2);  % Half the filter width (symmetric kernel)
fprintf('Savitzky-Golay filter width: %d bins, delay correction: %d bins (%.1f ms)\n', ...
    g.smoothvalue, filter_delay, filter_delay * g.bin_time * 1000);

for hmp = 1:4
    psth_spx_og = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{hmp}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Convert to Hz
    num_trials = size(cell_metrics.general.(ttl{hmp}){1}, 1);
    psth_hz = psth_spx_og / (num_trials * g.bin_time);
    psth_hz_smooth = smoothdata(psth_hz, 2, 'sgolay', g.smoothvalue);

    % Correct for filter delay by shifting forward (move data earlier in time)
    psth_hz_corrected = zeros(size(psth_hz_smooth));
    psth_hz_corrected(:, filter_delay+1:end) = psth_hz_smooth(:, 1:end-filter_delay);
    psth_hz_corrected(:, 1:filter_delay) = repmat(psth_hz_smooth(:, 1), 1, filter_delay);
    psthHz_full{hmp} = psth_hz_corrected;

    % Z-score using baseline period only
    baseline_mean = mean(psth_spx_og(:, 1:baseline_bins), 2);
    baseline_std = std(psth_spx_og(:, 1:baseline_bins), 0, 2);
    baseline_std(baseline_std == 0) = 1;
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;
    psth_spx_smooth = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);

    % Correct for filter delay by shifting forward (move data earlier in time)
    psth_spx_corrected = zeros(size(psth_spx_smooth));
    psth_spx_corrected(:, filter_delay+1:end) = psth_spx_smooth(:, 1:end-filter_delay);
    psth_spx_corrected(:, 1:filter_delay) = repmat(psth_spx_smooth(:, 1), 1, filter_delay);
    psthZ_full{hmp} = psth_spx_corrected;
end

%% Monosynaptic detection for each brain region and stimulus type
results_all = cell(2, 2);  % [region, stimulus_type] where stimulus_type: 1=CS, 2=US
monosyn_window_bins = round((g.pre_time)/g.bin_time+1 : (g.pre_time+g.monosyn_window)/g.bin_time);

for br = 1:2
    fprintf('\nProcessing %s...\n', brain_regions{br});
    idx_neurons = strcmp(cell_metrics.brainRegion, brain_regions{br});

    % For LA, separate PNs and INs
    if strcmp(brain_regions{br}, 'LA')
        idx_PN = idx_neurons & strcmp(cell_metrics.putativeCellType, 'PN');
        idx_IN = idx_neurons & strcmp(cell_metrics.putativeCellType, 'IN');
    else
        idx_PN = idx_neurons;
        idx_IN = false(size(idx_neurons));
    end

    n_neurons = sum(idx_neurons);
    n_PN = sum(idx_PN);
    n_IN = sum(idx_IN);

    if n_neurons == 0
        continue;
    end

    % Process CS and US separately
    for stim = 1:2  % 1=CS, 2=US
        stim_name = {'CS', 'US'};
        fprintf('  Processing %s...\n', stim_name{stim});

        % TTL indices: CS=1,2; US=3,4
        ttl_nolight_idx = stim*2 - 1;  % 1 for CS, 3 for US
        ttl_light_idx = stim*2;        % 2 for CS, 4 for US

        % Get number of trials
        num_trials_nolight = size(cell_metrics.general.(ttl{ttl_nolight_idx}){1}, 1);
        num_trials_light = size(cell_metrics.general.(ttl{ttl_light_idx}){1}, 1);

        % Extract PSTHs for PNs
        psth_nolight_PN = psthZ_full{ttl_nolight_idx}(idx_PN, :);
        psth_light_PN = psthZ_full{ttl_light_idx}(idx_PN, :);
        psth_nolight_Hz_PN = psthHz_full{ttl_nolight_idx}(idx_PN, :);
        psth_light_Hz_PN = psthHz_full{ttl_light_idx}(idx_PN, :);

        % Extract PSTHs for INs
        psth_nolight_IN = psthZ_full{ttl_nolight_idx}(idx_IN, :);
        psth_light_IN = psthZ_full{ttl_light_idx}(idx_IN, :);
        psth_nolight_Hz_IN = psthHz_full{ttl_nolight_idx}(idx_IN, :);
        psth_light_Hz_IN = psthHz_full{ttl_light_idx}(idx_IN, :);

        % Get trial-by-trial spike data
        [~, ~, postAP_norm_nolight] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{ttl_nolight_idx}, ...
            'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
        [~, ~, postAP_norm_light] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{ttl_light_idx}, ...
            'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

        % Process PNs
        fprintf('    Processing PNs (%d neurons)...\n', n_PN);
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

        % Responsiveness for PNs
        if g.use_two_rule
            responsive_nolight_PN = (peak_nolight_PN >= g.zscore_threshold_rule1 & prob_nolight_PN >= g.prob_threshold_rule1) | ...
                                    (peak_nolight_PN >= g.zscore_threshold_rule2 & prob_nolight_PN >= g.prob_threshold_rule2);
            responsive_light_PN = (peak_light_PN >= g.zscore_threshold_rule1 & prob_light_PN >= g.prob_threshold_rule1) | ...
                                  (peak_light_PN >= g.zscore_threshold_rule2 & prob_light_PN >= g.prob_threshold_rule2);
        else
            responsive_nolight_PN = peak_nolight_PN >= g.zscore_threshold_one_rule;
            responsive_light_PN = peak_light_PN >= g.zscore_threshold_one_rule;
        end

        responsive_either_PN = responsive_nolight_PN | responsive_light_PN;

        fprintf('    PNs - No-light: %d, Light: %d, Either: %d\n', ...
            sum(responsive_nolight_PN), sum(responsive_light_PN), sum(responsive_either_PN));

        % Process INs
        fprintf('    Processing INs (%d neurons)...\n', n_IN);
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

        % Responsiveness for INs
        if g.use_two_rule
            responsive_nolight_IN = (peak_nolight_IN >= g.zscore_threshold_rule1 & prob_nolight_IN >= g.prob_threshold_rule1) | ...
                                    (peak_nolight_IN >= g.zscore_threshold_rule2 & prob_nolight_IN >= g.prob_threshold_rule2);
            responsive_light_IN = (peak_light_IN >= g.zscore_threshold_rule1 & prob_light_IN >= g.prob_threshold_rule1) | ...
                                  (peak_light_IN >= g.zscore_threshold_rule2 & prob_light_IN >= g.prob_threshold_rule2);
        else
            responsive_nolight_IN = peak_nolight_IN >= g.zscore_threshold_one_rule;
            responsive_light_IN = peak_light_IN >= g.zscore_threshold_one_rule;
        end

        responsive_either_IN = responsive_nolight_IN | responsive_light_IN;

        fprintf('    INs - No-light: %d, Light: %d, Either: %d\n', ...
            sum(responsive_nolight_IN), sum(responsive_light_IN), sum(responsive_either_IN));

        % Store results
        results_all{br, stim}.responsive_nolight_PN = responsive_nolight_PN;
        results_all{br, stim}.responsive_light_PN = responsive_light_PN;
        results_all{br, stim}.responsive_either_PN = responsive_either_PN;
        results_all{br, stim}.responsive_nolight_IN = responsive_nolight_IN;
        results_all{br, stim}.responsive_light_IN = responsive_light_IN;
        results_all{br, stim}.responsive_either_IN = responsive_either_IN;

        results_all{br, stim}.psth_nolight_PN = psth_nolight_PN;
        results_all{br, stim}.psth_light_PN = psth_light_PN;
        results_all{br, stim}.psth_nolight_Hz_PN = psth_nolight_Hz_PN;
        results_all{br, stim}.psth_light_Hz_PN = psth_light_Hz_PN;

        results_all{br, stim}.psth_nolight_IN = psth_nolight_IN;
        results_all{br, stim}.psth_light_IN = psth_light_IN;
        results_all{br, stim}.psth_nolight_Hz_IN = psth_nolight_Hz_IN;
        results_all{br, stim}.psth_light_Hz_IN = psth_light_Hz_IN;

        results_all{br, stim}.peak_nolight_PN = peak_nolight_PN;
        results_all{br, stim}.peak_light_PN = peak_light_PN;
        results_all{br, stim}.peak_nolight_IN = peak_nolight_IN;
        results_all{br, stim}.peak_light_IN = peak_light_IN;

        results_all{br, stim}.prob_nolight_PN = prob_nolight_PN;
        results_all{br, stim}.prob_light_PN = prob_light_PN;
        results_all{br, stim}.prob_nolight_IN = prob_nolight_IN;
        results_all{br, stim}.prob_light_IN = prob_light_IN;

        results_all{br, stim}.n_PN = n_PN;
        results_all{br, stim}.n_IN = n_IN;
        results_all{br, stim}.n_neurons = n_neurons;

        results_all{br, stim}.postAP_norm_nolight = postAP_norm_nolight;
        results_all{br, stim}.postAP_norm_light = postAP_norm_light;
    end
end

fprintf('\nMonosynaptic detection complete.\n');

%% Prepare data for BAfc_monosyn_raster_ui
fprintf('\nPreparing data for raster UI...\n');

g_ui.cell_metrics = cell_metrics;
g_ui.params.monosyn_window = g.monosyn_window;
g_ui.params.pre_time_short = 0.1;
g_ui.params.post_time_short = 0.1;

if g.use_two_rule
    g_ui.params.zscore_threshold = g.zscore_threshold_rule1;
    g_ui.params.prob_threshold = g.prob_threshold_rule1;
    g_ui.params.zscore_threshold_strict = g.zscore_threshold_rule2;
    g_ui.params.prob_threshold_lenient = g.prob_threshold_rule2;
else
    g_ui.params.zscore_threshold = g.zscore_threshold_one_rule;
    g_ui.params.prob_threshold = 0;
    g_ui.params.zscore_threshold_strict = g.zscore_threshold_one_rule;
    g_ui.params.prob_threshold_lenient = 0;
end

% Create monosyn_results structure for all 4 TTLs
monosyn_results = struct();

for tt = 1:length(ttl)
    ttl_fn = matlab.lang.makeValidName(ttl{tt});
    stim_idx = ceil(tt/2);  % 1 or 2 (CS or US)

    all_neuron_idx = [];
    all_peak_zscore = [];
    all_response_probability = [];
    all_mean_latency = [];
    responsive_neuron_idx = [];

    for br = 1:2
        if isempty(results_all{br, stim_idx})
            continue;
        end

        res = results_all{br, stim_idx};

        % Get neuron indices
        idx_neurons = strcmp(cell_metrics.brainRegion, brain_regions{br});
        if strcmp(brain_regions{br}, 'LA')
            idx_PN = idx_neurons & strcmp(cell_metrics.putativeCellType, 'PN');
        else
            idx_PN = idx_neurons;
        end
        neuron_indices = find(idx_PN);

        % Select data based on light condition (odd=nolight, even=light)
        if mod(tt, 2) == 1  % No-light
            peak_zscore = res.peak_nolight_PN;
            prob = res.prob_nolight_PN;
        else  % Light
            peak_zscore = res.peak_light_PN;
            prob = res.prob_light_PN;
        end
        responsive = res.responsive_either_PN;

        mean_latency = nan(length(neuron_indices), 1);

        all_neuron_idx = [all_neuron_idx; neuron_indices(:)];
        all_peak_zscore = [all_peak_zscore; peak_zscore(:)];
        all_response_probability = [all_response_probability; prob(:)];
        all_mean_latency = [all_mean_latency; mean_latency(:)];
        responsive_idx = neuron_indices(responsive);
        responsive_neuron_idx = [responsive_neuron_idx; responsive_idx(:)];
    end

    monosyn_results.(ttl_fn).all_neuron_idx = all_neuron_idx;
    monosyn_results.(ttl_fn).all_peak_zscore = all_peak_zscore;
    monosyn_results.(ttl_fn).all_response_probability = all_response_probability;
    monosyn_results.(ttl_fn).all_mean_latency = all_mean_latency;
    monosyn_results.(ttl_fn).neuron_idx = responsive_neuron_idx;

    fprintf('  %s: %d neurons, %d responsive\n', ttl{tt}, length(all_neuron_idx), length(responsive_neuron_idx));
end

fprintf('\nData prepared. Launch UI with:\n');
fprintf('  BAfc_monosyn_raster_ui(g_ui, monosyn_results, ttl)\n');

%% Statistical comparison for CS and US separately
fprintf('\n=== Statistical Comparison: Light vs No-Light ===\n');

% Setup testing windows based on method
artifact_end = 0.012;  % 12ms artifact exclusion

if strcmp(g.testing_method, 'multi_window')
    % Fine-grained multi-window testing
    window_ends = 0.013:0.001:g.monosyn_window;  % 13ms to monosyn_window in 1ms steps
    n_windows = length(window_ends);
    fprintf('Using MULTI-WINDOW testing: %d windows from %.0f-%.0f ms\n', ...
        n_windows, artifact_end*1000, g.monosyn_window*1000);
    fprintf('WARNING: Family-wise error rate ~%.1f%%\n', (1 - 0.95^n_windows)*100);
elseif strcmp(g.testing_method, 'single_window')
    % Single pre-defined window testing
    window_ends = g.test_window_end;  % Single window endpoint
    n_windows = 1;
    fprintf('Using SINGLE-WINDOW testing: %.0f-%.0f ms\n', ...
        g.test_window_start*1000, g.test_window_end*1000);
    fprintf('No multiple comparisons correction needed.\n');
    % Override artifact_end for single window method
    artifact_end = g.test_window_start;
else
    error('Invalid testing_method. Must be ''multi_window'' or ''single_window''');
end

comparison_results = struct();

for br = 1:2
    region = brain_regions{br};

    for stim = 1:2
        stim_name = {'CS', 'US'};
        result_field = sprintf('%s_%s', region, stim_name{stim});

        fprintf('\n%s - %s:\n', region, stim_name{stim});

        if isempty(results_all{br, stim})
            comparison_results.(result_field).n_increased = 0;
            comparison_results.(result_field).n_decreased = 0;
            comparison_results.(result_field).n_unchanged = 0;
            comparison_results.(result_field).n_total = 0;
            continue;
        end

        res = results_all{br, stim};
        responsive = res.responsive_either_PN;
        n_responsive = sum(responsive);

        if n_responsive == 0
            comparison_results.(result_field).n_increased = 0;
            comparison_results.(result_field).n_decreased = 0;
            comparison_results.(result_field).n_unchanged = 0;
            comparison_results.(result_field).n_total = 0;
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

        % Get trial-by-trial data
        postAP_norm_nolight = res.postAP_norm_nolight;
        postAP_norm_light = res.postAP_norm_light;

        n_increased = 0;
        n_decreased = 0;
        n_unchanged = 0;
        increased_idx = [];
        decreased_idx = [];
        unchanged_idx = [];
        neuron_pvalues = cell(n_responsive, 1);

        for n = 1:n_responsive
            global_idx = responsive_neuron_indices(n);

            is_increased_any_window = false;
            is_decreased_any_window = false;
            p_vals = nan(1, n_windows);

            for w = 1:n_windows
                test_window = [artifact_end, window_ends(w)];

                % Get spike counts
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

                % Wilcoxon test
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
                    p_vals(w) = NaN;
                end
            end

            neuron_pvalues{n} = p_vals;

            % Classify neuron
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
        comparison_results.(result_field).n_increased = n_increased;
        comparison_results.(result_field).n_decreased = n_decreased;
        comparison_results.(result_field).n_unchanged = n_unchanged;
        comparison_results.(result_field).n_total = n_responsive;
        comparison_results.(result_field).increased_idx = increased_idx;
        comparison_results.(result_field).decreased_idx = decreased_idx;
        comparison_results.(result_field).unchanged_idx = unchanged_idx;
        comparison_results.(result_field).all_neuron_indices = responsive_neuron_indices;
        comparison_results.(result_field).neuron_pvalues = neuron_pvalues;

        fprintf('  %d responsive neurons (either condition)\n', n_responsive);
        fprintf('  Increased: %d (%.1f%%)\n', n_increased, n_increased/n_responsive*100);
        fprintf('  Decreased: %d (%.1f%%)\n', n_decreased, n_decreased/n_responsive*100);
        fprintf('  Unchanged: %d (%.1f%%)\n', n_unchanged, n_unchanged/n_responsive*100);
    end
end

%% Visualization: Main figure with nested pie charts and slope graphs
fig_comparison = figure('Units', 'pixels', 'Position', [100, 100, 1000, 500], 'Visible', 'on');
t_main = tiledlayout(fig_comparison, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Create nested tiledlayout for raster plots (rows 1-2, columns 1-2)
t_nested_left = tiledlayout(t_main, 3, 4, 'TileSpacing', 'compact', 'Padding', 'tight');
t_nested_left.Layout.Tile = 1;  % Start at tile 1 (row 1, column 1)
t_nested_left.Layout.TileSpan = [2 2];  % Span 2 rows, 2 columns
title(t_nested_left, 'Example neurons', 'FontSize', 11, 'FontWeight', 'bold');

% Create nested tiledlayout for pie charts and slope graphs (rows 1-2, columns 3-4)
t_nested_right = tiledlayout(t_main, 2, 4, 'TileSpacing', 'tight', 'Padding', 'compact');
t_nested_right.Layout.Tile = 3;  % Start at tile 3 (row 1, column 3)
t_nested_right.Layout.TileSpan = [2 2];  % Span 2 rows, 2 columns

region_names = {'LA', 'Astria'};
stim_names = {'CS', 'US'};
spaghetti_titles = {'LA CS', 'LA US', 'AStria CS', 'AStria US'};

for r = 1:2
    for s = 1:2
        result_field = sprintf('%s_%s', region_names{r}, stim_names{s});

        % Get data (store separately but combine for visualization)
        n_increased = comparison_results.(result_field).n_increased;
        n_decreased = comparison_results.(result_field).n_decreased;
        n_unchanged = comparison_results.(result_field).n_unchanged;
        n_total = comparison_results.(result_field).n_total;

        % Combine for visualization: Enhanced vs Non-enhanced
        n_enhanced = n_increased;
        n_non_enhanced = n_decreased + n_unchanged;

        % Pie chart (row 1 of nested layout)
        col_idx = (r-1)*2 + s;
        ax_pie = nexttile(t_nested_right, col_idx);
        pie_data = [n_enhanced, n_non_enhanced];
        colors = [0.8 0.2 0.2; 0.7 0.7 0.7];

        if sum(pie_data) > 0
            % Create custom labels with counts and percentages
            percentages = 100 * pie_data / sum(pie_data);
            labels = arrayfun(@(x, p) sprintf('%d (%.1f%%)', x, p), pie_data, percentages, 'UniformOutput', false);
            pie(ax_pie, pie_data, labels);
            colormap(ax_pie, colors);
        end

        % Add title using spaghetti_titles with manual positioning
        plot_idx = (r-1)*2 + s;
        title_handle = title(ax_pie, spaghetti_titles{plot_idx}, ...
            'FontSize', g.fontSize1, 'FontWeight', 'bold');
        % Move title up to avoid overlap with percentage labels
        title_handle.Position(2) = title_handle.Position(2) + 0.15;

        % Slope graph - Delta Peak FR for enhanced neurons (row 2 of nested layout)
        ax_slope_peak = nexttile(t_nested_right, col_idx + 4);

        if n_enhanced > 0
            % Get enhanced neuron indices
            increased_idx = comparison_results.(result_field).increased_idx;
            all_neuron_indices = comparison_results.(result_field).all_neuron_indices;
            neuron_pvalues = comparison_results.(result_field).neuron_pvalues;

            % Get region and stimulus indices
            region = region_names{r};
            res = results_all{r, s};

            % TTL indices for this stimulus
            ttl_nolight_idx = s*2 - 1;
            ttl_light_idx = s*2;

            % Get trial-by-trial data
            postAP_norm_nolight = res.postAP_norm_nolight;
            postAP_norm_light = res.postAP_norm_light;

            % Calculate mean spike counts for each enhanced neuron in their most significant window
            mean_spikes_nolight = zeros(n_enhanced, 1);
            mean_spikes_light = zeros(n_enhanced, 1);

            for i = 1:n_enhanced
                idx = increased_idx(i);

                % Define response window based on testing method
                if strcmp(g.testing_method, 'multi_window')
                    % Find this neuron's most significant window
                    neuron_pos = find(all_neuron_indices == idx);
                    if ~isempty(neuron_pos)
                        p_vals = neuron_pvalues{neuron_pos};
                        [~, min_window_idx] = min(p_vals);
                        window_end = window_ends(min_window_idx);
                    else
                        window_end = window_ends(1);  % Fallback
                    end
                    test_window = [0.012, window_end];  % Original artifact_end
                elseif strcmp(g.testing_method, 'single_window')
                    % Use the pre-defined window for all neurons
                    test_window = [g.test_window_start, g.test_window_end];
                end

                % Get spike counts for no-light trials
                spikes_nolight_trials = [];
                if ~isempty(postAP_norm_nolight{idx})
                    for trial = 1:length(postAP_norm_nolight{idx})
                        trial_spikes = postAP_norm_nolight{idx}{trial};
                        spike_count = sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                        spikes_nolight_trials = [spikes_nolight_trials; spike_count];
                    end
                end
                mean_spikes_nolight(i) = mean(spikes_nolight_trials);

                % Get spike counts for light trials
                spikes_light_trials = [];
                if ~isempty(postAP_norm_light{idx})
                    for trial = 1:length(postAP_norm_light{idx})
                        trial_spikes = postAP_norm_light{idx}{trial};
                        spike_count = sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                        spikes_light_trials = [spikes_light_trials; spike_count];
                    end
                end
                mean_spikes_light(i) = mean(spikes_light_trials);
            end

            % Create slope graph for spike counts
            hold(ax_slope_peak, 'on');

            % Plot mean trajectory with thicker line (underneath)
            mean_spikes_nolight_avg = mean(mean_spikes_nolight);
            mean_spikes_light_avg = mean(mean_spikes_light);
            plot(ax_slope_peak, [1 2], [mean_spikes_nolight_avg mean_spikes_light_avg], ...
                'k-', 'LineWidth', 3);

            % Plot individual neuron trajectories
            for i = 1:n_enhanced
                plot(ax_slope_peak, [1 2], [mean_spikes_nolight(i) mean_spikes_light(i)], ...
                    '-', 'Color', [0.8 0.2 0.2 0.3], 'LineWidth', 1);
            end

            % Add scatter points on top
            scatter(ax_slope_peak, ones(n_enhanced, 1), mean_spikes_nolight, 30, ...
                [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
            scatter(ax_slope_peak, 2*ones(n_enhanced, 1), mean_spikes_light, 30, ...
                [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.6);

            % Add mean markers
            scatter(ax_slope_peak, 1, mean_spikes_nolight_avg, 80, 'k', 'filled');
            scatter(ax_slope_peak, 2, mean_spikes_light_avg, 80, 'k', 'filled');

            % Population-level paired t-test on enhanced neurons
            [~, p_pop, ~, ~] = ttest(mean_spikes_nolight, mean_spikes_light);
            if p_pop < 0.001
                sig_str = '***';
            elseif p_pop < 0.01
                sig_str = '**';
            elseif p_pop < 0.05
                sig_str = '*';
            else
                sig_str = 'n.s.';
            end

            % Add significance line and text above the plot
            y_max = max([mean_spikes_nolight; mean_spikes_light]);
            y_pos = y_max + y_max * 0.15;
            plot(ax_slope_peak, [1 2], [y_pos y_pos], 'k-', 'LineWidth', 1.5);
            text(ax_slope_peak, 1.5, y_pos, sig_str, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);

            hold(ax_slope_peak, 'off');

            % Set common y-limits for all spaghetti plots
            ylim(ax_slope_peak, [0 2.4]);
            yticks(ax_slope_peak, [0 1.2 2.4]);

            xlim(ax_slope_peak, [0.5 2.5]);
            xticks(ax_slope_peak, [1 2]);
            xticklabels(ax_slope_peak, {'No light', 'Light'});
            if r == 1 && s == 1
                ylabel(ax_slope_peak, 'Mean spike count', 'FontSize', g.fontSize2);
            end

            % Add title using spaghetti_titles
            plot_idx = (r-1)*2 + s;
            title(ax_slope_peak, spaghetti_titles{plot_idx}, 'FontSize', g.fontSize2);
            set(ax_slope_peak, 'FontSize', g.fontSize2);
            box(ax_slope_peak, 'off');
        else
            % No enhanced neurons
            text(ax_slope_peak, 0.5, 0.5, 'No enhanced neurons', 'HorizontalAlignment', 'center', ...
                'FontSize', g.fontSize2);
            axis(ax_slope_peak, 'off');
        end

    end
end

%% Add example raster plots (rows 1-2 of left nested layout)
% Define example neurons: [animal, cellID, stimulus_type, stimulus_type_light]
examples = {
    'MD309_001', 20, 'triptest_sound_only', 'triptest_sound_only_light';   % CS-only
    'MD309_001', 20, 'triptest_shocks_only', 'triptest_shocks_only_light';  % US-only
    'MD318_001', 46, 'triptest_sound_only', 'triptest_sound_only_light';   % CS-only
    'MD317_001', 43, 'triptest_shocks_only', 'triptest_shocks_only_light'   % US-only
};

for ex = 1:4
    % Row 1: No-light condition
    ax_raster_nolight = nexttile(t_nested_left, ex);

    animal_name = examples{ex, 1};
    target_cellID = examples{ex, 2};
    ttl_type = examples{ex, 3};

    % Find neuron index
    animal_match = strcmp(cell_metrics.animal, animal_name);
    cellID_match = cell_metrics.cellID == target_cellID;
    neuron_idx = find(animal_match & cellID_match);

    if ~isempty(neuron_idx)
        % Get spike times using BAfc_psth_spx
        [~, preAP_norm, postAP_norm] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl_type, ...
            'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

        preAP_spikes = preAP_norm{neuron_idx};
        postAP_spikes = postAP_norm{neuron_idx};

        if ~isempty(preAP_spikes) || ~isempty(postAP_spikes)
            % Plot raster
            hold(ax_raster_nolight, 'on');
            n_trials = length(preAP_spikes);
            for trial = 1:n_trials
                % Plot pre-stimulus spikes (negative times)
                if ~isempty(preAP_spikes{trial})
                    valid_pre = preAP_spikes{trial}(preAP_spikes{trial} >= -0.05 & preAP_spikes{trial} <= 0);
                    plot(ax_raster_nolight, valid_pre*1000, trial*ones(size(valid_pre)), 'k.', 'MarkerSize', 6);
                end
                % Plot post-stimulus spikes (positive times)
                if ~isempty(postAP_spikes{trial})
                    valid_post = postAP_spikes{trial}(postAP_spikes{trial} >= 0 & postAP_spikes{trial} <= 0.05);
                    plot(ax_raster_nolight, valid_post*1000, trial*ones(size(valid_post)), 'k.', 'MarkerSize', 6);
                end
            end

            % Add stimulus line
            plot(ax_raster_nolight, [0 0], [0 n_trials+1], 'r-', 'LineWidth', 2);

            hold(ax_raster_nolight, 'off');

            xlim(ax_raster_nolight, [-50 50]);
            ylim(ax_raster_nolight, [0 n_trials+1]);
            set(ax_raster_nolight, 'XTickLabel', []);
            if ex == 1
                ylabel(ax_raster_nolight, 'Trial # (laser OFF)', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
            else
                set(ax_raster_nolight, 'YTickLabel', []);
            end

            % Title only for upper row
            raster_titles = {'LA CS', 'LA US', 'AStria CS', 'AStria US'};
            title(ax_raster_nolight, raster_titles{ex}, ...
                'FontSize', g.fontSize2);
            set(ax_raster_nolight, 'FontSize', g.fontSize2);
            box(ax_raster_nolight, 'off');
        else
            % No spike data
            text(ax_raster_nolight, 0.5, 0.5, 'No spike data', 'HorizontalAlignment', 'center');
            axis(ax_raster_nolight, 'off');
        end
    else
        % Neuron not found
        text(ax_raster_nolight, 0.5, 0.5, 'Neuron not found', 'HorizontalAlignment', 'center');
        axis(ax_raster_nolight, 'off');
    end

    % Row 2: Light condition
    ax_raster_light = nexttile(t_nested_left, ex + 4);

    ttl_type_light = examples{ex, 4};

    if ~isempty(neuron_idx)
        % Get spike times using BAfc_psth_spx for light condition
        [~, preAP_norm_light, postAP_norm_light] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl_type_light, ...
            'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

        preAP_spikes_light = preAP_norm_light{neuron_idx};
        postAP_spikes_light = postAP_norm_light{neuron_idx};

        if ~isempty(preAP_spikes_light) || ~isempty(postAP_spikes_light)
            % Plot raster
            hold(ax_raster_light, 'on');
            n_trials = length(preAP_spikes_light);
            for trial = 1:n_trials
                % Plot pre-stimulus spikes (negative times)
                if ~isempty(preAP_spikes_light{trial})
                    valid_pre = preAP_spikes_light{trial}(preAP_spikes_light{trial} >= -0.05 & preAP_spikes_light{trial} <= 0);
                    plot(ax_raster_light, valid_pre*1000, trial*ones(size(valid_pre)), 'k.', 'MarkerSize', 6);
                end
                % Plot post-stimulus spikes (positive times)
                if ~isempty(postAP_spikes_light{trial})
                    valid_post = postAP_spikes_light{trial}(postAP_spikes_light{trial} >= 0 & postAP_spikes_light{trial} <= 0.05);
                    plot(ax_raster_light, valid_post*1000, trial*ones(size(valid_post)), 'k.', 'MarkerSize', 6);
                end
            end

            % Add stimulus line
            plot(ax_raster_light, [0 0], [0 n_trials+1], 'r-', 'LineWidth', 2);

            hold(ax_raster_light, 'off');

            xlim(ax_raster_light, [-50 50]);
            ylim(ax_raster_light, [0 n_trials+1]);
            set(ax_raster_light, 'XTickLabel', []);
            if ex == 1
                ylabel(ax_raster_light, 'Trial # (laser ON)', 'Color', 'r', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
            else
                set(ax_raster_light, 'YTickLabel', []);
            end

            % No title for lower row (light condition)
            set(ax_raster_light, 'FontSize', g.fontSize2);
            box(ax_raster_light, 'off');
        else
            % No spike data
            text(ax_raster_light, 0.5, 0.5, 'No spike data', 'HorizontalAlignment', 'center');
            axis(ax_raster_light, 'off');
        end
    else
        % Neuron not found
        text(ax_raster_light, 0.5, 0.5, 'Neuron not found', 'HorizontalAlignment', 'center');
        axis(ax_raster_light, 'off');
    end

    % Row 3: Firing rate line plot (no-light vs light)
    ax_lineplot = nexttile(t_nested_left, ex + 8);

    if ~isempty(neuron_idx)
        % Find indices in ttl array for this stimulus type
        ttl_nolight_idx = find(strcmp(ttl, ttl_type));
        ttl_light_idx = find(strcmp(ttl, ttl_type_light));

        % Get PSTH data for both conditions
        if ~isempty(ttl_nolight_idx) && ~isempty(ttl_light_idx)
            psth_nolight = psthHz_full{ttl_nolight_idx};
            psth_light = psthHz_full{ttl_light_idx};
        else
            psth_nolight = [];
            psth_light = [];
        end

        if ~isempty(psth_nolight) && ~isempty(psth_light)
            % Extract this neuron's PSTH
            neuron_psth_nolight = psth_nolight(neuron_idx, :);
            neuron_psth_light = psth_light(neuron_idx, :);

            % Time axis in ms, centered at stimulus
            time_axis = ((-g.pre_time:g.bin_time:(g.post_time-g.bin_time)) * 1000);

            % Plot both conditions
            hold(ax_lineplot, 'on');
            plot(ax_lineplot, time_axis, neuron_psth_nolight, 'k-', 'LineWidth', 1.5);
            plot(ax_lineplot, time_axis, neuron_psth_light, 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);

            hold(ax_lineplot, 'off');

            xlim(ax_lineplot, [-50 50]);
            ylim(ax_lineplot, [-10 300]);

            % Add stimulus line (after setting ylim)
            hold(ax_lineplot, 'on');
            plot(ax_lineplot, [0 0], [-10 300], 'r-', 'LineWidth', 2);
            hold(ax_lineplot, 'off');
            xlabel(ax_lineplot, 'Time (ms)', 'FontSize', g.fontSize2);
            if ex == 1
                ylabel(ax_lineplot, 'FR (Hz)', 'FontSize', g.fontSize2);
            else
                set(ax_lineplot, 'YTickLabel', []);
            end

            % Add legend only for first plot
            if ex == 1
                leg = legend(ax_lineplot, {'No light', 'Light'}, 'Location', 'northeast', ...
                    'FontSize', g.fontSize2, 'Box', 'off');
                leg.ItemTokenSize = [10, 12];  % Shorter lines
                leg.Position(1) = leg.Position(1) + 0.06;  % Shift right
                leg.Position(2) = leg.Position(2) + 0.03;  % Shift up
            end

            set(ax_lineplot, 'FontSize', g.fontSize2);
            box(ax_lineplot, 'off');
        else
            % No PSTH data
            text(ax_lineplot, 0.5, 0.5, 'No PSTH data', 'HorizontalAlignment', 'center');
            axis(ax_lineplot, 'off');
        end
    else
        % Neuron not found
        text(ax_lineplot, 0.5, 0.5, 'Neuron not found', 'HorizontalAlignment', 'center');
        axis(ax_lineplot, 'off');
    end
end

% Add title for right section using annotation
annotation(fig_comparison, 'textbox', [0.55, 0.9, 0.4, 0.04], ...
    'String', 'Neurons with enhanced response', ...
    'FontSize', 11, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none');

%% Add panel labels
% A: Rasterplots (top-left of left panel)
annotation(fig_comparison, 'textbox', [0.01, 0.95, 0.03, 0.04], ...
    'String', 'A', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% B: FR comparison lineplots (bottom-left of left panel, row 3)
annotation(fig_comparison, 'textbox', [0.01, 0.35, 0.03, 0.04], ...
    'String', 'B', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% C: Pie charts (top-left of right panel)
annotation(fig_comparison, 'textbox', [0.51, 0.95, 0.03, 0.04], ...
    'String', 'C', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% D: Spaghetti plots (bottom-left of right panel)
annotation(fig_comparison, 'textbox', [0.51, 0.50, 0.03, 0.04], ...
    'String', 'D', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

fprintf('\nVisualization complete.\n');

%% Print detailed summary
fprintf('\n========================================\n');
fprintf('DETAILED SUMMARY BY CATEGORY\n');
fprintf('========================================\n');

for r = 1:2
    for s = 1:2
        result_field = sprintf('%s_%s', region_names{r}, stim_names{s});
        fprintf('\n--- %s %s ---\n', region_names{r}, stim_names{s});

        if comparison_results.(result_field).n_total == 0
            fprintf('  No responsive neurons\n');
            continue;
        end

        all_neuron_indices = comparison_results.(result_field).all_neuron_indices;
        neuron_pvalues = comparison_results.(result_field).neuron_pvalues;

        % Print increased
        fprintf('\nINCREASED (p<0.05, n=%d):\n', comparison_results.(result_field).n_increased);
        if ~isempty(comparison_results.(result_field).increased_idx)
            for i = 1:length(comparison_results.(result_field).increased_idx)
                idx = comparison_results.(result_field).increased_idx(i);
                cellID = cell_metrics.cellID(idx);
                animal = cell_metrics.animal{idx};

                neuron_pos = find(all_neuron_indices == idx);
                if ~isempty(neuron_pos)
                    p_vals = neuron_pvalues{neuron_pos};
                    if strcmp(g.testing_method, 'multi_window')
                        [min_p, min_idx] = min(p_vals);
                        window_end_ms = (0.013 + (min_idx-1)*0.001) * 1000;
                        fprintf('  idx=%d, cellID=%d, animal=%s, min_p=%.4f (12-%.0fms)\n', ...
                            idx, cellID, animal, min_p, window_end_ms);
                    elseif strcmp(g.testing_method, 'single_window')
                        fprintf('  idx=%d, cellID=%d, animal=%s, p=%.4f (%.0f-%.0fms)\n', ...
                            idx, cellID, animal, p_vals(1), g.test_window_start*1000, g.test_window_end*1000);
                    end
                end
            end
        end

        % Print decreased
        fprintf('\nDECREASED (p<0.05, n=%d):\n', comparison_results.(result_field).n_decreased);
        if ~isempty(comparison_results.(result_field).decreased_idx)
            for i = 1:length(comparison_results.(result_field).decreased_idx)
                idx = comparison_results.(result_field).decreased_idx(i);
                cellID = cell_metrics.cellID(idx);
                animal = cell_metrics.animal{idx};

                neuron_pos = find(all_neuron_indices == idx);
                if ~isempty(neuron_pos)
                    p_vals = neuron_pvalues{neuron_pos};
                    if strcmp(g.testing_method, 'multi_window')
                        [min_p, min_idx] = min(p_vals);
                        window_end_ms = (0.013 + (min_idx-1)*0.001) * 1000;
                        fprintf('  idx=%d, cellID=%d, animal=%s, min_p=%.4f (12-%.0fms)\n', ...
                            idx, cellID, animal, min_p, window_end_ms);
                    elseif strcmp(g.testing_method, 'single_window')
                        fprintf('  idx=%d, cellID=%d, animal=%s, p=%.4f (%.0f-%.0fms)\n', ...
                            idx, cellID, animal, p_vals(1), g.test_window_start*1000, g.test_window_end*1000);
                    end
                end
            end
        end

        % Print unchanged
        fprintf('\nUNCHANGED (p>=0.05, n=%d):\n', comparison_results.(result_field).n_unchanged);
        if ~isempty(comparison_results.(result_field).unchanged_idx)
            for i = 1:min(5, length(comparison_results.(result_field).unchanged_idx))  % Only first 5
                idx = comparison_results.(result_field).unchanged_idx(i);
                cellID = cell_metrics.cellID(idx);
                animal = cell_metrics.animal{idx};

                neuron_pos = find(all_neuron_indices == idx);
                if ~isempty(neuron_pos)
                    p_vals = neuron_pvalues{neuron_pos};
                    if strcmp(g.testing_method, 'multi_window')
                        [min_p, min_idx] = min(p_vals);
                        window_end_ms = (0.013 + (min_idx-1)*0.001) * 1000;
                        fprintf('  idx=%d, cellID=%d, animal=%s, min_p=%.4f (12-%.0fms)\n', ...
                            idx, cellID, animal, min_p, window_end_ms);
                    elseif strcmp(g.testing_method, 'single_window')
                        fprintf('  idx=%d, cellID=%d, animal=%s, p=%.4f (%.0f-%.0fms)\n', ...
                            idx, cellID, animal, p_vals(1), g.test_window_start*1000, g.test_window_end*1000);
                    end
                end
            end
            if length(comparison_results.(result_field).unchanged_idx) > 5
                fprintf('  ... and %d more\n', length(comparison_results.(result_field).unchanged_idx) - 5);
            end
        end
    end
end

fprintf('\n========================================\n');
fprintf('To check neurons in UI:\n');
fprintf('  BAfc_monosyn_raster_ui(g_ui, monosyn_results, ttl)\n');
fprintf('========================================\n');
