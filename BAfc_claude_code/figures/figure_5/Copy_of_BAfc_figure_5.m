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
g.artifact_end = 0.010;  % Start of test window (s) - after artifact
g.monosyn_window = 0.05;    % End of test window (s) - monosynaptic window
% Common options:
%   Full monosynaptic: artifact_end=0.012, monosyn_window=0.050
%   Early responses:   artifact_end=0.012, monosyn_window=0.020
%   Sustained:         artifact_end=0.020, monosyn_window=0.050

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
g.zscore_threshold_rule1 = 5;
g.prob_threshold_rule1 = 0.1;
g.zscore_threshold_rule2 = 5;
g.prob_threshold_rule2 = 0.1;

% One-rule responsiveness parameter
g.zscore_threshold_one_rule = 5;

%% Calculate PSTHs for all 4 conditions
psthZ_full = cell(1,4);
psthHz_full = cell(1,4);
baseline_bins = round(g.pre_time / g.bin_time);

% Calculate Savitzky-Golay filter delay correction
filter_delay = floor(g.smoothvalue / 2);

for hmp = 1:4
    [psth_spx_og, preAP_norm_all{hmp}, postAP_norm_all{hmp}] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{hmp}, ...
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

        % TTL indices: CS=1,2; US=3,4
        ttl_nolight_idx = stim*2 - 1;  % 1 for CS, 3 for US
        ttl_light_idx = stim*2;        % 2 for CS, 4 for US

        % Get number of trials
        
        

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

        % Get trial-by-trial spike data (from stored data)
        postAP_norm_nolight = postAP_norm_all{ttl_nolight_idx};
        postAP_norm_light = postAP_norm_all{ttl_light_idx};

        % Process PNs
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
                num_trials_nolight = size(cell_metrics.general.(ttl{ttl_nolight_idx}){global_idx}, 1);
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
                num_trials_light = size(cell_metrics.general.(ttl{ttl_light_idx}){global_idx}, 1);
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

        % Process INs
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
                num_trials_nolight = size(cell_metrics.general.(ttl{ttl_nolight_idx}){global_idx}, 1);
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
                num_trials_light = size(cell_metrics.general.(ttl{ttl_light_idx}){global_idx}, 1);
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


%% Prepare data for BAfc_monosyn_raster_ui

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

end

%% Statistical comparison for CS and US separately

% Setup testing windows based on method
if strcmp(g.testing_method, 'multi_window')
    % Fine-grained multi-window testing
    window_ends = (g.artifact_end + 0.001):0.001:g.monosyn_window;  % artifact+1ms to monosyn_window in 1ms steps
    n_windows = length(window_ends);
elseif strcmp(g.testing_method, 'single_window')
    % Single pre-defined window testing
    window_ends = g.monosyn_window;  % Single window endpoint
    n_windows = 1;
else
    error('Invalid testing_method. Must be ''multi_window'' or ''single_window''');
end

comparison_results = struct();

for br = 1:2
    region = brain_regions{br};

    for stim = 1:2
        stim_name = {'CS', 'US'};
        result_field = sprintf('%s_%s', region, stim_name{stim});


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
                test_window = [g.artifact_end, window_ends(w)];

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

    end
end

%% Visualization: Main figure with nested pie charts and slope graphs
fig_comparison = figure('Units', 'pixels', 'Position', [100, 100, 1000, 600], 'Visible', 'on');
t_main = tiledlayout(fig_comparison, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Create nested tiledlayout for raster plots (rows 1-2, columns 1-2)
t_nested_left = tiledlayout(t_main, 3, 4, 'TileSpacing', 'compact', 'Padding', 'tight');
t_nested_left.Layout.Tile = 1;  % Start at tile 1 (row 1, column 1)
t_nested_left.Layout.TileSpan = [2 2];  % Span 2 rows, 2 columns
title(t_nested_left, 'Example neurons', 'FontSize', 11, 'FontWeight', 'bold');

% Create nested tiledlayout for pie charts and slope graphs (rows 1-2, columns 3-4)
t_nested_right = tiledlayout(t_main, 3, 4, 'TileSpacing', 'tight', 'Padding', 'compact');
t_nested_right.Layout.Tile = 3;  % Start at tile 3 (row 1, column 3)
t_nested_right.Layout.TileSpan = [2 2];  % Span 2 rows, 2 columns

region_names = {'LA', 'Astria'};
stim_names = {'CS', 'US'};
spaghetti_titles = {'LA CS', 'LA US', 'AStria CS', 'AStria US'};

for r = 1:2
    for s = 1:2
        result_field = sprintf('%s_%s', region_names{r}, stim_names{s});

        % Get data
        n_increased = comparison_results.(result_field).n_increased;
        n_decreased = comparison_results.(result_field).n_decreased;
        n_unchanged = comparison_results.(result_field).n_unchanged;
        n_total = comparison_results.(result_field).n_total;

        % Pie chart (row 1 of nested layout) - show all three categories
        col_idx = (r-1)*2 + s;
        ax_pie = nexttile(t_nested_right, col_idx);
        pie_data = [n_increased, n_decreased, n_unchanged];
        colors = [0.8 0.2 0.2; 0.2 0.2 0.8; 0.7 0.7 0.7];

        if sum(pie_data) > 0
            % Create custom labels with percentages only (rounded to whole numbers)
            percentages = 100 * pie_data / sum(pie_data);
            labels = arrayfun(@(p) sprintf('%.0f%%', p), percentages, 'UniformOutput', false);
            h = pie(ax_pie, pie_data, labels);
            colormap(ax_pie, colors);

            % Extract text handles and adjust positions to avoid overlap
            text_handles = h(2:2:end);
            % Adjust positions based on which pie chart (different overlaps for each)
            plot_idx = (r-1)*2 + s;
            if plot_idx == 2  % LA US - adjust overlapping labels (11%, 34%, 54%)
                % Move 34% slightly lower (same distance as AStria CS 16%)
                if length(text_handles) >= 2 && ~isempty(text_handles(2).Position)
                    text_handles(2).Position = text_handles(2).Position + [0 -0.1 0];  % Move 34% slightly lower
                end
                % Move 54% much lower to avoid overlap with AStria CS 16%
                if length(text_handles) >= 3 && ~isempty(text_handles(3).Position)
                    text_handles(3).Position = text_handles(3).Position + [0 -0.3 0];  % Move 54% lower
                end
            elseif plot_idx == 3  % AStria CS - adjust overlapping labels
                % Move 16% slightly to avoid overlap
                if length(text_handles) >= 2 && ~isempty(text_handles(2).Position)
                    text_handles(2).Position = text_handles(2).Position + [0.05 0.1 0];  % Move 16% slightly right-up
                end
                if length(text_handles) >= 3 && ~isempty(text_handles(3).Position)
                    text_handles(3).Position = text_handles(3).Position + [0 -0.15 0];  % Move 66% down
                end
            end
        end

        % Add title using spaghetti_titles with manual positioning
        plot_idx = (r-1)*2 + s;
        title_handle = title(ax_pie, spaghetti_titles{plot_idx}, ...
            'FontSize', g.fontSize1, 'FontWeight', 'bold');
        % Move title up to avoid overlap with percentage labels
        title_handle.Position(2) = title_handle.Position(2) + 0.15;

        % Add n= text below the pie chart
        text(ax_pie, 0, -1.4, sprintf('n=%d', n_total), ...
            'HorizontalAlignment', 'center', 'FontSize', g.fontSize2, ...
            'FontWeight', 'bold');

        % Slope graph - FR for enhanced neurons (row 2 of nested layout)
        ax_slope_enhanced = nexttile(t_nested_right, col_idx + 4);
        n_enhanced = n_increased;

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

            % Calculate FR for each enhanced neuron in their most significant window
            fr_nolight = zeros(n_enhanced, 1);
            fr_light = zeros(n_enhanced, 1);

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
                    test_window = [g.artifact_end, window_end];
                elseif strcmp(g.testing_method, 'single_window')
                    % Use the pre-defined window for all neurons
                    test_window = [g.artifact_end, g.monosyn_window];
                end

                window_duration = test_window(2) - test_window(1);

                % Calculate FR for no-light trials
                spike_count_nolight = 0;
                trial_count_nolight = 0;
                if ~isempty(postAP_norm_nolight{idx})
                    for trial = 1:length(postAP_norm_nolight{idx})
                        trial_spikes = postAP_norm_nolight{idx}{trial};
                        spike_count_nolight = spike_count_nolight + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                    end
                    trial_count_nolight = length(postAP_norm_nolight{idx});
                end
                if trial_count_nolight > 0
                    fr_nolight(i) = (spike_count_nolight / trial_count_nolight) / window_duration;
                end

                % Calculate FR for light trials
                spike_count_light = 0;
                trial_count_light = 0;
                if ~isempty(postAP_norm_light{idx})
                    for trial = 1:length(postAP_norm_light{idx})
                        trial_spikes = postAP_norm_light{idx}{trial};
                        spike_count_light = spike_count_light + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                    end
                    trial_count_light = length(postAP_norm_light{idx});
                end
                if trial_count_light > 0
                    fr_light(i) = (spike_count_light / trial_count_light) / window_duration;
                end
            end

            % Create slope graph for FR (enhanced neurons)
            hold(ax_slope_enhanced, 'on');

            % Plot mean trajectory with thicker line (underneath)
            mean_fr_nolight = mean(fr_nolight);
            mean_fr_light = mean(fr_light);
            plot(ax_slope_enhanced, [1 2], [mean_fr_nolight mean_fr_light], ...
                'k-', 'LineWidth', 3);

            % Plot individual neuron trajectories
            for i = 1:n_enhanced
                plot(ax_slope_enhanced, [1 2], [fr_nolight(i) fr_light(i)], ...
                    '-', 'Color', [0.8 0.2 0.2 0.3], 'LineWidth', 1);
            end

            % Add scatter points on top
            scatter(ax_slope_enhanced, ones(n_enhanced, 1), fr_nolight, 30, ...
                [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
            scatter(ax_slope_enhanced, 2*ones(n_enhanced, 1), fr_light, 30, ...
                [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.6);

            % Add mean markers
            scatter(ax_slope_enhanced, 1, mean_fr_nolight, 80, 'k', 'filled');
            scatter(ax_slope_enhanced, 2, mean_fr_light, 80, 'k', 'filled');

            % Population-level paired t-test on enhanced neurons
            [~, p_pop, ~, ~] = ttest(fr_nolight, fr_light);
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
            y_max = max([fr_nolight; fr_light]);
            y_pos = y_max + y_max * 0.15;
            plot(ax_slope_enhanced, [1 2], [y_pos y_pos], 'k-', 'LineWidth', 1.5);
            text(ax_slope_enhanced, 1.5, y_pos, sig_str, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);

            hold(ax_slope_enhanced, 'off');

            % Set common y-limits for all spaghetti plots
            ylim(ax_slope_enhanced, [0 60]);
            yticks(ax_slope_enhanced, [0 30 60]);

            xlim(ax_slope_enhanced, [0.5 2.5]);
            xticks(ax_slope_enhanced, [1 2]);
            xticklabels(ax_slope_enhanced, {''});
            if r == 1 && s == 1
                ylabel(ax_slope_enhanced, 'FR (Hz)', 'FontSize', g.fontSize2);
            end

            % Add title using spaghetti_titles
            plot_idx = (r-1)*2 + s;
            title(ax_slope_enhanced, spaghetti_titles{plot_idx}, 'FontSize', g.fontSize2);
            set(ax_slope_enhanced, 'FontSize', g.fontSize2);
            box(ax_slope_enhanced, 'off');
        else
            % No enhanced neurons
            text(ax_slope_enhanced, 0.5, 0.5, 'No enhanced neurons', 'HorizontalAlignment', 'center', ...
                'FontSize', g.fontSize2);
            axis(ax_slope_enhanced, 'off');
        end

        % Slope graph - FR for decreased neurons (row 3 of nested layout)
        ax_slope_decreased = nexttile(t_nested_right, col_idx + 8);
        n_decreased = comparison_results.(result_field).n_decreased;

        if n_decreased > 0
            % Get decreased neuron indices
            decreased_idx = comparison_results.(result_field).decreased_idx;
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

            % Calculate FR for each decreased neuron in their most significant window
            fr_nolight_dec = zeros(n_decreased, 1);
            fr_light_dec = zeros(n_decreased, 1);

            for i = 1:n_decreased
                idx = decreased_idx(i);

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
                    test_window = [g.artifact_end, window_end];
                elseif strcmp(g.testing_method, 'single_window')
                    % Use the pre-defined window for all neurons
                    test_window = [g.artifact_end, g.monosyn_window];
                end

                window_duration = test_window(2) - test_window(1);

                % Calculate FR for no-light trials
                spike_count_nolight = 0;
                trial_count_nolight = 0;
                if ~isempty(postAP_norm_nolight{idx})
                    for trial = 1:length(postAP_norm_nolight{idx})
                        trial_spikes = postAP_norm_nolight{idx}{trial};
                        spike_count_nolight = spike_count_nolight + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                    end
                    trial_count_nolight = length(postAP_norm_nolight{idx});
                end
                if trial_count_nolight > 0
                    fr_nolight_dec(i) = (spike_count_nolight / trial_count_nolight) / window_duration;
                end

                % Calculate FR for light trials
                spike_count_light = 0;
                trial_count_light = 0;
                if ~isempty(postAP_norm_light{idx})
                    for trial = 1:length(postAP_norm_light{idx})
                        trial_spikes = postAP_norm_light{idx}{trial};
                        spike_count_light = spike_count_light + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                    end
                    trial_count_light = length(postAP_norm_light{idx});
                end
                if trial_count_light > 0
                    fr_light_dec(i) = (spike_count_light / trial_count_light) / window_duration;
                end
            end

            % Create slope graph for FR (decreased neurons)
            hold(ax_slope_decreased, 'on');

            % Plot mean trajectory with thicker line (underneath)
            mean_fr_nolight_dec = mean(fr_nolight_dec);
            mean_fr_light_dec = mean(fr_light_dec);
            plot(ax_slope_decreased, [1 2], [mean_fr_nolight_dec mean_fr_light_dec], ...
                'k-', 'LineWidth', 3);

            % Plot individual neuron trajectories
            for i = 1:n_decreased
                plot(ax_slope_decreased, [1 2], [fr_nolight_dec(i) fr_light_dec(i)], ...
                    '-', 'Color', [0.2 0.2 0.8 0.3], 'LineWidth', 1);
            end

            % Add scatter points on top
            scatter(ax_slope_decreased, ones(n_decreased, 1), fr_nolight_dec, 30, ...
                [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
            scatter(ax_slope_decreased, 2*ones(n_decreased, 1), fr_light_dec, 30, ...
                [0.2 0.2 0.8], 'filled', 'MarkerFaceAlpha', 0.6);

            % Add mean markers
            scatter(ax_slope_decreased, 1, mean_fr_nolight_dec, 80, 'k', 'filled');
            scatter(ax_slope_decreased, 2, mean_fr_light_dec, 80, 'k', 'filled');

            % Population-level paired t-test on decreased neurons
            [~, p_pop_dec, ~, ~] = ttest(fr_nolight_dec, fr_light_dec);
            if p_pop_dec < 0.001
                sig_str_dec = '***';
            elseif p_pop_dec < 0.01
                sig_str_dec = '**';
            elseif p_pop_dec < 0.05
                sig_str_dec = '*';
            else
                sig_str_dec = 'n.s.';
            end

            % Add significance line and text above the plot
            y_max_dec = max([fr_nolight_dec; fr_light_dec]);
            y_pos_dec = y_max_dec + y_max_dec * 0.15;
            plot(ax_slope_decreased, [1 2], [y_pos_dec y_pos_dec], 'k-', 'LineWidth', 1.5);
            text(ax_slope_decreased, 1.5, y_pos_dec, sig_str_dec, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);

            hold(ax_slope_decreased, 'off');

            % Set y-limits for decreased spaghetti plots
            ylim(ax_slope_decreased, [0 80]);
            yticks(ax_slope_decreased, [0 40 80]);

            xlim(ax_slope_decreased, [0.5 2.5]);
            xticks(ax_slope_decreased, [1 2]);
            xticklabels(ax_slope_decreased, {'No light', 'Light'});
            if r == 1 && s == 1
                ylabel(ax_slope_decreased, 'FR (Hz)', 'FontSize', g.fontSize2);
            end

            set(ax_slope_decreased, 'FontSize', g.fontSize2);
            box(ax_slope_decreased, 'off');
        else
            % No decreased neurons
            text(ax_slope_decreased, 0.5, 0.5, 'No decreased neurons', 'HorizontalAlignment', 'center', ...
                'FontSize', g.fontSize2);
            axis(ax_slope_decreased, 'off');
        end

    end
end

%% Add example raster plots (rows 1-2 of left nested layout)
% Define example neurons: [animal, cellID, stimulus_type, stimulus_type_light]
% examples = {
%     'MD309_001', 20, 'triptest_sound_only', 'triptest_sound_only_light';   % CS-only
%     'MD309_001', 20, 'triptest_shocks_only', 'triptest_shocks_only_light';  % US-only
%     'MD318_001', 46, 'triptest_sound_only', 'triptest_sound_only_light';   % CS-only
%     'MD317_001', 43, 'triptest_shocks_only', 'triptest_shocks_only_light'   % US-only
% };
examples = {
    'MD309_001', 20, 'triptest_shocks_only', 'triptest_shocks_only_light';  % LA US UP
    'MD317_001', 43, 'triptest_shocks_only', 'triptest_shocks_only_light';   % AStria US UP
    'MD307_001', 16, 'triptest_shocks_only', 'triptest_shocks_only_light'   % LA US DOWN
    'MD315_001', 15, 'triptest_shocks_only', 'triptest_shocks_only_light'   % AStria US DOWN
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
        % Get spike times from stored data
        ttl_idx = find(strcmp(ttl, ttl_type));
        preAP_norm = preAP_norm_all{ttl_idx};
        postAP_norm = postAP_norm_all{ttl_idx};

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
            raster_titles = {'LA US ↑', 'AStria US ↑', 'LA US ↓', 'AStria US ↓'};
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
        % Get spike times from stored data for light condition
        ttl_idx_light = find(strcmp(ttl, ttl_type_light));
        preAP_norm_light = preAP_norm_all{ttl_idx_light};
        postAP_norm_light = postAP_norm_all{ttl_idx_light};

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

            % Apply 11ms smoothing for visualization only
            smoothvalue_vis = 11;
            filter_delay_vis = floor(smoothvalue_vis / 2);

            % Smooth no-light data
            neuron_psth_nolight_smooth = smoothdata(neuron_psth_nolight, 2, 'sgolay', smoothvalue_vis);
            neuron_psth_nolight_corrected = zeros(size(neuron_psth_nolight_smooth));
            neuron_psth_nolight_corrected(filter_delay_vis+1:end) = neuron_psth_nolight_smooth(1:end-filter_delay_vis);
            neuron_psth_nolight_corrected(1:filter_delay_vis) = repmat(neuron_psth_nolight_smooth(1), 1, filter_delay_vis);

            % Smooth light data
            neuron_psth_light_smooth = smoothdata(neuron_psth_light, 2, 'sgolay', smoothvalue_vis);
            neuron_psth_light_corrected = zeros(size(neuron_psth_light_smooth));
            neuron_psth_light_corrected(filter_delay_vis+1:end) = neuron_psth_light_smooth(1:end-filter_delay_vis);
            neuron_psth_light_corrected(1:filter_delay_vis) = repmat(neuron_psth_light_smooth(1), 1, filter_delay_vis);

            % Time axis in ms, centered at stimulus
            time_axis = ((-g.pre_time:g.bin_time:(g.post_time-g.bin_time)) * 1000);

            % Plot both conditions
            hold(ax_lineplot, 'on');
            % Use different colors for decreased examples (ex 3,4) vs increased (ex 1,2)
            if ex >= 3
                % Decreased: black and blue
                plot(ax_lineplot, time_axis, neuron_psth_nolight_corrected, 'k-', 'LineWidth', 1.5);
                plot(ax_lineplot, time_axis, neuron_psth_light_corrected, 'Color', [0.2 0.2 0.8], 'LineWidth', 1.5);
            else
                % Increased: black and red
                plot(ax_lineplot, time_axis, neuron_psth_nolight_corrected, 'k-', 'LineWidth', 1.5);
                plot(ax_lineplot, time_axis, neuron_psth_light_corrected, 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
            end

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

            % Add legend for first and third plots (increased and decreased)
            if ex == 1
                leg = legend(ax_lineplot, {'No light', 'Laser on'}, 'Location', 'northeast', ...
                    'FontSize', g.fontSize2, 'Box', 'off');
                leg.ItemTokenSize = [10, 12];  % Shorter lines
                leg.Position(1) = leg.Position(1) + 0.06;  % Shift right
                leg.Position(2) = leg.Position(2) + 0.03;  % Shift up
            elseif ex == 3
                leg = legend(ax_lineplot, {'No light', 'Laser on'}, 'Location', 'northeast', ...
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
annotation(fig_comparison, 'textbox', [0.55, 0.93, 0.4, 0.04], ...
    'String', 'Light-modulated neurons', ...
    'FontSize', 11, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none');

% Add legend below the title using square markers
ax_legend = axes(fig_comparison, 'Position', [0.55, 0.90, 0.4, 0.04], 'Visible', 'off');
hold(ax_legend, 'on');
% Create invisible plot objects with square markers for legend
h1 = plot(ax_legend, NaN, NaN, 's', 'MarkerSize', 12, 'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerEdgeColor', 'k', 'LineWidth', 1);
h2 = plot(ax_legend, NaN, NaN, 's', 'MarkerSize', 12, 'MarkerFaceColor', [0.2 0.2 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 1);
h3 = plot(ax_legend, NaN, NaN, 's', 'MarkerSize', 12, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'k', 'LineWidth', 1);
leg = legend(ax_legend, [h1 h2 h3], {'Increased', 'Decreased', 'Unchanged'}, ...
    'Orientation', 'horizontal', 'FontSize', g.fontSize2, ...
    'Location', 'north', 'Box', 'off');
leg.ItemTokenSize = [30, 30];  % Increase marker size in legend
hold(ax_legend, 'off');

%% Add panel labels
% C: Rasterplots (top-left of left panel)
annotation(fig_comparison, 'textbox', [0.01, 0.97, 0.03, 0.04], ...
    'String', 'C', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% D: FR comparison lineplots (bottom-left of left panel, row 3)
annotation(fig_comparison, 'textbox', [0.01, 0.35, 0.03, 0.04], ...
    'String', 'D', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% E: Pie charts (top-right of right panel)
annotation(fig_comparison, 'textbox', [0.51, 0.97, 0.03, 0.04], ...
    'String', 'E', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% F: Enhanced spaghetti plots (middle-right of right panel)
annotation(fig_comparison, 'textbox', [0.51, 0.63, 0.03, 0.04], ...
    'String', 'F', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% G: Decreased spaghetti plots (bottom-right of right panel)
annotation(fig_comparison, 'textbox', [0.51, 0.29, 0.03, 0.04], ...
    'String', 'G', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');


%% Export data
export_figure5_to_excel(comparison_results, cell_metrics, results_all, region_names, stim_names, g);
exportgraphics(gcf, 'figure_5.png', 'Resolution', 300);
%% Excel Export Function
function export_figure5_to_excel(comparison_results, cell_metrics, results_all, region_names, stim_names, g)
    output_filename = 'figure_5_data.xlsx';

    if exist(output_filename, 'file')
        delete(output_filename);
    end

    % Helper function for significance markers
    function sig_str = format_significance(p_val)
        if p_val < 0.001
            sig_str = '***';
        elseif p_val < 0.01
            sig_str = '**';
        elseif p_val < 0.05
            sig_str = '*';
        else
            sig_str = 'n.s.';
        end
    end

    % Window ends for multi-window testing
    if strcmp(g.testing_method, 'multi_window')
        window_ends = (g.artifact_end + 0.001):0.001:g.monosyn_window;
    end

    %% Sheet 1: Panel E - Pie Charts (Proportions)
    sheet_data = {};
    sheet_data{1,1} = 'PANEL E: PIE CHARTS - PROPORTIONS OF ENHANCED NEURONS';
    sheet_data{2,1} = sprintf('Testing method: %s', g.testing_method);
    if strcmp(g.testing_method, 'single_window')
        sheet_data{3,1} = sprintf('Test window: %.0f-%.0f ms', g.artifact_end*1000, g.monosyn_window*1000);
    end
    sheet_data{5,1} = 'Region';
    sheet_data{5,2} = 'Stimulus';
    sheet_data{5,3} = 'Total Monosynaptic';
    sheet_data{5,4} = 'Enhanced';
    sheet_data{5,5} = 'Enhanced (%)';
    sheet_data{5,6} = 'Non-Enhanced';
    sheet_data{5,7} = 'Non-Enhanced (%)';

    row_idx = 6;
    for r = 1:2
        for s = 1:2
            result_field = sprintf('%s_%s', region_names{r}, stim_names{s});
            result = comparison_results.(result_field);

            n_enhanced = result.n_increased;
            n_non_enhanced = result.n_decreased + result.n_unchanged;

            sheet_data{row_idx, 1} = region_names{r};
            sheet_data{row_idx, 2} = stim_names{s};
            sheet_data{row_idx, 3} = result.n_total;
            sheet_data{row_idx, 4} = n_enhanced;
            if result.n_total > 0
                sheet_data{row_idx, 5} = 100 * n_enhanced / result.n_total;
            else
                sheet_data{row_idx, 5} = 0;
            end
            sheet_data{row_idx, 6} = n_non_enhanced;
            if result.n_total > 0
                sheet_data{row_idx, 7} = 100 * n_non_enhanced / result.n_total;
            else
                sheet_data{row_idx, 7} = 0;
            end
            row_idx = row_idx + 1;
        end
    end

    writecell(sheet_data, output_filename, 'Sheet', 'PanelE_PieCharts');

    %% Sheet 2: Panel F - Spaghetti Plots (FR)
    sheet_data = {};
    sheet_data{1,1} = 'PANEL F: SPAGHETTI PLOTS - FR FOR ENHANCED NEURONS';
    sheet_data{2,1} = 'FR calculated as (nSpikes/trial) / window_duration (no baseline subtraction)';
    sheet_data{4,1} = 'Region';
    sheet_data{4,2} = 'Stimulus';
    sheet_data{4,3} = 'N Enhanced';
    sheet_data{4,4} = 'No Light Mean (Hz)';
    sheet_data{4,5} = 'No Light SEM (Hz)';
    sheet_data{4,6} = 'No Light Median (Hz)';
    sheet_data{4,7} = 'No Light SD (Hz)';
    sheet_data{4,8} = '';
    sheet_data{4,9} = 'With Light Mean (Hz)';
    sheet_data{4,10} = 'With Light SEM (Hz)';
    sheet_data{4,11} = 'With Light Median (Hz)';
    sheet_data{4,12} = 'With Light SD (Hz)';
    sheet_data{4,13} = '';
    sheet_data{4,14} = 'Paired t-test p';
    sheet_data{4,15} = 'Significance';

    row_idx = 5;
    for r = 1:2
        for s = 1:2
            result_field = sprintf('%s_%s', region_names{r}, stim_names{s});
            result = comparison_results.(result_field);

            sheet_data{row_idx, 1} = region_names{r};
            sheet_data{row_idx, 2} = stim_names{s};
            sheet_data{row_idx, 3} = result.n_increased;

            if result.n_increased > 0
                res = results_all{r, s};
                postAP_norm_nolight = res.postAP_norm_nolight;
                postAP_norm_light = res.postAP_norm_light;

                increased_idx = result.increased_idx;
                all_neuron_indices = result.all_neuron_indices;
                neuron_pvalues = result.neuron_pvalues;
                n_enhanced = result.n_increased;

                fr_nolight = zeros(n_enhanced, 1);
                fr_light = zeros(n_enhanced, 1);

                for i = 1:n_enhanced
                    idx = increased_idx(i);

                    if strcmp(g.testing_method, 'multi_window')
                        neuron_pos = find(all_neuron_indices == idx);
                        if ~isempty(neuron_pos)
                            p_vals = neuron_pvalues{neuron_pos};
                            [~, min_window_idx] = min(p_vals);
                            window_end = window_ends(min_window_idx);
                        else
                            window_end = window_ends(1);
                        end
                        test_window = [g.artifact_end, window_end];
                    elseif strcmp(g.testing_method, 'single_window')
                        test_window = [g.artifact_end, g.monosyn_window];
                    end

                    window_duration = test_window(2) - test_window(1);

                    % Calculate FR for no-light
                    spike_count_nolight = 0;
                    trial_count_nolight = 0;
                    if ~isempty(postAP_norm_nolight{idx})
                        for trial = 1:length(postAP_norm_nolight{idx})
                            trial_spikes = postAP_norm_nolight{idx}{trial};
                            spike_count_nolight = spike_count_nolight + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                        end
                        trial_count_nolight = length(postAP_norm_nolight{idx});
                    end
                    if trial_count_nolight > 0
                        fr_nolight(i) = (spike_count_nolight / trial_count_nolight) / window_duration;
                    end

                    % Calculate FR for light
                    spike_count_light = 0;
                    trial_count_light = 0;
                    if ~isempty(postAP_norm_light{idx})
                        for trial = 1:length(postAP_norm_light{idx})
                            trial_spikes = postAP_norm_light{idx}{trial};
                            spike_count_light = spike_count_light + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                        end
                        trial_count_light = length(postAP_norm_light{idx});
                    end
                    if trial_count_light > 0
                        fr_light(i) = (spike_count_light / trial_count_light) / window_duration;
                    end
                end

                [~, p_pop, ~, ~] = ttest(fr_nolight, fr_light);

                sheet_data{row_idx, 4} = mean(fr_nolight);
                sheet_data{row_idx, 5} = std(fr_nolight) / sqrt(length(fr_nolight));
                sheet_data{row_idx, 6} = median(fr_nolight);
                sheet_data{row_idx, 7} = std(fr_nolight);
                sheet_data{row_idx, 9} = mean(fr_light);
                sheet_data{row_idx, 10} = std(fr_light) / sqrt(length(fr_light));
                sheet_data{row_idx, 11} = median(fr_light);
                sheet_data{row_idx, 12} = std(fr_light);
                sheet_data{row_idx, 14} = p_pop;
                sheet_data{row_idx, 15} = format_significance(p_pop);
            end

            row_idx = row_idx + 1;
        end
    end

    writecell(sheet_data, output_filename, 'Sheet', 'PanelF_SpaghetttiPlots');

    %% Sheet 3: Raw Data - Individual Neurons
    sheet_data = {};
    sheet_data{1,1} = 'RAW DATA: INDIVIDUAL NEURON CLASSIFICATIONS';
    sheet_data{3,1} = 'Global Index';
    sheet_data{3,2} = 'Cell ID';
    sheet_data{3,3} = 'Animal ID';
    sheet_data{3,4} = 'Region';
    sheet_data{3,5} = 'Cell Type';
    sheet_data{3,6} = 'Stimulus';
    sheet_data{3,7} = 'Classification';
    sheet_data{3,8} = 'p-value';
    sheet_data{3,9} = 'Test Window (ms)';
    sheet_data{3,10} = '';
    sheet_data{3,11} = 'FR No Light (Hz)';
    sheet_data{3,12} = 'FR Light (Hz)';
    sheet_data{3,13} = 'FR Diff (Hz)';

    % Window ends for multi-window testing
    if strcmp(g.testing_method, 'multi_window')
        window_ends = (g.artifact_end + 0.001):0.001:g.monosyn_window;
    end

    row_idx = 4;
    for r = 1:2
        for s = 1:2
            result_field = sprintf('%s_%s', region_names{r}, stim_names{s});
            result = comparison_results.(result_field);

            if result.n_total == 0
                continue;
            end

            % Get data from results_all
            res = results_all{r, s};
            postAP_norm_nolight = res.postAP_norm_nolight;
            postAP_norm_light = res.postAP_norm_light;

            % Process increased
            for i = 1:length(result.increased_idx)
                idx = result.increased_idx(i);
                neuron_pos = find(result.all_neuron_indices == idx);

                sheet_data{row_idx, 1} = idx;
                sheet_data{row_idx, 2} = cell_metrics.cellID(idx);
                sheet_data{row_idx, 3} = cell_metrics.animal{idx};
                sheet_data{row_idx, 4} = region_names{r};
                sheet_data{row_idx, 5} = cell_metrics.putativeCellType{idx};
                sheet_data{row_idx, 6} = stim_names{s};
                sheet_data{row_idx, 7} = 'Increased';

                if ~isempty(neuron_pos)
                    p_vals = result.neuron_pvalues{neuron_pos};
                    if strcmp(g.testing_method, 'multi_window')
                        [min_p, min_idx_val] = min(p_vals);
                        window_end = window_ends(min_idx_val);
                        window_end_ms = window_end * 1000;
                        test_window = [g.artifact_end, window_end];
                        sheet_data{row_idx, 8} = min_p;
                        sheet_data{row_idx, 9} = sprintf('%.0f-%.0f', g.artifact_end*1000, window_end_ms);
                    elseif strcmp(g.testing_method, 'single_window')
                        test_window = [g.artifact_end, g.monosyn_window];
                        sheet_data{row_idx, 8} = p_vals(1);
                        sheet_data{row_idx, 9} = sprintf('%.0f-%.0f', g.artifact_end*1000, g.monosyn_window*1000);
                    end

                    % Calculate FR
                    window_duration = test_window(2) - test_window(1);

                    spike_count_nolight = 0;
                    trial_count_nolight = 0;
                    if ~isempty(postAP_norm_nolight{idx})
                        for trial = 1:length(postAP_norm_nolight{idx})
                            trial_spikes = postAP_norm_nolight{idx}{trial};
                            spike_count_nolight = spike_count_nolight + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                        end
                        trial_count_nolight = length(postAP_norm_nolight{idx});
                    end
                    fr_nolight = 0;
                    if trial_count_nolight > 0
                        fr_nolight = (spike_count_nolight / trial_count_nolight) / window_duration;
                    end

                    spike_count_light = 0;
                    trial_count_light = 0;
                    if ~isempty(postAP_norm_light{idx})
                        for trial = 1:length(postAP_norm_light{idx})
                            trial_spikes = postAP_norm_light{idx}{trial};
                            spike_count_light = spike_count_light + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                        end
                        trial_count_light = length(postAP_norm_light{idx});
                    end
                    fr_light = 0;
                    if trial_count_light > 0
                        fr_light = (spike_count_light / trial_count_light) / window_duration;
                    end

                    sheet_data{row_idx, 11} = fr_nolight;
                    sheet_data{row_idx, 12} = fr_light;
                    sheet_data{row_idx, 13} = fr_light - fr_nolight;
                end
                row_idx = row_idx + 1;
            end

            % Process decreased
            for i = 1:length(result.decreased_idx)
                idx = result.decreased_idx(i);
                neuron_pos = find(result.all_neuron_indices == idx);

                sheet_data{row_idx, 1} = idx;
                sheet_data{row_idx, 2} = cell_metrics.cellID(idx);
                sheet_data{row_idx, 3} = cell_metrics.animal{idx};
                sheet_data{row_idx, 4} = region_names{r};
                sheet_data{row_idx, 5} = cell_metrics.putativeCellType{idx};
                sheet_data{row_idx, 6} = stim_names{s};
                sheet_data{row_idx, 7} = 'Decreased';

                if ~isempty(neuron_pos)
                    p_vals = result.neuron_pvalues{neuron_pos};
                    if strcmp(g.testing_method, 'multi_window')
                        [min_p, min_idx_val] = min(p_vals);
                        window_end = window_ends(min_idx_val);
                        window_end_ms = window_end * 1000;
                        test_window = [g.artifact_end, window_end];
                        sheet_data{row_idx, 8} = min_p;
                        sheet_data{row_idx, 9} = sprintf('%.0f-%.0f', g.artifact_end*1000, window_end_ms);
                    elseif strcmp(g.testing_method, 'single_window')
                        test_window = [g.artifact_end, g.monosyn_window];
                        sheet_data{row_idx, 8} = p_vals(1);
                        sheet_data{row_idx, 9} = sprintf('%.0f-%.0f', g.artifact_end*1000, g.monosyn_window*1000);
                    end

                    % Calculate FR
                    window_duration = test_window(2) - test_window(1);

                    spike_count_nolight = 0;
                    trial_count_nolight = 0;
                    if ~isempty(postAP_norm_nolight{idx})
                        for trial = 1:length(postAP_norm_nolight{idx})
                            trial_spikes = postAP_norm_nolight{idx}{trial};
                            spike_count_nolight = spike_count_nolight + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                        end
                        trial_count_nolight = length(postAP_norm_nolight{idx});
                    end
                    fr_nolight = 0;
                    if trial_count_nolight > 0
                        fr_nolight = (spike_count_nolight / trial_count_nolight) / window_duration;
                    end

                    spike_count_light = 0;
                    trial_count_light = 0;
                    if ~isempty(postAP_norm_light{idx})
                        for trial = 1:length(postAP_norm_light{idx})
                            trial_spikes = postAP_norm_light{idx}{trial};
                            spike_count_light = spike_count_light + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                        end
                        trial_count_light = length(postAP_norm_light{idx});
                    end
                    fr_light = 0;
                    if trial_count_light > 0
                        fr_light = (spike_count_light / trial_count_light) / window_duration;
                    end

                    sheet_data{row_idx, 11} = fr_nolight;
                    sheet_data{row_idx, 12} = fr_light;
                    sheet_data{row_idx, 13} = fr_light - fr_nolight;
                end
                row_idx = row_idx + 1;
            end

            % Process unchanged
            for i = 1:length(result.unchanged_idx)
                idx = result.unchanged_idx(i);
                neuron_pos = find(result.all_neuron_indices == idx);

                sheet_data{row_idx, 1} = idx;
                sheet_data{row_idx, 2} = cell_metrics.cellID(idx);
                sheet_data{row_idx, 3} = cell_metrics.animal{idx};
                sheet_data{row_idx, 4} = region_names{r};
                sheet_data{row_idx, 5} = cell_metrics.putativeCellType{idx};
                sheet_data{row_idx, 6} = stim_names{s};
                sheet_data{row_idx, 7} = 'Unchanged';

                if ~isempty(neuron_pos)
                    p_vals = result.neuron_pvalues{neuron_pos};
                    if strcmp(g.testing_method, 'multi_window')
                        [min_p, min_idx_val] = min(p_vals);
                        window_end = window_ends(min_idx_val);
                        window_end_ms = window_end * 1000;
                        test_window = [g.artifact_end, window_end];
                        sheet_data{row_idx, 8} = min_p;
                        sheet_data{row_idx, 9} = sprintf('%.0f-%.0f', g.artifact_end*1000, window_end_ms);
                    elseif strcmp(g.testing_method, 'single_window')
                        test_window = [g.artifact_end, g.monosyn_window];
                        sheet_data{row_idx, 8} = p_vals(1);
                        sheet_data{row_idx, 9} = sprintf('%.0f-%.0f', g.artifact_end*1000, g.monosyn_window*1000);
                    end

                    % Calculate FR
                    window_duration = test_window(2) - test_window(1);

                    spike_count_nolight = 0;
                    trial_count_nolight = 0;
                    if ~isempty(postAP_norm_nolight{idx})
                        for trial = 1:length(postAP_norm_nolight{idx})
                            trial_spikes = postAP_norm_nolight{idx}{trial};
                            spike_count_nolight = spike_count_nolight + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                        end
                        trial_count_nolight = length(postAP_norm_nolight{idx});
                    end
                    fr_nolight = 0;
                    if trial_count_nolight > 0
                        fr_nolight = (spike_count_nolight / trial_count_nolight) / window_duration;
                    end

                    spike_count_light = 0;
                    trial_count_light = 0;
                    if ~isempty(postAP_norm_light{idx})
                        for trial = 1:length(postAP_norm_light{idx})
                            trial_spikes = postAP_norm_light{idx}{trial};
                            spike_count_light = spike_count_light + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
                        end
                        trial_count_light = length(postAP_norm_light{idx});
                    end
                    fr_light = 0;
                    if trial_count_light > 0
                        fr_light = (spike_count_light / trial_count_light) / window_duration;
                    end

                    sheet_data{row_idx, 11} = fr_nolight;
                    sheet_data{row_idx, 12} = fr_light;
                    sheet_data{row_idx, 13} = fr_light - fr_nolight;
                end
                row_idx = row_idx + 1;
            end
        end
    end

    writecell(sheet_data, output_filename, 'Sheet', 'RawData_Classifications');

    %% Sheet 4: Chi-square tests - LA vs AStria comparisons
    sheet_data = {};
    sheet_data{1,1} = 'CHI-SQUARE TESTS: LA VS ASTRIA';
    sheet_data{2,1} = 'Comparison of classification distributions between brain regions';

    % Get data for both tests
    LA_CS = comparison_results.LA_CS;
    LA_US = comparison_results.LA_US;
    Astria_CS = comparison_results.Astria_CS;
    Astria_US = comparison_results.Astria_US;

    % Test 1: LA CS vs AStria CS
    sheet_data{4,1} = 'TEST 1: LA CS vs AStria CS';
    sheet_data{5,1} = '';
    sheet_data{6,1} = 'Contingency Table:';
    sheet_data{7,1} = '';
    sheet_data{7,2} = 'Increased';
    sheet_data{7,3} = 'Decreased';
    sheet_data{7,4} = 'Unchanged';
    sheet_data{7,5} = 'Total';

    sheet_data{8,1} = 'LA CS';
    sheet_data{8,2} = LA_CS.n_increased;
    sheet_data{8,3} = LA_CS.n_decreased;
    sheet_data{8,4} = LA_CS.n_unchanged;
    sheet_data{8,5} = LA_CS.n_total;

    sheet_data{9,1} = 'AStria CS';
    sheet_data{9,2} = Astria_CS.n_increased;
    sheet_data{9,3} = Astria_CS.n_decreased;
    sheet_data{9,4} = Astria_CS.n_unchanged;
    sheet_data{9,5} = Astria_CS.n_total;

    % Perform chi-square test for CS
    contingency_table_CS = [LA_CS.n_increased, LA_CS.n_decreased, LA_CS.n_unchanged;
                            Astria_CS.n_increased, Astria_CS.n_decreased, Astria_CS.n_unchanged];
    [chi2_CS, p_CS] = calculate_chi_square(contingency_table_CS);
    df_CS = (size(contingency_table_CS,1)-1) * (size(contingency_table_CS,2)-1);

    sheet_data{11,1} = 'Chi-square statistic:';
    sheet_data{11,2} = chi2_CS;
    sheet_data{12,1} = 'Degrees of freedom:';
    sheet_data{12,2} = df_CS;
    sheet_data{13,1} = 'p-value:';
    sheet_data{13,2} = p_CS;
    sheet_data{13,3} = format_significance(p_CS);

    % Test 2: LA US vs AStria US
    sheet_data{16,1} = 'TEST 2: LA US vs AStria US';
    sheet_data{17,1} = '';
    sheet_data{18,1} = 'Contingency Table:';
    sheet_data{19,1} = '';
    sheet_data{19,2} = 'Increased';
    sheet_data{19,3} = 'Decreased';
    sheet_data{19,4} = 'Unchanged';
    sheet_data{19,5} = 'Total';

    sheet_data{20,1} = 'LA US';
    sheet_data{20,2} = LA_US.n_increased;
    sheet_data{20,3} = LA_US.n_decreased;
    sheet_data{20,4} = LA_US.n_unchanged;
    sheet_data{20,5} = LA_US.n_total;

    sheet_data{21,1} = 'AStria US';
    sheet_data{21,2} = Astria_US.n_increased;
    sheet_data{21,3} = Astria_US.n_decreased;
    sheet_data{21,4} = Astria_US.n_unchanged;
    sheet_data{21,5} = Astria_US.n_total;

    % Perform chi-square test for US only using contingency table
    contingency_table_US = [LA_US.n_increased, LA_US.n_decreased, LA_US.n_unchanged;
                            Astria_US.n_increased, Astria_US.n_decreased, Astria_US.n_unchanged];
    [chi2_US, p_US] = calculate_chi_square(contingency_table_US);
    df_US = (size(contingency_table_US,1)-1) * (size(contingency_table_US,2)-1);

    sheet_data{23,1} = 'Chi-square statistic:';
    sheet_data{23,2} = chi2_US;
    sheet_data{24,1} = 'Degrees of freedom:';
    sheet_data{24,2} = df_US;
    sheet_data{25,1} = 'p-value:';
    sheet_data{25,2} = p_US;
    sheet_data{25,3} = format_significance(p_US);

    writecell(sheet_data, output_filename, 'Sheet', 'ChiSquare_RegionComparisons');
end
