% BAfc_figure_optogenetics.m
% Publication-quality figure showing optogenetic PV silencing effects on monosynaptic responses
%
% This version uses BAfc_identify_responsive_neurons() function for neuron identification
%
% PURPOSE:
% Generate publication-style figure showing how optogenetic silencing of PV interneurons
% affects monosynaptic responses (0-25ms) to shock (US) and tone (CS) stimuli.
%
% FIGURE LAYOUT:
% Figure 1 (Original): Example neurons and light inhibition verification
% Figure 2 (Short Latency): Population-level short latency response analysis
%
% Author: Claude Code
% Date: 2025-10-19

clear all
%% ===== USER CONFIGURATION =====

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

% Stimulus selection
ttl = {'triptest_sound_only', 'triptest_sound_only_light', ...
    'triptest_shocks_only', 'triptest_shocks_only_light'};

%% ===== LOAD DATA =====

fprintf('\nLoading neural data...\n');
g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
g.cell_metrics = BAfc_match_celltypes('cell_metrics', g.cell_metrics);

%% ===== PARAMETERS =====

fprintf('\nLoading default parameters...\n');

% Get default parameters (modify BAfc_optogenetics_params.m to change defaults)
params = BAfc_optogenetics_params();
g.params = params;
% Override specific parameters here if needed for this figure
% params.zscore_threshold = 8;  % Example: use lower threshold
% params.bR = 'LA';             % Example: only analyze LA neurons

fprintf('\nIdentifying responsive neurons...\n');
responsive_results = BAfc_identify_responsive_neurons(g.cell_metrics, ttl, params);

%% ===== EXTRACT RESPONSIVE NEURON INDICES =====

% TTL field names
ttl_CS_nolight = matlab.lang.makeValidName(ttl{1});     % triptest_sound_only
ttl_CS_light = matlab.lang.makeValidName(ttl{2});       % triptest_sound_only_light
ttl_US_nolight = matlab.lang.makeValidName(ttl{3});     % triptest_shocks_only
ttl_US_light = matlab.lang.makeValidName(ttl{4});       % triptest_shocks_only_light

% Get responsive neurons for each condition
responsive_CS_nolight_idx = responsive_results.(ttl_CS_nolight).neuron_idx;
responsive_CS_light_idx = responsive_results.(ttl_CS_light).neuron_idx;
responsive_US_nolight_idx = responsive_results.(ttl_US_nolight).neuron_idx;
responsive_US_light_idx = responsive_results.(ttl_US_light).neuron_idx;

% Determine CS, US, and CSUS responsive neurons
% CS responsive: responsive to CS during no-light OR light
CS_responsive_idx_unique = unique([responsive_CS_nolight_idx; responsive_CS_light_idx]);
idx_CS_responsive = CS_responsive_idx_unique;

% US responsive: responsive to US during no-light OR light
US_responsive_idx_unique = unique([responsive_US_nolight_idx; responsive_US_light_idx]);
idx_US_responsive = US_responsive_idx_unique;

% CSUS responsive: responsive to both CS and US
CSUS_responsive_logical = ismember(1:numel(g.cell_metrics.cellID), idx_CS_responsive) & ...
                          ismember(1:numel(g.cell_metrics.cellID), idx_US_responsive);
idx_CSUS_responsive = find(CSUS_responsive_logical);

% CS-only and US-only
idx_CS_only = setdiff(idx_CS_responsive, idx_US_responsive);
idx_US_only = setdiff(idx_US_responsive, idx_CS_responsive);

fprintf('\n===== RESPONSIVE NEURON SUMMARY =====\n');
fprintf('CS responsive: %d neurons\n', length(idx_CS_responsive));
fprintf('US responsive: %d neurons\n', length(idx_US_responsive));
fprintf('CSUS responsive: %d neurons\n', length(idx_CSUS_responsive));
fprintf('CS-only: %d neurons\n', length(idx_CS_only));
fprintf('US-only: %d neurons\n', length(idx_US_only));

%% ===== SELECT EXAMPLE NEURONS =====

% Use manually selected example neurons
idx_example_IN = 97;
idx_example_PN = 423;

fprintf('\nExample neurons (manually selected):\n');
fprintf('  PN: neuron %d (%s, %s)\n', idx_example_PN, ...
    g.cell_metrics.brainRegion{idx_example_PN}, g.cell_metrics.putativeCellType{idx_example_PN});
fprintf('  IN: neuron %d (%s, %s)\n', idx_example_IN, ...
    g.cell_metrics.brainRegion{idx_example_IN}, g.cell_metrics.putativeCellType{idx_example_IN});

%% ===== SETUP FIGURE PARAMETERS FROM PARAMS =====

% Copy parameters to g structure for convenience
g.pre_time = params.pre_time_short;
g.post_time = params.post_time_short;
g.bin_time = params.bin_time;
g.monosyn_window = params.monosyn_window;
g.short_latency_window = params.short_latency_window;
g.smoothvalue = params.smoothvalue;
g.zscore_baseline_time = params.zscore_baseline_time;
g.pre_time_baseline = params.pre_time_baseline;
g.change_threshold = params.change_threshold;
g.pre_time_long = params.pre_time_long;
g.post_time_long = params.post_time_long;

%% ===== CALCULATE PSTH DATA FOR ALL NEURONS (DO THIS ONCE) =====

fprintf('\nCalculating PSTH data for all neurons...\n');
psth_data = struct();
psth_baseline_data = struct();

for tt = 1:length(ttl)

    % Get PSTH and spike data in one call (for raster plots)
    [psth_spx, preAP_norm, postAP_norm] = BAfc_psth_spx('cell_metrics', g.cell_metrics, ...
        'ttl', ttl{tt}, ...
        'pre_time', g.pre_time, ...
        'post_time', g.post_time, ...
        'bin_time', g.bin_time);

    % Store raw PSTH for firing rate calculations
    psth_data(tt).psth_spx_raw = psth_spx;
    psth_data(tt).preAP_norm = preAP_norm;
    psth_data(tt).postAP_norm = postAP_norm;
    psth_data(tt).ttl_name = ttl{tt};

    % Get PSTH with longer baseline for z-score calculation
    [psth_spx_baseline, ~, ~] = BAfc_psth_spx('cell_metrics', g.cell_metrics, ...
        'ttl', ttl{tt}, ...
        'pre_time', g.pre_time_baseline, ...
        'post_time', g.post_time, ...
        'bin_time', g.bin_time);

    % Calculate z-scored PSTH first
    baseline_bins_long = round(g.zscore_baseline_time / g.bin_time);
    psth_z = zeros(size(psth_spx));

    for ii = 1:size(psth_spx, 1)
        baseline_psth = psth_spx_baseline(ii, 1:baseline_bins_long);
        baseline_mean = mean(baseline_psth);
        baseline_std = std(baseline_psth);
        if baseline_std == 0
            baseline_std = 1;
        end
        psth_z(ii, :) = (psth_spx(ii, :) - baseline_mean) / baseline_std;
    end

    % Detect stimulus artifact window (for shock conditions)
    baseline_bins_short = round(g.pre_time / g.bin_time);
    artifact_start_bin_short = baseline_bins_short - round(params.artifact_pre / g.bin_time);
    artifact_end_bin_short = baseline_bins_short + round(params.artifact_post / g.bin_time);

    % Check if this is a shock condition
    is_shock_condition = contains(lower(ttl{tt}), 'shock');

    % Smooth z-scored PSTH if requested
    if g.smoothvalue > 0
        psth_z_smooth = zeros(size(psth_z));
        for ii = 1:size(psth_z, 1)
            % If shock condition, handle artifact window specially
            if is_shock_condition
                % Interpolate artifact window before smoothing
                psth_with_interp = psth_z(ii, :);

                % Linear interpolation across artifact window
                if artifact_start_bin_short > 0 && artifact_end_bin_short <= size(psth_z, 2)
                    value_before = psth_z(ii, artifact_start_bin_short);
                    value_after = psth_z(ii, artifact_end_bin_short);
                    n_artifact_bins = artifact_end_bin_short - artifact_start_bin_short + 1;
                    interp_values = linspace(value_before, value_after, n_artifact_bins);
                    psth_with_interp(artifact_start_bin_short:artifact_end_bin_short) = interp_values;
                end

                % Smooth the interpolated z-scored PSTH
                psth_z_smooth(ii, :) = smoothdata(psth_with_interp, 'sgolay', g.smoothvalue);

                % Restore zeros in artifact window after smoothing
                psth_z_smooth(ii, artifact_start_bin_short:artifact_end_bin_short) = 0;
            else
                % No artifact, smooth normally
                psth_z_smooth(ii, :) = smoothdata(psth_z(ii, :), 'sgolay', g.smoothvalue);
            end
        end
        psth_z = psth_z_smooth;
    end

    % Store z-scored (and smoothed) PSTH
    psth_data(tt).psth_z = psth_z;

    % Store baseline data for later use
    psth_baseline_data(tt).psth_spx = psth_spx_baseline;
    psth_baseline_data(tt).ttl_name = ttl{tt};
end

%% ===== CALCULATE LONG WINDOW DATA FOR IN LIGHT INHIBITION VERIFICATION =====

fprintf('\nCalculating long window data for light inhibition verification...\n');

% Get long window data for light condition only (US light = ttl{4})
[psth_spx_long, preAP_norm_long, postAP_norm_long] = BAfc_psth_spx('cell_metrics', g.cell_metrics, ...
    'ttl', ttl{4}, ...  % US light condition (triptest_shocks_only_light)
    'pre_time', g.pre_time_long, ...
    'post_time', g.post_time_long, ...
    'bin_time', g.bin_time);

% Artifact handling for long window (shock condition)
baseline_bins_long_window = round(g.pre_time_long / g.bin_time);
artifact_start_bin_long_window = baseline_bins_long_window - round(params.artifact_pre / g.bin_time);
artifact_end_bin_long_window = baseline_bins_long_window + round(params.artifact_post / g.bin_time);

if g.smoothvalue > 0
    % Handle artifact window for this specific neuron
    psth_with_interp = psth_spx_long(idx_example_IN, :);

    % Linear interpolation across artifact window
    if artifact_start_bin_long_window > 0 && artifact_end_bin_long_window <= size(psth_spx_long, 2)
        value_before = psth_spx_long(idx_example_IN, artifact_start_bin_long_window);
        value_after = psth_spx_long(idx_example_IN, artifact_end_bin_long_window);
        n_artifact_bins = artifact_end_bin_long_window - artifact_start_bin_long_window + 1;
        interp_values = linspace(value_before, value_after, n_artifact_bins);
        psth_with_interp(artifact_start_bin_long_window:artifact_end_bin_long_window) = interp_values;
    end

    % Smooth the interpolated PSTH
    psth_spx_long(idx_example_IN, :) = smoothdata(psth_with_interp, 'sgolay', g.smoothvalue);

    % Restore zeros in artifact window after smoothing
    psth_spx_long(idx_example_IN, artifact_start_bin_long_window:artifact_end_bin_long_window) = 0;
end

% Store for IN example
IN_long_data.psth_spx = psth_spx_long;
IN_long_data.preAP_norm = preAP_norm_long;
IN_long_data.postAP_norm = postAP_norm_long;

%% ===== SHORT LATENCY RESPONSE ANALYSIS FOR CS AND US RESPONSIVE NEURONS =====

fprintf('\nAnalyzing short latency responses...\n');

% Initialize structure for short latency results
short_latency_results = struct();
group_names = {'CS', 'US'};
group_indices = {idx_CS_responsive, idx_US_responsive};

% TTL mapping: CS uses indices 1,2 (sound_only, sound_only_light)
%              US uses indices 3,4 (shocks_only, shocks_only_light)
ttl_map = struct();
ttl_map.CS = [1, 2];   % triptest_sound_only, triptest_sound_only_light
ttl_map.US = [3, 4];   % triptest_shocks_only, triptest_shocks_only_light

for grp = 1:2
    grp_name = group_names{grp};
    neuron_indices = group_indices{grp};
    ttl_idx_nolight = ttl_map.(grp_name)(1);
    ttl_idx_light = ttl_map.(grp_name)(2);

    % Initialize arrays
    n_neurons = length(neuron_indices);
    short_latency_results.(grp_name).neuron_idx = neuron_indices;
    short_latency_results.(grp_name).cellType = cell(1, n_neurons);
    short_latency_results.(grp_name).brainRegion = cell(1, n_neurons);
    short_latency_results.(grp_name).animal = cell(1, n_neurons);

    % No-light metrics
    short_latency_results.(grp_name).peak_zscore_nolight = zeros(1, n_neurons);
    short_latency_results.(grp_name).mean_zscore_nolight = zeros(1, n_neurons);
    short_latency_results.(grp_name).response_prob_nolight = zeros(1, n_neurons);
    short_latency_results.(grp_name).latency_nolight = zeros(1, n_neurons);

    % Light metrics
    short_latency_results.(grp_name).peak_zscore_light = zeros(1, n_neurons);
    short_latency_results.(grp_name).mean_zscore_light = zeros(1, n_neurons);
    short_latency_results.(grp_name).response_prob_light = zeros(1, n_neurons);
    short_latency_results.(grp_name).latency_light = zeros(1, n_neurons);

    % Process each neuron
    for idx = 1:n_neurons
        neuron_idx = neuron_indices(idx);

        % Store metadata
        short_latency_results.(grp_name).cellType{idx} = g.cell_metrics.putativeCellType{neuron_idx};
        short_latency_results.(grp_name).brainRegion{idx} = g.cell_metrics.brainRegion{neuron_idx};
        short_latency_results.(grp_name).animal{idx} = g.cell_metrics.animal{neuron_idx};

        % Analyze both conditions (no-light, light)
        for cond = 1:2
            if cond == 1
                tt = ttl_idx_nolight;
            else
                tt = ttl_idx_light;
            end

            % Get pre-calculated z-scored data
            psth_z = psth_data(tt).psth_z(neuron_idx, :);
            postAP_norm = psth_data(tt).postAP_norm;

            % Short latency window analysis (0-25ms post-stimulus)
            short_window_baseline_bins = round(g.pre_time / g.bin_time);
            short_latency_bins = round(g.short_latency_window / g.bin_time);
            window_start = short_window_baseline_bins + 1;
            window_end = short_window_baseline_bins + short_latency_bins;

            % Response magnitude (peak and mean z-score)
            window_psth_z = psth_z(window_start:window_end);
            peak_zscore = max(window_psth_z);
            mean_zscore = mean(window_psth_z);

            % Calculate response probability and latency
            num_trials = length(postAP_norm{neuron_idx});
            responsive_trials = 0;
            all_spike_latencies = [];

            for trial = 1:num_trials
                % Get spikes in short latency window for this trial
                trial_spikes = postAP_norm{neuron_idx}{trial};
                spikes_in_window = trial_spikes(trial_spikes > 0 & trial_spikes <= g.short_latency_window);

                if ~isempty(spikes_in_window)
                    responsive_trials = responsive_trials + 1;
                    % Store latency of first spike
                    all_spike_latencies = [all_spike_latencies; spikes_in_window(1)];
                end
            end

            response_prob = responsive_trials / num_trials;

            % Calculate mean latency (in ms)
            if ~isempty(all_spike_latencies)
                mean_latency_ms = mean(all_spike_latencies) * 1000;
            else
                mean_latency_ms = NaN;
            end

            % Store results
            if cond == 1  % No light
                short_latency_results.(grp_name).peak_zscore_nolight(idx) = peak_zscore;
                short_latency_results.(grp_name).mean_zscore_nolight(idx) = mean_zscore;
                short_latency_results.(grp_name).response_prob_nolight(idx) = response_prob;
                short_latency_results.(grp_name).latency_nolight(idx) = mean_latency_ms;
            else  % With light
                short_latency_results.(grp_name).peak_zscore_light(idx) = peak_zscore;
                short_latency_results.(grp_name).mean_zscore_light(idx) = mean_zscore;
                short_latency_results.(grp_name).response_prob_light(idx) = response_prob;
                short_latency_results.(grp_name).latency_light(idx) = mean_latency_ms;
            end
        end
    end

    % Calculate deltas
    short_latency_results.(grp_name).delta_peak_zscore = ...
        short_latency_results.(grp_name).peak_zscore_light - short_latency_results.(grp_name).peak_zscore_nolight;
    short_latency_results.(grp_name).delta_mean_zscore = ...
        short_latency_results.(grp_name).mean_zscore_light - short_latency_results.(grp_name).mean_zscore_nolight;
    short_latency_results.(grp_name).delta_prob = ...
        short_latency_results.(grp_name).response_prob_light - short_latency_results.(grp_name).response_prob_nolight;
    short_latency_results.(grp_name).delta_latency = ...
        short_latency_results.(grp_name).latency_light - short_latency_results.(grp_name).latency_nolight;
end

% Display summary
fprintf('\n===== SHORT LATENCY RESPONSE ANALYSIS (0-25ms) =====\n');
for grp = 1:2
    grp_name = group_names{grp};

    % Calculate response change categories
    delta_peak_z = short_latency_results.(grp_name).peak_zscore_light - short_latency_results.(grp_name).peak_zscore_nolight;
    n_increased = sum(delta_peak_z > 0);
    n_decreased = sum(delta_peak_z < 0);
    n_nochange = sum(delta_peak_z == 0);
    n_total = length(delta_peak_z);

    fprintf('\n%s Responsive Neurons (n=%d):\n', grp_name, n_total);
    fprintf('  Response Change: Increased=%d (%.1f%%), Decreased=%d (%.1f%%), No change=%d (%.1f%%)\n', ...
        n_increased, n_increased/n_total*100, n_decreased, n_decreased/n_total*100, n_nochange, n_nochange/n_total*100);
    fprintf('  Peak Z-score: No-light = %.2f±%.2f, Light = %.2f±%.2f\n', ...
        mean(short_latency_results.(grp_name).peak_zscore_nolight), std(short_latency_results.(grp_name).peak_zscore_nolight), ...
        mean(short_latency_results.(grp_name).peak_zscore_light), std(short_latency_results.(grp_name).peak_zscore_light));
    fprintf('  Response Prob: No-light = %.2f±%.2f, Light = %.2f±%.2f\n', ...
        mean(short_latency_results.(grp_name).response_prob_nolight), std(short_latency_results.(grp_name).response_prob_nolight), ...
        mean(short_latency_results.(grp_name).response_prob_light), std(short_latency_results.(grp_name).response_prob_light));
    fprintf('  Latency (ms): No-light = %.2f±%.2f, Light = %.2f±%.2f\n', ...
        nanmean(short_latency_results.(grp_name).latency_nolight), nanstd(short_latency_results.(grp_name).latency_nolight), ...
        nanmean(short_latency_results.(grp_name).latency_light), nanstd(short_latency_results.(grp_name).latency_light));
end

%% ===== PREPARE DATA FOR POPULATION ANALYSIS (US RESPONSIVE NEURONS) =====

fprintf('\nPreparing data for population analysis...\n');

% Store results for monosynaptic population analysis (US responsive neurons)
pop_results = struct();
pop_results.neuron_idx = [];
pop_results.cellType = {};
pop_results.brainRegion = {};
pop_results.animal = {};
pop_results.peak_zscore_nolight = [];
pop_results.peak_zscore_light = [];
pop_results.response_prob_nolight = [];
pop_results.response_prob_light = [];
pop_results.latency_nolight = [];
pop_results.latency_light = [];
pop_results.delta_zscore = [];
pop_results.delta_prob = [];

% Process all US responsive neurons
for idx = 1:length(idx_US_responsive)
    neuron_idx = idx_US_responsive(idx);

    % Store metadata
    pop_results.neuron_idx(idx) = neuron_idx;
    pop_results.cellType{idx} = g.cell_metrics.putativeCellType{neuron_idx};
    pop_results.brainRegion{idx} = g.cell_metrics.brainRegion{neuron_idx};
    pop_results.animal{idx} = g.cell_metrics.animal{neuron_idx};

    % Process both conditions (no light, with light)
    % US conditions are at indices 3 and 4
    for cond = 1:2
        if cond == 1
            tt = 3;  % triptest_shocks_only (US no-light)
        else
            tt = 4;  % triptest_shocks_only_light (US light)
        end

        % Get pre-calculated z-scored data
        psth_z = psth_data(tt).psth_z(neuron_idx, :);
        postAP_norm = psth_data(tt).postAP_norm;

        % Monosynaptic window analysis (0-25ms post-stimulus)
        % Note: psth_spx has short pre_time (0.1s), so stimulus onset is at that bin
        short_window_baseline_bins = round(g.pre_time / g.bin_time);
        monosyn_bins = round(g.monosyn_window / g.bin_time);
        monosyn_window_start = short_window_baseline_bins + 1;
        monosyn_window_end = short_window_baseline_bins + monosyn_bins;

        % Peak z-score in monosynaptic window
        monosyn_psth_z = psth_z(monosyn_window_start:monosyn_window_end);
        peak_zscore = max(monosyn_psth_z);

        % Calculate response probability and latency
        num_trials = length(postAP_norm{neuron_idx});
        responsive_trials = 0;
        all_spike_latencies = [];

        for trial = 1:num_trials
            % Get spikes in monosynaptic window for this trial
            trial_spikes = postAP_norm{neuron_idx}{trial};
            spikes_in_window = trial_spikes(trial_spikes > 0 & trial_spikes <= g.monosyn_window);

            if ~isempty(spikes_in_window)
                responsive_trials = responsive_trials + 1;
                % Store latency of first spike
                all_spike_latencies = [all_spike_latencies; spikes_in_window(1)];
            end
        end

        response_prob = responsive_trials / num_trials;

        % Calculate mean latency to first spike in responsive trials
        if ~isempty(all_spike_latencies)
            mean_latency_ms = mean(all_spike_latencies) * 1000;
        else
            mean_latency_ms = NaN;
        end

        % Store results
        if cond == 1  % No light
            pop_results.peak_zscore_nolight(idx) = peak_zscore;
            pop_results.response_prob_nolight(idx) = response_prob;
            pop_results.latency_nolight(idx) = mean_latency_ms;
        else  % With light
            pop_results.peak_zscore_light(idx) = peak_zscore;
            pop_results.response_prob_light(idx) = response_prob;
            pop_results.latency_light(idx) = mean_latency_ms;
        end
    end

    % Calculate changes
    pop_results.delta_zscore(idx) = pop_results.peak_zscore_light(idx) - pop_results.peak_zscore_nolight(idx);
    pop_results.delta_prob(idx) = pop_results.response_prob_light(idx) - pop_results.response_prob_nolight(idx);
end

%% ===== VISUALIZE SHORT LATENCY RESPONSES (FIGURE 2) =====

fprintf('\nGenerating Figure 2: Short latency population responses...\n');

fig_short_latency = figure('Position', [100, 100, 1600, 800]);
fig_short_latency.Color = 'w';

for grp = 1:2
    grp_name = group_names{grp};

    % Get data
    peak_z_nolight = short_latency_results.(grp_name).peak_zscore_nolight;
    peak_z_light = short_latency_results.(grp_name).peak_zscore_light;
    prob_nolight = short_latency_results.(grp_name).response_prob_nolight;
    prob_light = short_latency_results.(grp_name).response_prob_light;
    lat_nolight = short_latency_results.(grp_name).latency_nolight;
    lat_light = short_latency_results.(grp_name).latency_light;

    % Calculate change in peak z-score
    delta_peak_z = peak_z_light - peak_z_nolight;

    % Classify neurons based on change threshold
    percent_change = (peak_z_light - peak_z_nolight) ./ peak_z_nolight;

    idx_increased = percent_change > g.change_threshold;
    idx_decreased = percent_change < -g.change_threshold;
    idx_nochange = abs(percent_change) <= g.change_threshold;

    % Count neurons in each category
    n_increased = sum(idx_increased);
    n_decreased = sum(idx_decreased);
    n_nochange = sum(idx_nochange);
    n_total = length(delta_peak_z);

    % Base row for this group
    base_row = (grp - 1) * 2 + 1;

    % --- Panel A: Peak Z-score comparison (color-coded by change) ---
    subplot(4, 4, (base_row-1)*4 + 1);
    hold on;

    % Define colors for different response categories
    color_increased = [0.8, 0.2, 0.2];  % Red - increased response
    color_decreased = [0.2, 0.2, 0.8];  % Blue - decreased response
    color_nochange = [0.5, 0.5, 0.5];   % Gray - no change

    % Plot individual neurons with colors based on change
    for ii = 1:length(peak_z_nolight)
        if idx_increased(ii)
            plot([1, 2], [peak_z_nolight(ii), peak_z_light(ii)], 'o-', ...
                'Color', color_increased, 'MarkerSize', 6, 'LineWidth', 1.5);
        elseif idx_decreased(ii)
            plot([1, 2], [peak_z_nolight(ii), peak_z_light(ii)], 'o-', ...
                'Color', color_decreased, 'MarkerSize', 6, 'LineWidth', 1.5);
        else
            plot([1, 2], [peak_z_nolight(ii), peak_z_light(ii)], 'o-', ...
                'Color', color_nochange, 'MarkerSize', 6, 'LineWidth', 1);
        end
    end

    % Population mean (offset horizontally to avoid hiding data points)
    errorbar(0.85, mean(peak_z_nolight), std(peak_z_nolight)/sqrt(length(peak_z_nolight)), ...
        'ko', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'k');
    errorbar(2.15, mean(peak_z_light), std(peak_z_light)/sqrt(length(peak_z_light)), ...
        'o', 'Color', [0, 0.4, 0.8], 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', [0, 0.4, 0.8]);

    xlim([0.5, 2.5]);
    xticks([1, 2]);
    xticklabels({'No light', 'PV silencing'});
    ylabel('Peak Z-score (log scale)', 'FontSize', 11);
    title(sprintf('%s: Peak Response Magnitude', grp_name), 'FontSize', 12, 'FontWeight', 'bold');
    box on;
    set(gca, 'TickDir', 'out', 'YScale', 'log');

    % Statistical test (only neurons with change)
    idx_changed = idx_increased | idx_decreased;
    if sum(idx_changed) > 0
        [~, p_val] = ttest(peak_z_nolight(idx_changed), peak_z_light(idx_changed));
        y_lim = ylim;
        text(1.5, exp(log(y_lim(2)) * 0.95), sprintf('p = %.4f', p_val), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10, ...
            'BackgroundColor', 'w', 'EdgeColor', 'k');
    end

    % --- Panel B: Pie chart showing response change distribution ---
    subplot(4, 4, (base_row-1)*4 + 2);

    pie_data = [n_increased, n_decreased, n_nochange];
    pie_colors = [color_increased; color_decreased; color_nochange];

    if sum(pie_data) > 0
        p = pie(pie_data);

        % Set colors for pie slices
        slice_idx = 1;
        for ii = 1:2:length(p)
            if slice_idx <= size(pie_colors, 1)
                p(ii).FaceColor = pie_colors(slice_idx, :);
                slice_idx = slice_idx + 1;
            end
        end

        legend({'Increased', 'Decreased', 'No change'}, 'Location', 'southoutside', 'FontSize', 9);
        title(sprintf('%s: Response Change\n(n=%d, %d, %d)', grp_name, n_increased, n_decreased, n_nochange), ...
            'FontSize', 11, 'FontWeight', 'bold');
    end

    % --- Panel C: Response probability comparison ---
    subplot(4, 4, (base_row-1)*4 + 3);
    hold on;

    % Individual neurons with colors based on change
    for ii = 1:length(prob_nolight)
        if idx_increased(ii)
            plot([1, 2], [prob_nolight(ii), prob_light(ii)], 'o-', ...
                'Color', color_increased, 'MarkerSize', 6, 'LineWidth', 1.5);
        elseif idx_decreased(ii)
            plot([1, 2], [prob_nolight(ii), prob_light(ii)], 'o-', ...
                'Color', color_decreased, 'MarkerSize', 6, 'LineWidth', 1.5);
        else
            plot([1, 2], [prob_nolight(ii), prob_light(ii)], 'o-', ...
                'Color', color_nochange, 'MarkerSize', 6, 'LineWidth', 1);
        end
    end

    % Population mean (offset horizontally to avoid hiding data points)
    errorbar(0.85, mean(prob_nolight), std(prob_nolight)/sqrt(length(prob_nolight)), ...
        'ko', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'k');
    errorbar(2.15, mean(prob_light), std(prob_light)/sqrt(length(prob_light)), ...
        'o', 'Color', [0, 0.4, 0.8], 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', [0, 0.4, 0.8]);

    xlim([0.5, 2.5]);
    ylim([0, 1]);
    xticks([1, 2]);
    xticklabels({'No light', 'PV silencing'});
    ylabel('Response Probability', 'FontSize', 11);
    title(sprintf('%s: Response Probability', grp_name), 'FontSize', 12, 'FontWeight', 'bold');
    box on;
    set(gca, 'TickDir', 'out');

    % Statistical test (only neurons with change)
    if sum(idx_changed) > 0
        [~, p_val] = ttest(prob_nolight(idx_changed), prob_light(idx_changed));
        text(1.5, 0.95, sprintf('p = %.4f', p_val), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10, ...
            'BackgroundColor', 'w', 'EdgeColor', 'k');
    end

    % --- Panel D: Latency comparison ---
    subplot(4, 4, (base_row-1)*4 + 4);
    hold on;

    % Filter out NaN values for plotting
    valid_idx = ~isnan(lat_nolight) & ~isnan(lat_light);

    % Individual neurons with colors based on change
    for ii = 1:length(lat_nolight)
        if valid_idx(ii)
            if idx_increased(ii)
                plot([1, 2], [lat_nolight(ii), lat_light(ii)], 'o-', ...
                    'Color', color_increased, 'MarkerSize', 6, 'LineWidth', 1.5);
            elseif idx_decreased(ii)
                plot([1, 2], [lat_nolight(ii), lat_light(ii)], 'o-', ...
                    'Color', color_decreased, 'MarkerSize', 6, 'LineWidth', 1.5);
            else
                plot([1, 2], [lat_nolight(ii), lat_light(ii)], 'o-', ...
                    'Color', color_nochange, 'MarkerSize', 6, 'LineWidth', 1);
            end
        end
    end

    % Population mean (offset horizontally to avoid hiding data points)
    errorbar(0.85, nanmean(lat_nolight), nanstd(lat_nolight)/sqrt(sum(~isnan(lat_nolight))), ...
        'ko', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'k');
    errorbar(2.15, nanmean(lat_light), nanstd(lat_light)/sqrt(sum(~isnan(lat_light))), ...
        'o', 'Color', [0, 0.4, 0.8], 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', [0, 0.4, 0.8]);

    xlim([0.5, 2.5]);
    xticks([1, 2]);
    xticklabels({'No light', 'PV silencing'});
    ylabel('Latency (ms)', 'FontSize', 11);
    title(sprintf('%s: Response Latency', grp_name), 'FontSize', 12, 'FontWeight', 'bold');
    box on;
    set(gca, 'TickDir', 'out');

    % Statistical test (only neurons with change and valid latency)
    idx_changed_valid = idx_changed & valid_idx;
    if sum(idx_changed_valid) > 0
        [~, p_val] = ttest(lat_nolight(idx_changed_valid), lat_light(idx_changed_valid));
        y_lim = ylim;
        text(1.5, y_lim(2) * 0.95, sprintf('p = %.4f', p_val), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10, ...
            'BackgroundColor', 'w', 'EdgeColor', 'k');
    end
end

% Add main title
sgtitle('Short Latency Responses (0-25ms): Optogenetic PV Silencing Effects', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ===== PREPARE DETAILED DATA FOR EXAMPLE NEURONS =====

fprintf('\nPreparing data for example neurons...\n');

example_neurons = [idx_example_PN, idx_example_IN];
example_labels = {'Example PN', 'Example IN'};
example_data = struct();

for ex = 1:length(example_neurons)
    neuron_idx = example_neurons(ex);

    example_data(ex).neuron_idx = neuron_idx;
    example_data(ex).label = example_labels{ex};
    example_data(ex).cellType = g.cell_metrics.putativeCellType{neuron_idx};
    example_data(ex).brainRegion = g.cell_metrics.brainRegion{neuron_idx};

    % Process both conditions using pre-calculated data (US no-light and US light)
    % US conditions are at indices 3 and 4
    for cond = 1:2
        if cond == 1
            tt = 3;  % triptest_shocks_only (US no-light)
        else
            tt = 4;  % triptest_shocks_only_light (US light)
        end

        % Get pre-calculated data
        psth_spx_raw = psth_data(tt).psth_spx_raw;
        psth_z = psth_data(tt).psth_z(neuron_idx, :);
        preAP_norm = psth_data(tt).preAP_norm;
        postAP_norm = psth_data(tt).postAP_norm;

        % Convert raw PSTH to firing rate (Hz)
        num_trials = length(postAP_norm{neuron_idx});
        psth_hz = (psth_spx_raw(neuron_idx, :) / num_trials) / g.bin_time;

        % Store data (use cond, not tt, so we get cond1 and cond2)
        cond_name = sprintf('cond%d', cond);
        example_data(ex).(cond_name).psth_hz = psth_hz;
        example_data(ex).(cond_name).psth_z = psth_z;
        example_data(ex).(cond_name).preAP_norm = preAP_norm{neuron_idx};
        example_data(ex).(cond_name).postAP_norm = postAP_norm{neuron_idx};
        example_data(ex).(cond_name).num_trials = num_trials;
    end
end

%% ===== CREATE FIGURE 1: EXAMPLE NEURONS =====

fprintf('\nGenerating Figure 1: Example neurons...\n');

% Create main figure
fig1 = figure('Position', [50, 50, 1400, 600]);
fig1.Color = 'w';

% Create time axis
time_axis = -g.pre_time:g.bin_time:(g.post_time - g.bin_time);
time_axis_long = -g.pre_time_long:g.bin_time:(g.post_time_long - g.bin_time);

% Monosynaptic window boundaries
monosyn_window_times = [0, g.monosyn_window];

% Color scheme
color_nolight = [0, 0, 0];  % Black
color_light = [0, 0.4, 0.8];  % Blue

%% EXAMPLE NEURONS (PN and IN) - 2 rows × 5 columns

for ex = 1:2
    % Calculate base column for this neuron
    if ex == 1  % PN
        base_col = 1;
    else  % IN
        base_col = 3;
    end

    % Get data
    data_nolight = example_data(ex).cond1;
    data_light = example_data(ex).cond2;

    % Find population index for this neuron
    pop_idx = find(pop_results.neuron_idx == example_data(ex).neuron_idx);

    % --- RASTER PLOTS (No light vs With light) ---
    for cond = 1:2
        if cond == 1
            data = data_nolight;
            title_suffix = 'US (no light)';
            col_idx = base_col;
        else
            data = data_light;
            title_suffix = 'US + PV silencing';
            col_idx = base_col + 1;
        end

        subplot(2, 5, col_idx);
        hold on;

        % Plot monosynaptic window
        patch([monosyn_window_times(1), monosyn_window_times(2), monosyn_window_times(2), monosyn_window_times(1)], ...
            [0, 0, data.num_trials+1, data.num_trials+1], ...
            [1, 0.85, 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

        % Plot raster
        for trial = 1:data.num_trials
            % Pre-stimulus spikes
            if ~isempty(data.preAP_norm{trial})
                spike_times = data.preAP_norm{trial};
                plot(spike_times, trial * ones(size(spike_times)), 'k.', 'MarkerSize', 8);
            end

            % Post-stimulus spikes
            if ~isempty(data.postAP_norm{trial})
                spike_times = data.postAP_norm{trial};
                plot(spike_times, trial * ones(size(spike_times)), 'k.', 'MarkerSize', 8);
            end
        end

        % Plot stimulus onset
        plot([0, 0], [0, data.num_trials+1], 'r-', 'LineWidth', 1.5);

        xlim([-g.pre_time, g.post_time]);
        ylim([0, data.num_trials+1]);
        xlabel('Time from US onset (s)', 'FontSize', 10);
        ylabel('Trial', 'FontSize', 10);
        title(sprintf('%s\n%s', title_suffix, example_data(ex).label), 'FontSize', 11, 'FontWeight', 'bold');
        box on;
        set(gca, 'TickDir', 'out');
    end

    % --- PSTH OVERLAY (bottom) ---
    if ex == 1
        subplot(2, 5, [6, 7]);  % PN PSTH spans columns 1-2, row 2
    else
        subplot(2, 5, [8, 9]);  % IN PSTH spans columns 3-4, row 2
    end
    hold on;

    % Plot monosynaptic window
    y_max = max([max(data_nolight.psth_hz), max(data_light.psth_hz)]) * 1.2;
    patch([monosyn_window_times(1), monosyn_window_times(2), monosyn_window_times(2), monosyn_window_times(1)], ...
        [0, 0, y_max, y_max], ...
        [1, 0.85, 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % Plot both PSTHs
    plot(time_axis, data_nolight.psth_hz, '-', 'Color', color_nolight, 'LineWidth', 2, 'DisplayName', 'No light');
    plot(time_axis, data_light.psth_hz, '-', 'Color', color_light, 'LineWidth', 2, 'DisplayName', 'PV silencing');

    % Plot stimulus onset
    plot([0, 0], [0, y_max], 'r-', 'LineWidth', 1.5);

    xlim([-g.pre_time, g.post_time]);
    ylim([0, y_max]);
    xlabel('Time from US onset (s)', 'FontSize', 10);
    ylabel('Firing rate (Hz)', 'FontSize', 10);
    title(sprintf('%s (neuron %d, %s)', example_data(ex).label, example_data(ex).neuron_idx, example_data(ex).brainRegion), ...
        'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    box on;
    set(gca, 'TickDir', 'out');

    % Add statistics text
    if ~isempty(pop_idx)
        text(0.02, 0.98, sprintf('\\DeltaZ = %.1f\n\\DeltaProb = %.2f', ...
            pop_results.delta_zscore(pop_idx), pop_results.delta_prob(pop_idx)), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9, ...
            'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 2);
    end

    % --- ADD LIGHT INHIBITION RASTER FOR IN (column 5) ---
    if ex == 2  % Only for IN
        % Long window raster for light condition
        subplot(2, 5, 5);  % Column 5, row 1
        hold on;

        data_IN_long = IN_long_data.postAP_norm{idx_example_IN};
        num_trials_long = length(data_IN_long);

        % Plot raster
        for trial = 1:num_trials_long
            % Pre-stimulus spikes
            if ~isempty(IN_long_data.preAP_norm{idx_example_IN}{trial})
                spike_times = IN_long_data.preAP_norm{idx_example_IN}{trial};
                plot(spike_times, trial * ones(size(spike_times)), 'k.', 'MarkerSize', 4);
            end

            % Post-stimulus spikes
            if ~isempty(data_IN_long{trial})
                spike_times = data_IN_long{trial};
                plot(spike_times, trial * ones(size(spike_times)), 'k.', 'MarkerSize', 4);
            end
        end

        % Plot stimulus onset
        plot([0, 0], [0, num_trials_long+1], 'r-', 'LineWidth', 1.5);

        xlim([-g.pre_time_long, g.post_time_long]);
        ylim([0, num_trials_long+1]);
        xlabel('Time from US onset (s)', 'FontSize', 10);
        ylabel('Trial', 'FontSize', 10);
        title(sprintf('IN Light Inhibition\n(-2 to +2s)'), 'FontSize', 11, 'FontWeight', 'bold');
        box on;
        set(gca, 'TickDir', 'out');
    end
end

% Add main title
sgtitle('Optogenetic PV Silencing Effects on Monosynaptic US Responses', 'FontSize', 14, 'FontWeight', 'bold');

%% ===== CREATE FIGURE 3: RESPONSE CHANGE VISUALIZATION =====

fprintf('\nGenerating Figure 3: Response change visualization...\n');

fig3 = figure('Position', [150, 100, 1400, 600]);
fig3.Color = 'w';

% Define colors
color_increased = [0.8, 0.2, 0.2];  % Red
color_decreased = [0.2, 0.2, 0.8];  % Blue
color_nochange = [0.5, 0.5, 0.5];   % Gray

for grp = 1:2
    grp_name = group_names{grp};

    % Get data
    peak_z_nolight = short_latency_results.(grp_name).peak_zscore_nolight;
    peak_z_light = short_latency_results.(grp_name).peak_zscore_light;

    % Calculate percent change
    percent_change = (peak_z_light - peak_z_nolight) ./ peak_z_nolight * 100;

    % Classify neurons
    idx_increased = percent_change > g.change_threshold * 100;
    idx_decreased = percent_change < -g.change_threshold * 100;
    idx_nochange = abs(percent_change) <= g.change_threshold * 100;

    n_increased = sum(idx_increased);
    n_decreased = sum(idx_decreased);
    n_nochange = sum(idx_nochange);
    n_total = length(percent_change);

    % Base row
    base_row = grp;

    % --- Column 1: Paired scatter plot (no-light vs light) ---
    subplot(2, 3, (base_row-1)*3 + 1);
    hold on;

    % Unity line
    max_val = max([peak_z_nolight; peak_z_light]);
    plot([0, max_val], [0, max_val], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Unity');

    % Plot individual neurons color-coded by change direction
    for ii = 1:length(peak_z_nolight)
        if idx_increased(ii)
            scatter(peak_z_nolight(ii), peak_z_light(ii), 80, color_increased, 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        elseif idx_decreased(ii)
            scatter(peak_z_nolight(ii), peak_z_light(ii), 80, color_decreased, 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        else
            scatter(peak_z_nolight(ii), peak_z_light(ii), 60, color_nochange, 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
    end

    xlabel('Peak Z-score (No light)', 'FontSize', 11);
    ylabel('Peak Z-score (PV silencing)', 'FontSize', 11);
    title(sprintf('%s: Paired Response Comparison\n(n=%d)', grp_name, n_total), ...
        'FontSize', 12, 'FontWeight', 'bold');
    axis square;
    box on;
    set(gca, 'TickDir', 'out');

    % Add legend
    h1 = scatter(NaN, NaN, 80, color_increased, 'filled', 'MarkerEdgeColor', 'k');
    h2 = scatter(NaN, NaN, 80, color_decreased, 'filled', 'MarkerEdgeColor', 'k');
    h3 = scatter(NaN, NaN, 60, color_nochange, 'filled', 'MarkerEdgeColor', 'k');
    legend([h1, h2, h3], {'Increased', 'Decreased', 'No change'}, ...
        'Location', 'northwest', 'FontSize', 9);

    % --- Column 2: Bar chart showing proportions ---
    subplot(2, 3, (base_row-1)*3 + 2);

    bar_data = [n_increased, n_decreased, n_nochange];
    bar_colors = [color_increased; color_decreased; color_nochange];

    b = bar(bar_data, 'FaceColor', 'flat');
    b.CData = bar_colors;

    % Add percentage labels on bars
    for ii = 1:3
        text(ii, bar_data(ii) + 0.5, sprintf('%.1f%%', bar_data(ii)/n_total*100), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 10, 'FontWeight', 'bold');
    end

    xticks(1:3);
    xticklabels({'Increased', 'Decreased', 'No change'});
    xtickangle(45);
    ylabel('Number of Neurons', 'FontSize', 11);
    title(sprintf('%s: Response Change Distribution', grp_name), ...
        'FontSize', 12, 'FontWeight', 'bold');
    ylim([0, max(bar_data) * 1.15]);
    box on;
    set(gca, 'TickDir', 'out');

    % --- Column 3: Magnitude of change (violin/distribution plot) ---
    subplot(2, 3, (base_row-1)*3 + 3);
    hold on;

    % Separate percent changes by category
    pct_increased = percent_change(idx_increased);
    pct_decreased = percent_change(idx_decreased);
    pct_nochange = percent_change(idx_nochange);

    % Create violin-like distribution using histogram
    if n_increased > 0
        % Increased group
        [counts_inc, edges_inc] = histcounts(pct_increased, 10);
        centers_inc = (edges_inc(1:end-1) + edges_inc(2:end)) / 2;
        % Normalize to width of 0.3
        counts_inc_norm = counts_inc / max(counts_inc) * 0.3;

        for ii = 1:length(centers_inc)
            patch([1 - counts_inc_norm(ii), 1 + counts_inc_norm(ii), ...
                   1 + counts_inc_norm(ii), 1 - counts_inc_norm(ii)], ...
                  [centers_inc(ii), centers_inc(ii), centers_inc(ii), centers_inc(ii)], ...
                  color_increased, 'EdgeColor', 'k', 'LineWidth', 1);
        end

        % Add median line
        plot([0.7, 1.3], [median(pct_increased), median(pct_increased)], ...
            'k-', 'LineWidth', 2);
    end

    if n_decreased > 0
        % Decreased group
        [counts_dec, edges_dec] = histcounts(pct_decreased, 10);
        centers_dec = (edges_dec(1:end-1) + edges_dec(2:end)) / 2;
        counts_dec_norm = counts_dec / max(counts_dec) * 0.3;

        for ii = 1:length(centers_dec)
            patch([2 - counts_dec_norm(ii), 2 + counts_dec_norm(ii), ...
                   2 + counts_dec_norm(ii), 2 - counts_dec_norm(ii)], ...
                  [centers_dec(ii), centers_dec(ii), centers_dec(ii), centers_dec(ii)], ...
                  color_decreased, 'EdgeColor', 'k', 'LineWidth', 1);
        end

        plot([1.7, 2.3], [median(pct_decreased), median(pct_decreased)], ...
            'k-', 'LineWidth', 2);
    end

    if n_nochange > 0
        % No change group
        [counts_nc, edges_nc] = histcounts(pct_nochange, 10);
        centers_nc = (edges_nc(1:end-1) + edges_nc(2:end)) / 2;
        counts_nc_norm = counts_nc / max(counts_nc) * 0.3;

        for ii = 1:length(centers_nc)
            patch([3 - counts_nc_norm(ii), 3 + counts_nc_norm(ii), ...
                   3 + counts_nc_norm(ii), 3 - counts_nc_norm(ii)], ...
                  [centers_nc(ii), centers_nc(ii), centers_nc(ii), centers_nc(ii)], ...
                  color_nochange, 'EdgeColor', 'k', 'LineWidth', 1);
        end

        plot([2.7, 3.3], [median(pct_nochange), median(pct_nochange)], ...
            'k-', 'LineWidth', 2);
    end

    % Add threshold lines
    plot([0.5, 3.5], [g.change_threshold * 100, g.change_threshold * 100], ...
        'r--', 'LineWidth', 1.5);
    plot([0.5, 3.5], [-g.change_threshold * 100, -g.change_threshold * 100], ...
        'r--', 'LineWidth', 1.5);

    xlim([0.5, 3.5]);
    xticks(1:3);
    xticklabels({'Increased', 'Decreased', 'No change'});
    xtickangle(45);
    ylabel('Percent Change (%)', 'FontSize', 11);
    title(sprintf('%s: Magnitude of Change', grp_name), ...
        'FontSize', 12, 'FontWeight', 'bold');
    box on;
    set(gca, 'TickDir', 'out');

    % Add statistics summary
    text(0.02, 0.98, sprintf('Median:\nInc: %.1f%%\nDec: %.1f%%', ...
        median(pct_increased), median(pct_decreased)), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9, ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 2);
end

% Add main title
sgtitle('Response Change Analysis: PV Silencing Effects on Short Latency Responses (0-25ms)', ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('\n========================================\n');
fprintf('Analysis complete!\n');
fprintf('Results stored in workspace variables:\n');
fprintf('  - responsive_results (from BAfc_identify_responsive_neurons)\n');
fprintf('  - short_latency_results\n');
fprintf('  - pop_results\n');
fprintf('  - example_data\n');
fprintf('========================================\n\n');
