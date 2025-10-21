function results = BAfc_identify_responsive_neurons(cell_metrics, ttl_list, params)
% BAfc_identify_responsive_neurons - Identifies monosynaptically responsive neurons
%
% SYNTAX:
%   results = BAfc_identify_responsive_neurons(cell_metrics, ttl_list)
%   results = BAfc_identify_responsive_neurons(cell_metrics, ttl_list, params)
%
% INPUTS:
%   cell_metrics - Cell metrics structure from BAfc_load_neurons
%   ttl_list     - Cell array of TTL event types to analyze
%   params       - (Optional) Structure with analysis parameters. If not provided or
%                  empty, default parameters from BAfc_optogenetics_params() are used.
%                  To modify specific parameters, call BAfc_optogenetics_params() first,
%                  modify desired fields, then pass to this function.
%
% OUTPUTS:
%   results - Structure containing:
%            .params            - Copy of parameters used
%            .ttl_types         - Copy of TTL list
%            .num_analyzed      - Number of neurons analyzed
%            .[ttl_name]        - For each TTL type, contains responsive neuron data
%
% DESCRIPTION:
%   Identifies neurons with monosynaptic responses using a two-rule system:
%
%   Rule 1: Z-score >= zscore_threshold (5) AND probability >= prob_threshold (0.25)
%   Rule 2: Z-score >= zscore_threshold_strict (10) AND probability >= prob_threshold_lenient (0.1)
%
%   A neuron is classified as responsive if it meets EITHER rule.
%
%   Z-score calculation:
%   - If params.zscore_method = 'peak': uses peak z-score in monosynaptic window
%   - If params.zscore_method = 'mean': uses mean z-score in artifact-free period
%     (for shock stimuli, excludes first 12ms after stimulus)
%
%   Additional metric: Mean latency calculated from responsive trials
%
% EXAMPLE:
%   % Use default parameters
%   results = BAfc_identify_responsive_neurons(cell_metrics, ttl);
%
%   % Modify specific parameters
%   params = BAfc_optogenetics_params();
%   params.zscore_threshold = 8;
%   params.bR = 'LA';
%   results = BAfc_identify_responsive_neurons(cell_metrics, ttl, params);
%
% SEE ALSO: BAfc_optogenetics_params
%
% AUTHOR: Claude Code
% DATE: 2025-10-19

%% Validate inputs
if ~isstruct(cell_metrics)
    error('cell_metrics must be a structure');
end

if ~iscell(ttl_list)
    error('ttl_list must be a cell array');
end

% Use default parameters if not provided
if nargin < 3 || isempty(params)
    params = BAfc_optogenetics_params();
    fprintf('\n===== USING DEFAULT PARAMETERS =====\n');
else
    % Fill in any missing fields with defaults
    default_params = BAfc_optogenetics_params();
    param_fields = fieldnames(default_params);
    for i = 1:length(param_fields)
        if ~isfield(params, param_fields{i})
            params.(param_fields{i}) = default_params.(param_fields{i});
            fprintf('Using default value for %s: %s\n', param_fields{i}, ...
                mat2str(default_params.(param_fields{i})));
        end
    end
end

%% Filter neurons by cell type and brain region

fprintf('\n===== FILTERING NEURONS =====\n');

% Filter for cell type
if strcmp(params.cellType, 'All')
    idx_cellType = strcmp(cell_metrics.putativeCellType, 'PN') | ...
                   strcmp(cell_metrics.putativeCellType, 'IN');
else
    idx_cellType = strcmp(cell_metrics.putativeCellType, params.cellType);
end

% Filter for brain region
if strcmp(params.bR, 'All')
    idx_region = strcmp(cell_metrics.brainRegion, 'LA') | ...
                 strcmp(cell_metrics.brainRegion, 'BA');
else
    idx_region = strcmp(cell_metrics.brainRegion, params.bR);
end

% Combine filters
idx_selected = idx_cellType & idx_region;
selected_indices = find(idx_selected);
num_selected = sum(idx_selected);

fprintf('Found %d %s neurons in %s out of %d total neurons\n', ...
    num_selected, params.cellType, params.bR, numel(cell_metrics.cellID));

%% Initialize results structure

results = struct();
results.params = params;
results.ttl_types = ttl_list;
results.num_analyzed = num_selected;

%% Detect light-inhibited neurons (pre-stimulus inhibition during light trials)

fprintf('\n===== DETECTING LIGHT-INHIBITED NEURONS =====\n');
fprintf('Method: %s\n', params.light_inhib_method);

% Find light trials (TTL names containing 'light')
light_trial_idx = [];
for tt = 1:length(ttl_list)
    if contains(lower(ttl_list{tt}), 'light')
        light_trial_idx = [light_trial_idx, tt];
    end
end

% Initialize light-inhibited neurons storage
light_inhibited_neurons = [];
light_inhibited_zscores = nan(length(selected_indices), 1);  % Store z-scores for all neurons
light_inhibited_pvalues = nan(length(selected_indices), 1);  % Store p-values for all neurons

if ~isempty(light_trial_idx)
    % Analyze pre-stimulus inhibition for each light condition
    for tt = light_trial_idx
        fprintf('Checking light inhibition in: %s\n', ttl_list{tt});

        % Get PSTH data for this light condition
        [psth_spx_light, ~, ~] = BAfc_psth_spx('cell_metrics', cell_metrics, ...
            'ttl', ttl_list{tt}, ...
            'pre_time', params.pre_time, ...
            'post_time', params.post_time, ...
            'bin_time', params.bin_time);

        % Define pre-stimulus windows for comparison using parameters
        baseline_bins = round(params.pre_time / params.bin_time);
        recent_window_bins = round(params.light_inhib_window_recent / params.bin_time);
        prestim_recent_bins = baseline_bins - recent_window_bins + 1:baseline_bins;

        baseline_window_bins = round(params.light_inhib_window_baseline / params.bin_time);
        prestim_baseline_start = baseline_bins - recent_window_bins - baseline_window_bins + 1;
        prestim_baseline_end = baseline_bins - recent_window_bins;
        prestim_baseline_bins = prestim_baseline_start:prestim_baseline_end;

        % Test each selected neuron for light inhibition
        for idx = 1:length(selected_indices)
            ii = selected_indices(idx);

            % Get firing rates in both windows
            fr_recent = mean(psth_spx_light(ii, prestim_recent_bins)) / params.bin_time;  % Hz
            fr_baseline = mean(psth_spx_light(ii, prestim_baseline_bins)) / params.bin_time;  % Hz

            % Detect inhibition based on method
            is_inhibited = false;

            if strcmp(params.light_inhib_method, 'ttest')
                % T-test method: statistical test for significant decrease
                [~, p_val] = ttest2(psth_spx_light(ii, prestim_recent_bins), ...
                                    psth_spx_light(ii, prestim_baseline_bins), ...
                                    'Tail', 'left');

                % Store p-value
                light_inhibited_pvalues(idx) = p_val;

                % Also calculate z-score for reference
                baseline_mean = mean(psth_spx_light(ii, prestim_baseline_bins));
                baseline_std = std(psth_spx_light(ii, prestim_baseline_bins));
                if baseline_std > 0
                    recent_mean = mean(psth_spx_light(ii, prestim_recent_bins));
                    light_inhibited_zscores(idx) = (recent_mean - baseline_mean) / baseline_std;
                end

                if p_val < params.light_inhib_p_threshold && fr_recent < fr_baseline
                    is_inhibited = true;
                end

            elseif strcmp(params.light_inhib_method, 'zscore')
                % Z-score method: compare recent period to baseline mean/std
                baseline_mean = mean(psth_spx_light(ii, prestim_baseline_bins));
                baseline_std = std(psth_spx_light(ii, prestim_baseline_bins));

                if baseline_std > 0
                    recent_mean = mean(psth_spx_light(ii, prestim_recent_bins));
                    z_score = (recent_mean - baseline_mean) / baseline_std;

                    % Store z-score
                    light_inhibited_zscores(idx) = z_score;

                    if z_score < params.light_inhib_zscore_threshold
                        is_inhibited = true;
                    end
                end
            else
                error('Invalid light_inhib_method: %s. Must be ''ttest'' or ''zscore''.', ...
                    params.light_inhib_method);
            end

            if is_inhibited
                light_inhibited_neurons = [light_inhibited_neurons, ii];
            end
        end
    end
end

% Get unique light-inhibited neurons (in case detected in multiple conditions)
light_inhibited_neurons = unique(light_inhibited_neurons);

fprintf('Found %d light-inhibited neurons (%.1f%% of %d analyzed)\n', ...
    length(light_inhibited_neurons), ...
    length(light_inhibited_neurons)/num_selected*100, ...
    num_selected);

% Store in results
results.light_inhibited_neurons = light_inhibited_neurons;
results.light_inhibited_zscores = light_inhibited_zscores;  % Z-scores for all analyzed neurons
results.light_inhibited_pvalues = light_inhibited_pvalues;  % P-values for all analyzed neurons (if ttest method)
results.light_inhibited_neuron_indices = selected_indices;  % Corresponding neuron indices

% Summary statistics for light-inhibited neurons
if ~isempty(light_inhibited_neurons)
    inhib_idx = ismember(selected_indices, light_inhibited_neurons);
    fprintf('Light-inhibited neuron z-scores: mean=%.2f, median=%.2f, range=[%.2f, %.2f]\n', ...
        nanmean(light_inhibited_zscores(inhib_idx)), ...
        nanmedian(light_inhibited_zscores(inhib_idx)), ...
        min(light_inhibited_zscores(inhib_idx)), ...
        max(light_inhibited_zscores(inhib_idx)));
end

% Create exclusion mask for responsive neuron detection
if params.exclude_light_inhibited
    exclude_light_inhibited_mask = ismember(selected_indices, light_inhibited_neurons);
    fprintf('Light-inhibited neurons will be EXCLUDED from responsive neuron detection.\n');
else
    exclude_light_inhibited_mask = false(size(selected_indices));
    fprintf('Light-inhibited neurons will NOT be excluded (params.exclude_light_inhibited = false).\n');
end

%% Analyze each TTL type

for tt = 1:length(ttl_list)
    fprintf('\n========================================\n');
    fprintf('Processing TTL type: %s\n', ttl_list{tt});
    fprintf('========================================\n');

    % Get PSTH data with fine temporal resolution
    [psth_spx, ~, postAP_norm] = BAfc_psth_spx('cell_metrics', cell_metrics, ...
        'ttl', ttl_list{tt}, ...
        'pre_time', params.pre_time, ...
        'post_time', params.post_time, ...
        'bin_time', params.bin_time);

    % Calculate baseline firing rate and std for each neuron (before z-scoring)
    baseline_bins = round(params.pre_time / params.bin_time);
    baseline_psth = psth_spx(:, 1:baseline_bins);

    % Z-score the PSTH first
    psth_z = zeros(size(psth_spx));
    for ii = 1:size(psth_spx, 1)
        baseline_mean = mean(baseline_psth(ii, :));
        baseline_std = std(baseline_psth(ii, :));
        if baseline_std == 0
            baseline_std = 1;
        end
        psth_z(ii, :) = (psth_spx(ii, :) - baseline_mean) / baseline_std;
    end

    % Detect stimulus artifact window (for shock conditions)
    artifact_start_bin = baseline_bins - round(params.artifact_pre / params.bin_time);
    artifact_end_bin = baseline_bins + round(params.artifact_post / params.bin_time);

    % Check if this is a shock condition (contains 'shock' in TTL name)
    is_shock_condition = contains(lower(ttl_list{tt}), 'shock');

    % Smooth z-scored PSTH if requested
    if params.smoothvalue > 0
        psth_z_smooth = zeros(size(psth_z));
        for ii = 1:size(psth_z, 1)
            % If shock condition, handle artifact window specially
            if is_shock_condition
                % Interpolate artifact window before smoothing
                psth_with_interp = psth_z(ii, :);

                % Linear interpolation across artifact window
                if artifact_start_bin > 0 && artifact_end_bin <= size(psth_z, 2)
                    value_before = psth_z(ii, artifact_start_bin);
                    value_after = psth_z(ii, artifact_end_bin);
                    n_artifact_bins = artifact_end_bin - artifact_start_bin + 1;
                    interp_values = linspace(value_before, value_after, n_artifact_bins);
                    psth_with_interp(artifact_start_bin:artifact_end_bin) = interp_values;
                end

                % Smooth the interpolated z-scored PSTH
                psth_z_smooth(ii, :) = smoothdata(psth_with_interp, 'sgolay', params.smoothvalue);

                % Restore zeros in artifact window after smoothing
                psth_z_smooth(ii, artifact_start_bin:artifact_end_bin) = 0;
            else
                % No artifact, smooth normally
                psth_z_smooth(ii, :) = smoothdata(psth_z(ii, :), 'sgolay', params.smoothvalue);
            end
        end
        psth_z = psth_z_smooth;
    end

    % Monosynaptic window bins
    monosyn_bins = round(params.monosyn_window / params.bin_time);
    monosyn_window_start = baseline_bins + 1;
    monosyn_window_end = baseline_bins + monosyn_bins;

    % For mean z-score calculation with artifact exclusion
    if strcmp(params.zscore_method, 'mean')
        % If shock condition, exclude artifact period from mean calculation
        if is_shock_condition
            artifact_bins_in_window = round(params.artifact_post / params.bin_time);
            % Calculate artifact-free window start (after artifact ends)
            artifact_free_start = baseline_bins + artifact_bins_in_window + 1;
        else
            % No artifact, use full monosynaptic window
            artifact_free_start = monosyn_window_start;
        end
    end

    % Initialize storage arrays for ALL selected neurons
    all_neuron_idx = zeros(length(selected_indices), 1);
    all_cellID = zeros(length(selected_indices), 1);
    all_animal = cell(length(selected_indices), 1);
    all_batchID = zeros(length(selected_indices), 1);
    all_peak_zscore = nan(length(selected_indices), 1);
    all_response_probability = nan(length(selected_indices), 1);
    all_mean_latency = nan(length(selected_indices), 1);
    all_num_trials = zeros(length(selected_indices), 1);

    % Loop through selected neurons
    for idx = 1:length(selected_indices)
        ii = selected_indices(idx);

        % Store neuron identifiers
        all_neuron_idx(idx) = ii;
        all_cellID(idx) = cell_metrics.cellID(ii);
        all_animal{idx} = cell_metrics.animal{ii};
        all_batchID(idx) = cell_metrics.batchIDs(ii);

        % Get number of trials for this neuron
        num_trials = length(postAP_norm{ii});
        all_num_trials(idx) = num_trials;

        if num_trials == 0
            continue;
        end

        % CRITERION 1: Z-score in monosynaptic window (peak or mean)
        if strcmp(params.zscore_method, 'peak')
            % Peak z-score in monosynaptic window
            monosyn_psth_z = psth_z(ii, monosyn_window_start:monosyn_window_end);
            monosyn_peak_zscore = max(monosyn_psth_z);
            all_peak_zscore(idx) = monosyn_peak_zscore;
        elseif strcmp(params.zscore_method, 'mean')
            % Mean z-score in artifact-free portion of monosynaptic window
            if artifact_free_start <= monosyn_window_end
                monosyn_psth_z = psth_z(ii, artifact_free_start:monosyn_window_end);
                monosyn_mean_zscore = mean(monosyn_psth_z);
                all_peak_zscore(idx) = monosyn_mean_zscore;  % Store in same field for consistency
            else
                % Artifact window exceeds monosynaptic window, no valid data
                all_peak_zscore(idx) = NaN;
            end
        else
            error('Invalid zscore_method: %s. Must be ''peak'' or ''mean''.', params.zscore_method);
        end

        % CRITERION 2: Response probability (trial-by-trial analysis)
        responsive_trials = 0;
        all_spike_latencies = [];

        for trial = 1:num_trials
            % Get spikes in monosynaptic window for this trial
            trial_spikes = postAP_norm{ii}{trial};
            spikes_in_window = trial_spikes(trial_spikes > 0 & trial_spikes <= params.monosyn_window);

            if ~isempty(spikes_in_window)
                responsive_trials = responsive_trials + 1;
                % Store latency of first spike
                all_spike_latencies = [all_spike_latencies; spikes_in_window(1)];
            end
        end

        response_prob = responsive_trials / num_trials;
        all_response_probability(idx) = response_prob;

        % CRITERION 3: Mean latency to first spike in responsive trials
        if ~isempty(all_spike_latencies)
            mean_latency = mean(all_spike_latencies) * 1000; % Convert to ms
            all_mean_latency(idx) = mean_latency;
        end
    end

    % Find indices of responsive neurons using two-rule system
    % Rule 1: zscore >= threshold (5) AND prob >= threshold (0.25)
    % Rule 2: zscore >= strict threshold (10) AND prob >= lenient threshold (0.1)
    % Note: all_peak_zscore may contain mean z-scores if zscore_method='mean'

    rule1_idx = all_peak_zscore >= params.zscore_threshold & ...
                all_response_probability >= params.prob_threshold & ...
                ~isnan(all_peak_zscore);

    rule2_idx = all_peak_zscore >= params.zscore_threshold_strict & ...
                all_response_probability >= params.prob_threshold_lenient & ...
                ~isnan(all_peak_zscore);

    % Combine both rules (logical OR)
    combined_responsive_idx = rule1_idx | rule2_idx;

    % EXCLUDE light-inhibited neurons from responsive neurons
    combined_responsive_idx(exclude_light_inhibited_mask) = false;

    responsive_idx = find(combined_responsive_idx);

    % Report which rule(s) detected each responsive neuron
    n_rule1_only = sum(rule1_idx & ~rule2_idx & ~exclude_light_inhibited_mask');
    n_rule2_only = sum(~rule1_idx & rule2_idx & ~exclude_light_inhibited_mask');
    n_both_rules = sum(rule1_idx & rule2_idx & ~exclude_light_inhibited_mask');
    n_excluded = sum((rule1_idx | rule2_idx) & exclude_light_inhibited_mask');

    fprintf('Rule 1 only (z>=%.1f, p>=%.2f): %d neurons\n', ...
        params.zscore_threshold, params.prob_threshold, n_rule1_only);
    fprintf('Rule 2 only (z>=%.1f, p>=%.2f): %d neurons\n', ...
        params.zscore_threshold_strict, params.prob_threshold_lenient, n_rule2_only);
    fprintf('Both rules: %d neurons\n', n_both_rules);
    fprintf('Excluded (light-inhibited): %d neurons\n', n_excluded);

    % Store ALL neuron results for this TTL type
    ttl_fieldname = matlab.lang.makeValidName(ttl_list{tt});
    results.(ttl_fieldname).all_neuron_idx = all_neuron_idx;
    results.(ttl_fieldname).all_cellID = all_cellID;
    results.(ttl_fieldname).all_animal = all_animal;
    results.(ttl_fieldname).all_batchID = all_batchID;
    results.(ttl_fieldname).all_peak_zscore = all_peak_zscore;
    results.(ttl_fieldname).all_response_probability = all_response_probability;
    results.(ttl_fieldname).all_mean_latency = all_mean_latency;
    results.(ttl_fieldname).all_num_trials = all_num_trials;

    % Store indices of RESPONSIVE neurons
    results.(ttl_fieldname).responsive_idx = responsive_idx;

    % Store responsive neuron data (for convenience)
    results.(ttl_fieldname).neuron_idx = all_neuron_idx(responsive_idx);
    results.(ttl_fieldname).cellID = all_cellID(responsive_idx);
    results.(ttl_fieldname).animal = all_animal(responsive_idx);
    results.(ttl_fieldname).batchID = all_batchID(responsive_idx);
    results.(ttl_fieldname).peak_zscore = all_peak_zscore(responsive_idx);
    results.(ttl_fieldname).response_probability = all_response_probability(responsive_idx);
    results.(ttl_fieldname).mean_latency = all_mean_latency(responsive_idx);
    results.(ttl_fieldname).num_trials = all_num_trials(responsive_idx);

    fprintf('Found %d monosynaptically activated neurons (%.1f%% of %d analyzed)\n', ...
        length(responsive_idx), 100*length(responsive_idx)/num_selected, num_selected);
end

fprintf('\n========================================\n');
fprintf('Responsiveness detection complete!\n');
fprintf('========================================\n\n');

end
