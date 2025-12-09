function results = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl_list, params)
% BAfc_identify_responsive_neurons_v2 - Identifies monosynaptically responsive neurons
%
% FIXED VERSION: Properly handles all brain regions including Astria, CeA, etc.
%
% SYNTAX:
%   results = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl_list)
%   results = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl_list, params)
%
% INPUTS:
%   cell_metrics - Cell metrics structure from BAfc_load_neurons
%   ttl_list     - Cell array of TTL event types to analyze
%   params       - (Optional) Structure with analysis parameters. If not provided,
%                  defaults from BAfc_monosyn_params() are used.
%
% OUTPUTS:
%   results - Structure containing:
%            .params            - Copy of parameters used
%            .ttl_types         - Copy of TTL list
%            .num_analyzed      - Number of neurons analyzed
%            .analyzed_regions  - List of brain regions included
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
% CHANGES FROM ORIGINAL:
%   - Fixed brain region filtering to support all regions (not just LA/BA)
%   - Uses BAfc_monosyn_params() instead of BAfc_optogenetics_params()
%   - Cleaner region handling logic
%
% EXAMPLE:
%   % Use default parameters
%   params = BAfc_monosyn_params();
%   results = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl, params);
%
%   % Modify specific parameters
%   params = BAfc_monosyn_params();
%   params.zscore_threshold = 8;
%   params.brain_regions = {'LA', 'Astria'};
%   results = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl, params);
%
% SEE ALSO: BAfc_monosyn_params
%
% AUTHOR: Claude Code
% DATE: 2025-10-28

%% Validate inputs
if ~isstruct(cell_metrics)
    error('cell_metrics must be a structure');
end

if ~iscell(ttl_list)
    error('ttl_list must be a cell array');
end

% Use default parameters if not provided
if nargin < 3 || isempty(params)
    params = BAfc_monosyn_params();
    fprintf('\n===== USING DEFAULT PARAMETERS =====\n');
else
    % Fill in any missing fields with defaults
    default_params = BAfc_monosyn_params();
    param_fields = fieldnames(default_params);
    for i = 1:length(param_fields)
        if ~isfield(params, param_fields{i})
            params.(param_fields{i}) = default_params.(param_fields{i});
            fprintf('Using default value for %s\n', param_fields{i});
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

% Filter for brain region - FIXED VERSION
if ischar(params.brain_regions) && strcmp(params.brain_regions, 'All')
    % Get all unique brain regions present in data
    unique_regions = unique(cell_metrics.brainRegion);
    fprintf('Brain regions found in data: %s\n', strjoin(unique_regions, ', '));

    % Include all neurons with valid brain region labels
    idx_region = false(size(cell_metrics.brainRegion));
    for i = 1:length(unique_regions)
        idx_region = idx_region | strcmp(cell_metrics.brainRegion, unique_regions{i});
    end

    analyzed_regions = unique_regions;

elseif iscell(params.brain_regions)
    % Specific list of regions provided
    idx_region = false(size(cell_metrics.brainRegion));
    for i = 1:length(params.brain_regions)
        idx_region = idx_region | strcmp(cell_metrics.brainRegion, params.brain_regions{i});
    end

    analyzed_regions = params.brain_regions;
    fprintf('Analyzing specific brain regions: %s\n', strjoin(analyzed_regions, ', '));

else
    % Single region as string
    idx_region = strcmp(cell_metrics.brainRegion, params.brain_regions);
    analyzed_regions = {params.brain_regions};
    fprintf('Analyzing brain region: %s\n', params.brain_regions);
end

% Combine filters
idx_selected = idx_cellType & idx_region;
selected_indices = find(idx_selected);
num_selected = sum(idx_selected);

fprintf('Found %d %s neurons in specified regions out of %d total neurons\n', ...
    num_selected, params.cellType, numel(cell_metrics.cellID));

%% Initialize results structure

results = struct();
results.params = params;
results.ttl_types = ttl_list;
results.num_analyzed = num_selected;
results.analyzed_regions = analyzed_regions;

%% Detect light-inhibited neurons (pre-stimulus inhibition during light trials)

if params.detect_light_inhibited
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
    light_inhibited_zscores = nan(length(selected_indices), 1);
    light_inhibited_pvalues = nan(length(selected_indices), 1);

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

            % Define pre-stimulus windows for comparison
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
                fr_recent = mean(psth_spx_light(ii, prestim_recent_bins)) / params.bin_time;
                fr_baseline = mean(psth_spx_light(ii, prestim_baseline_bins)) / params.bin_time;

                % Detect inhibition based on method
                is_inhibited = false;

                if strcmp(params.light_inhib_method, 'ttest')
                    % T-test method
                    [~, p_val] = ttest2(psth_spx_light(ii, prestim_recent_bins), ...
                                        psth_spx_light(ii, prestim_baseline_bins), ...
                                        'Tail', 'left');

                    light_inhibited_pvalues(idx) = p_val;

                    % Calculate z-score for reference
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
                    % Z-score method
                    baseline_mean = mean(psth_spx_light(ii, prestim_baseline_bins));
                    baseline_std = std(psth_spx_light(ii, prestim_baseline_bins));

                    if baseline_std > 0
                        recent_mean = mean(psth_spx_light(ii, prestim_recent_bins));
                        z_score = (recent_mean - baseline_mean) / baseline_std;
                        light_inhibited_zscores(idx) = z_score;

                        if z_score < params.light_inhib_zscore_threshold
                            is_inhibited = true;
                        end
                    end
                else
                    error('Invalid light_inhib_method: %s', params.light_inhib_method);
                end

                if is_inhibited
                    light_inhibited_neurons = [light_inhibited_neurons, ii];
                end
            end
        end
    end

    % Get unique light-inhibited neurons
    light_inhibited_neurons = unique(light_inhibited_neurons);

    fprintf('Found %d light-inhibited neurons (%.1f%% of %d analyzed)\n', ...
        length(light_inhibited_neurons), ...
        length(light_inhibited_neurons)/num_selected*100, ...
        num_selected);

    % Store in results
    results.light_inhibited_neurons = light_inhibited_neurons;
    results.light_inhibited_zscores = light_inhibited_zscores;
    results.light_inhibited_pvalues = light_inhibited_pvalues;
    results.light_inhibited_neuron_indices = selected_indices;

    % Create exclusion mask
    if params.exclude_light_inhibited
        exclude_light_inhibited_mask = ismember(selected_indices, light_inhibited_neurons);
        fprintf('Light-inhibited neurons will be EXCLUDED from responsive detection.\n');
    else
        exclude_light_inhibited_mask = false(size(selected_indices));
        fprintf('Light-inhibited neurons will NOT be excluded.\n');
    end
else
    % No light-inhibited neuron detection
    results.light_inhibited_neurons = [];
    results.light_inhibited_zscores = [];
    results.light_inhibited_pvalues = [];
    results.light_inhibited_neuron_indices = [];
    exclude_light_inhibited_mask = false(length(selected_indices), 1);
end

%% Analyze each TTL type

for tt = 1:length(ttl_list)
    fprintf('\n========================================\n');
    fprintf('Processing TTL type: %s\n', ttl_list{tt});
    fprintf('========================================\n');

    % Get PSTH data
    [psth_spx, ~, postAP_norm] = BAfc_psth_spx('cell_metrics', cell_metrics, ...
        'ttl', ttl_list{tt}, ...
        'pre_time', params.pre_time, ...
        'post_time', params.post_time, ...
        'bin_time', params.bin_time);

    % Calculate baseline and z-score
    baseline_bins = round(params.pre_time / params.bin_time);
    baseline_psth = psth_spx(:, 1:baseline_bins);

    psth_z = zeros(size(psth_spx));
    for ii = 1:size(psth_spx, 1)
        baseline_mean = mean(baseline_psth(ii, :));
        baseline_std = std(baseline_psth(ii, :));
        if baseline_std == 0
            baseline_std = 1;
        end
        psth_z(ii, :) = (psth_spx(ii, :) - baseline_mean) / baseline_std;
    end

    % Detect artifact window for shock stimuli
    artifact_start_bin = baseline_bins - round(params.artifact_pre / params.bin_time);
    artifact_end_bin = baseline_bins + round(params.artifact_post / params.bin_time);
    is_shock_condition = contains(lower(ttl_list{tt}), 'shock');

    % Smooth z-scored PSTH
    if params.smoothvalue > 0
        psth_z_smooth = zeros(size(psth_z));
        for ii = 1:size(psth_z, 1)
            if is_shock_condition
                % Interpolate artifact window before smoothing
                psth_with_interp = psth_z(ii, :);
                if artifact_start_bin > 0 && artifact_end_bin <= size(psth_z, 2)
                    value_before = psth_z(ii, artifact_start_bin);
                    value_after = psth_z(ii, artifact_end_bin);
                    n_artifact_bins = artifact_end_bin - artifact_start_bin + 1;
                    interp_values = linspace(value_before, value_after, n_artifact_bins);
                    psth_with_interp(artifact_start_bin:artifact_end_bin) = interp_values;
                end
                psth_z_smooth(ii, :) = smoothdata(psth_with_interp, 'sgolay', params.smoothvalue);
                psth_z_smooth(ii, artifact_start_bin:artifact_end_bin) = 0;
            else
                psth_z_smooth(ii, :) = smoothdata(psth_z(ii, :), 'sgolay', params.smoothvalue);
            end
        end
        psth_z = psth_z_smooth;
    end

    % Monosynaptic window
    monosyn_bins = round(params.monosyn_window / params.bin_time);
    monosyn_window_start = baseline_bins + 1;
    monosyn_window_end = baseline_bins + monosyn_bins;

    % Artifact-free window for mean z-score calculation
    if strcmp(params.zscore_method, 'mean') && is_shock_condition
        artifact_bins_in_window = round(params.artifact_post / params.bin_time);
        artifact_free_start = baseline_bins + artifact_bins_in_window + 1;
    else
        artifact_free_start = monosyn_window_start;
    end

    % Initialize storage arrays
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

        all_neuron_idx(idx) = ii;
        all_cellID(idx) = cell_metrics.cellID(ii);
        all_animal{idx} = cell_metrics.animal{ii};
        all_batchID(idx) = cell_metrics.batchIDs(ii);

        num_trials = length(postAP_norm{ii});
        all_num_trials(idx) = num_trials;

        if num_trials == 0
            continue;
        end

        % Calculate z-score (peak or mean)
        if strcmp(params.zscore_method, 'peak')
            monosyn_psth_z = psth_z(ii, monosyn_window_start:monosyn_window_end);
            all_peak_zscore(idx) = max(monosyn_psth_z);
        elseif strcmp(params.zscore_method, 'mean')
            if artifact_free_start <= monosyn_window_end
                monosyn_psth_z = psth_z(ii, artifact_free_start:monosyn_window_end);
                all_peak_zscore(idx) = mean(monosyn_psth_z);
            else
                all_peak_zscore(idx) = NaN;
            end
        else
            error('Invalid zscore_method: %s', params.zscore_method);
        end

        % Calculate response probability
        responsive_trials = 0;
        all_spike_latencies = [];

        for trial = 1:num_trials
            trial_spikes = postAP_norm{ii}{trial};
            spikes_in_window = trial_spikes(trial_spikes > 0 & trial_spikes <= params.monosyn_window);

            if ~isempty(spikes_in_window)
                responsive_trials = responsive_trials + 1;
                all_spike_latencies = [all_spike_latencies; spikes_in_window(1)];
            end
        end

        all_response_probability(idx) = responsive_trials / num_trials;

        if ~isempty(all_spike_latencies)
            all_mean_latency(idx) = mean(all_spike_latencies) * 1000; % ms
        end
    end

    % Apply two-rule system
    rule1_idx = all_peak_zscore >= params.zscore_threshold & ...
                all_response_probability >= params.prob_threshold & ...
                ~isnan(all_peak_zscore);

    rule2_idx = all_peak_zscore >= params.zscore_threshold_strict & ...
                all_response_probability >= params.prob_threshold_lenient & ...
                ~isnan(all_peak_zscore);

    combined_responsive_idx = rule1_idx | rule2_idx;
    combined_responsive_idx(exclude_light_inhibited_mask) = false;

    responsive_idx = find(combined_responsive_idx);

    % Report results
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

    % Store results
    ttl_fieldname = matlab.lang.makeValidName(ttl_list{tt});
    results.(ttl_fieldname).all_neuron_idx = all_neuron_idx;
    results.(ttl_fieldname).all_cellID = all_cellID;
    results.(ttl_fieldname).all_animal = all_animal;
    results.(ttl_fieldname).all_batchID = all_batchID;
    results.(ttl_fieldname).all_peak_zscore = all_peak_zscore;
    results.(ttl_fieldname).all_response_probability = all_response_probability;
    results.(ttl_fieldname).all_mean_latency = all_mean_latency;
    results.(ttl_fieldname).all_num_trials = all_num_trials;
    results.(ttl_fieldname).responsive_idx = responsive_idx;
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
