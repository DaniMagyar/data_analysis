% BAfc_figure_5_supp_light_inhibited.m
% Supplementary figure: Light-inhibited neurons during optogenetic manipulation
% Identifies and visualizes neurons inhibited by pre-stimulus light
% NEW LAYOUT: Row 1=LA heatmaps, Row 2=LA lineplots, Row 3=Astria heatmaps, Row 4=Astria lineplots

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

% TTL types: No-light and light conditions for CS and US
ttl_light = {'triptest_sound_only_light', 'triptest_shocks_only_light'};
ttl_nolight = {'triptest_sound_only', 'triptest_shocks_only'};
hmptitles = {'CS', 'US'};

% Combine all TTL types for loading
ttl_all = [ttl_nolight, ttl_light];
cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl_all);
cell_metrics = BAfc_match_celltypes('cell_metrics', cell_metrics);

g.colors = BAfc_colors;
g.fontSize1 = 10;
g.fontSize2 = 10;
g.xlinewidth = 2;
g.pre_time = 5;
g.post_time = 4;
g.test_time = 1;

% Clustering Parameters (same as figure 3)
g.alpha = 0.5;
g.excitation_threshold = 2;
g.inhibition_fr_drop = 0.50;
g.use_percentile = true;
g.clim_percentile = 99;
g.onset_threshold = g.excitation_threshold;
g.bin_time = 0.001;
g.min_consec_bins = max(1, round(0.020 / g.bin_time));
g.smoothvalue = 201;
g.plotwin = [2 2];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;

% Light inhibition detection parameters
g.light_inhib_window_recent = 0.5;  % 0-0.5s before stimulus (recent)
g.light_inhib_window_baseline = 4.5;  % 4.5s baseline period
g.light_inhib_p_threshold = 0.05;  % p-value threshold for Wilcoxon
g.light_inhib_fr_drop = 0.5;  % 50% FR drop threshold (same as BAfc_figure_2)

% Brain regions to analyze
brain_regions = {'LA', 'Astria'};

%% Calculate PSTHs for all conditions
psthZ_nolight = cell(1, 2);  % CS, US (no-light)
psthHz_nolight = cell(1, 2);
psthZ_light = cell(1, 2);  % CS, US (light)
psthHz_light = cell(1, 2);

% Calculate Savitzky-Golay filter delay correction
filter_delay = floor(g.smoothvalue / 2);  % Half the filter width (symmetric kernel)

% Process no-light conditions
for hmp = 1:2  % CS, US
    psth_spx_og = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl_nolight{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Convert to Hz
    num_trials = size(cell_metrics.general.(ttl_nolight{hmp}){1}, 1);
    psth_hz = psth_spx_og / (num_trials * g.bin_time);
    psth_hz_smooth = smoothdata(psth_hz, 2, 'sgolay', g.smoothvalue);

    % Correct for filter delay by shifting forward (move data earlier in time)
    psth_hz_corrected = zeros(size(psth_hz_smooth));
    psth_hz_corrected(:, filter_delay+1:end) = psth_hz_smooth(:, 1:end-filter_delay);
    psth_hz_corrected(:, 1:filter_delay) = repmat(psth_hz_smooth(:, 1), 1, filter_delay);
    psthHz_nolight{hmp} = psth_hz_corrected;

    % Z-score using baseline period only
    baseline_idx = 1:(g.pre_time / g.bin_time);
    baseline_mean = mean(psth_spx_og(:, baseline_idx), 2);
    baseline_std = std(psth_spx_og(:, baseline_idx), 0, 2);
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;
    psth_spx_smooth = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);

    % Correct for filter delay by shifting forward (move data earlier in time)
    psth_spx_corrected = zeros(size(psth_spx_smooth));
    psth_spx_corrected(:, filter_delay+1:end) = psth_spx_smooth(:, 1:end-filter_delay);
    psth_spx_corrected(:, 1:filter_delay) = repmat(psth_spx_smooth(:, 1), 1, filter_delay);
    psthZ_nolight{hmp} = psth_spx_corrected;
end

% Process light conditions
for hmp = 1:2  % CS, US
    psth_spx_og = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl_light{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Convert to Hz
    num_trials = size(cell_metrics.general.(ttl_light{hmp}){1}, 1);
    psth_hz = psth_spx_og / (num_trials * g.bin_time);
    psth_hz_smooth = smoothdata(psth_hz, 2, 'sgolay', g.smoothvalue);

    % Correct for filter delay by shifting forward (move data earlier in time)
    psth_hz_corrected = zeros(size(psth_hz_smooth));
    psth_hz_corrected(:, filter_delay+1:end) = psth_hz_smooth(:, 1:end-filter_delay);
    psth_hz_corrected(:, 1:filter_delay) = repmat(psth_hz_smooth(:, 1), 1, filter_delay);
    psthHz_light{hmp} = psth_hz_corrected;

    % Z-score using baseline period only
    baseline_idx = 1:(g.pre_time / g.bin_time);
    baseline_mean = mean(psth_spx_og(:, baseline_idx), 2);
    baseline_std = std(psth_spx_og(:, baseline_idx), 0, 2);
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;
    psth_spx_smooth = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);

    % Correct for filter delay by shifting forward (move data earlier in time)
    psth_spx_corrected = zeros(size(psth_spx_smooth));
    psth_spx_corrected(:, filter_delay+1:end) = psth_spx_smooth(:, 1:end-filter_delay);
    psth_spx_corrected(:, 1:filter_delay) = repmat(psth_spx_smooth(:, 1), 1, filter_delay);
    psthZ_light{hmp} = psth_spx_corrected;
end

%% Identify light-inhibited neurons for each brain region
results_all = cell(1, 2);

for br = 1:2
    % Get ALL neuron indices for this brain region (no cell type filter)
    idx_neurons = strcmp(cell_metrics.brainRegion, brain_regions{br});
    n_neurons = sum(idx_neurons);

    if n_neurons == 0
        continue;
    end

    % Define pre-stimulus windows for comparison
    baseline_bins = round(g.pre_time / g.bin_time);
    recent_window_bins = round(g.light_inhib_window_recent / g.bin_time);
    prestim_recent_bins = baseline_bins - recent_window_bins + 1:baseline_bins;

    baseline_window_bins = round(g.light_inhib_window_baseline / g.bin_time);
    prestim_baseline_start = baseline_bins - recent_window_bins - baseline_window_bins + 1;
    prestim_baseline_end = baseline_bins - recent_window_bins;
    prestim_baseline_bins = prestim_baseline_start:prestim_baseline_end;

    % Initialize light-inhibited detection arrays
    neuron_indices = find(idx_neurons);
    light_inhibited_global = [];  % Store global neuron indices
    light_inhib_pvalues = nan(n_neurons, 1);
    light_inhib_fr_drops = nan(n_neurons, 1);

    % Test each neuron across all light conditions
    % If inhibited in ANY light condition, mark as light-inhibited
    for hmp = 1:numel(ttl_light)
        psth_spx_light = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl_light{hmp}, ...
            'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

        for n = 1:n_neurons
            ii = neuron_indices(n);

            % Get firing rates in both windows
            fr_recent = mean(psth_spx_light(ii, prestim_recent_bins)) / g.bin_time;
            fr_baseline = mean(psth_spx_light(ii, prestim_baseline_bins)) / g.bin_time;

            % Calculate FR drop (same as BAfc_figure_2)
            fr_drop = (fr_baseline - fr_recent) / (fr_baseline + eps);

            % Wilcoxon ranksum test (non-parametric, unpaired)
            [p_val, ~] = ranksum(psth_spx_light(ii, prestim_baseline_bins), ...
                                 psth_spx_light(ii, prestim_recent_bins), ...
                                 'tail', 'right');  % Test if baseline > recent

            % Store p-value if not already stored or if more significant
            if isnan(light_inhib_pvalues(n)) || p_val < light_inhib_pvalues(n)
                light_inhib_pvalues(n) = p_val;
            end

            % Store FR drop if not already stored or if larger
            if isnan(light_inhib_fr_drops(n)) || fr_drop > light_inhib_fr_drops(n)
                light_inhib_fr_drops(n) = fr_drop;
            end

            % Detect inhibition: p<0.05 AND FR drop >= 50%
            is_inhibited = false;
            if p_val < g.light_inhib_p_threshold && fr_drop >= g.light_inhib_fr_drop
                is_inhibited = true;
            end

            if is_inhibited
                light_inhibited_global = [light_inhibited_global, ii];
            end
        end
    end

    % Get unique light-inhibited neurons
    light_inhibited_global = unique(light_inhibited_global);
    n_light_inhibited = length(light_inhibited_global);

    if n_light_inhibited == 0
        continue;
    end

    % Extract PSTHs for light-inhibited neurons only
    % No-light conditions
    psth_CS_nolight = psthZ_nolight{1}(light_inhibited_global, :);
    psth_US_nolight = psthZ_nolight{2}(light_inhibited_global, :);
    psth_CS_Hz_nolight = psthHz_nolight{1}(light_inhibited_global, :);
    psth_US_Hz_nolight = psthHz_nolight{2}(light_inhibited_global, :);

    % Light conditions
    psth_CS_light = psthZ_light{1}(light_inhibited_global, :);
    psth_US_light = psthZ_light{2}(light_inhibited_global, :);
    psth_CS_Hz_light = psthHz_light{1}(light_inhibited_global, :);
    psth_US_Hz_light = psthHz_light{2}(light_inhibited_global, :);

    % Get local indices for light-inhibited neurons
    [~, local_idx] = ismember(light_inhibited_global, neuron_indices);
    light_inhib_fr_drops_inhibited = light_inhib_fr_drops(local_idx);
    light_inhib_pvalues_inhibited = light_inhib_pvalues(local_idx);

    % Calculate responses based on CS and US no-light (same as figure 3)
    CS_peak = max(psth_CS_nolight(:, g.roi), [], 2);
    US_peak = max(psth_US_nolight(:, g.roi), [], 2);

    baseline_idx = 1:(g.pre_time / g.bin_time);
    CS_baseline_fr = mean(psth_CS_Hz_nolight(:, baseline_idx), 2);
    US_baseline_fr = mean(psth_US_Hz_nolight(:, baseline_idx), 2);
    CS_test_fr = mean(psth_CS_Hz_nolight(:, g.roi), 2);
    US_test_fr = mean(psth_US_Hz_nolight(:, g.roi), 2);
    CS_fr_drop = (CS_baseline_fr - CS_test_fr) ./ (CS_baseline_fr + eps);
    US_fr_drop = (US_baseline_fr - US_test_fr) ./ (US_baseline_fr + eps);

    % Classify neurons based on CS and US responses (same as figure 2/3)
    CS_excited = CS_peak >= g.excitation_threshold;
    US_excited = US_peak >= g.excitation_threshold;
    CS_inhibited = CS_fr_drop >= g.inhibition_fr_drop;
    US_inhibited = US_fr_drop >= g.inhibition_fr_drop;

    Clusters = zeros(n_light_inhibited, 1);
    Clusters(CS_excited & ~US_excited) = 1;  % CS-selective
    Clusters(US_excited & ~CS_excited) = 2;  % US-selective
    Clusters(CS_excited & US_excited) = 3;   % Multisensory
    Clusters(~CS_excited & ~US_excited & (CS_inhibited | US_inhibited)) = 5;  % Inhibited
    Clusters(~CS_excited & ~US_excited & ~CS_inhibited & ~US_inhibited) = 4;  % Non-responsive

    % Compute latencies for sorting and visualization (using no-light conditions)
    CS_onset_lat_nolight = nan(n_light_inhibited, 1);
    CS_offset_lat_nolight = nan(n_light_inhibited, 1);
    US_onset_lat_nolight = nan(n_light_inhibited, 1);
    US_offset_lat_nolight = nan(n_light_inhibited, 1);

    % Also compute for light conditions
    CS_onset_lat_light = nan(n_light_inhibited, 1);
    CS_offset_lat_light = nan(n_light_inhibited, 1);
    US_onset_lat_light = nan(n_light_inhibited, 1);
    US_offset_lat_light = nan(n_light_inhibited, 1);

    for n = 1:n_light_inhibited
        [CS_onset_lat_nolight(n), CS_offset_lat_nolight(n)] = compute_onset_offset_latency(psth_CS_nolight(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [US_onset_lat_nolight(n), US_offset_lat_nolight(n)] = compute_onset_offset_latency(psth_US_nolight(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [CS_onset_lat_light(n), CS_offset_lat_light(n)] = compute_onset_offset_latency(psth_CS_light(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [US_onset_lat_light(n), US_offset_lat_light(n)] = compute_onset_offset_latency(psth_US_light(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
    end

    % Sort neurons (same as figure 2/3)
    leafOrder = [];
    cluster_order = [1, 2, 3, 4, 5];

    for c = cluster_order
        clust_idx = find(Clusters == c);
        if isempty(clust_idx)
            continue;
        end

        if c == 1  % CS-selective
            onset_c = CS_onset_lat_nolight(clust_idx);
            offset_c = CS_offset_lat_nolight(clust_idx);
        elseif c == 2  % US-selective
            onset_c = US_onset_lat_nolight(clust_idx);
            offset_c = US_offset_lat_nolight(clust_idx);
        elseif c == 3  % Multisensory
            CS_duration = CS_offset_lat_nolight(clust_idx) - CS_onset_lat_nolight(clust_idx);
            US_duration = US_offset_lat_nolight(clust_idx) - US_onset_lat_nolight(clust_idx);
            CS_rank_score = CS_onset_lat_nolight(clust_idx) + g.alpha * CS_duration;
            US_rank_score = US_onset_lat_nolight(clust_idx) + g.alpha * US_duration;
            rank_score = mean([CS_rank_score, US_rank_score], 2, 'omitnan');
            sort_matrix = [isnan(rank_score), rank_score];
            [~, sort_idx] = sortrows(sort_matrix, [1 2]);
            leafOrder = [leafOrder; clust_idx(sort_idx)];
            continue;
        elseif c == 4  % Non-responsive
            mean_zscore = mean([mean(psth_CS_nolight(clust_idx, g.roi), 2), mean(psth_US_nolight(clust_idx, g.roi), 2)], 2);
            [~, sort_idx] = sort(mean_zscore, 'descend');
            leafOrder = [leafOrder; clust_idx(sort_idx)];
            continue;
        else  % Inhibited
            mean_zscore = mean([mean(psth_CS_nolight(clust_idx, g.roi), 2), mean(psth_US_nolight(clust_idx, g.roi), 2)], 2);
            [~, sort_idx] = sort(mean_zscore, 'descend');
            leafOrder = [leafOrder; clust_idx(sort_idx)];
            continue;
        end

        duration_c = offset_c - onset_c;
        rank_score = onset_c + g.alpha * duration_c;
        sort_matrix = [isnan(rank_score), rank_score];
        [~, sort_idx] = sortrows(sort_matrix, [1 2]);
        leafOrder = [leafOrder; clust_idx(sort_idx)];
    end

    % Store results
    results_all{br}.Clusters = Clusters;
    results_all{br}.leafOrder = leafOrder;

    % No-light PSTHs
    results_all{br}.psth_CS_nolight = psth_CS_nolight;
    results_all{br}.psth_US_nolight = psth_US_nolight;
    results_all{br}.psth_CS_Hz_nolight = psth_CS_Hz_nolight;
    results_all{br}.psth_US_Hz_nolight = psth_US_Hz_nolight;

    % Light PSTHs
    results_all{br}.psth_CS_light = psth_CS_light;
    results_all{br}.psth_US_light = psth_US_light;
    results_all{br}.psth_CS_Hz_light = psth_CS_Hz_light;
    results_all{br}.psth_US_Hz_light = psth_US_Hz_light;

    % Latencies (no-light)
    results_all{br}.CS_onset_lat_nolight = CS_onset_lat_nolight;
    results_all{br}.CS_offset_lat_nolight = CS_offset_lat_nolight;
    results_all{br}.US_onset_lat_nolight = US_onset_lat_nolight;
    results_all{br}.US_offset_lat_nolight = US_offset_lat_nolight;

    % Latencies (light)
    results_all{br}.CS_onset_lat_light = CS_onset_lat_light;
    results_all{br}.CS_offset_lat_light = CS_offset_lat_light;
    results_all{br}.US_onset_lat_light = US_onset_lat_light;
    results_all{br}.US_offset_lat_light = US_offset_lat_light;

    results_all{br}.n_neurons = n_light_inhibited;
    results_all{br}.light_inhibited_idx = light_inhibited_global;
    results_all{br}.light_inhib_fr_drops = light_inhib_fr_drops_inhibited;
    results_all{br}.light_inhib_pvalues = light_inhib_pvalues_inhibited;
end

%% Check if any light-inhibited neurons found
has_data = false;
for br = 1:2
    if ~isempty(results_all{br})
        has_data = true;
        break;
    end
end

if ~has_data
    return;
end

%% Create figure
% NEW LAYOUT: Row 1=Example rasters (4 cols), Row 2=LA heatmaps (4 cols), Row 3=LA lineplots (4 cols), Row 4=Astria heatmaps (4 cols), Row 5=Astria lineplots (4 cols)
fig = figure('Units', 'pixels', 'Position', [100, 100, 1000, 1000], 'Visible', 'on');
t = tiledlayout(fig, 5, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Determine global color limits across all heatmaps
all_values = [];
for br = 1:2
    if ~isempty(results_all{br})
        for stim = 1:2  % CS, US
            for light_cond = 1:2  % no-light, light
                if stim == 1 && light_cond == 1
                    psth_sorted = results_all{br}.psth_CS_nolight(results_all{br}.leafOrder, :);
                elseif stim == 1 && light_cond == 2
                    psth_sorted = results_all{br}.psth_CS_light(results_all{br}.leafOrder, :);
                elseif stim == 2 && light_cond == 1
                    psth_sorted = results_all{br}.psth_US_nolight(results_all{br}.leafOrder, :);
                else
                    psth_sorted = results_all{br}.psth_US_light(results_all{br}.leafOrder, :);
                end
                matrix = psth_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
                all_values = [all_values; matrix(:)];
            end
        end
    end
end

clim_min = prctile(all_values, (100 - g.clim_percentile) / 2);
clim_max = prctile(all_values, 100 - (100 - g.clim_percentile) / 2);

% Determine global y-limits for lineplots (Hz)
all_lineplot_values = [];
for br = 1:2
    if ~isempty(results_all{br})
        res = results_all{br};
        CS_nolight_hz = res.psth_CS_Hz_nolight(res.leafOrder, :);
        CS_light_hz = res.psth_CS_Hz_light(res.leafOrder, :);
        US_nolight_hz = res.psth_US_Hz_nolight(res.leafOrder, :);
        US_light_hz = res.psth_US_Hz_light(res.leafOrder, :);

        plot_idx = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        all_lineplot_values = [all_lineplot_values; mean(CS_nolight_hz(:, plot_idx), 1)'; ...
                                                     mean(CS_light_hz(:, plot_idx), 1)'; ...
                                                     mean(US_nolight_hz(:, plot_idx), 1)'; ...
                                                     mean(US_light_hz(:, plot_idx), 1)'];
    end
end

ylim_lineplot = [0, max(all_lineplot_values) * 1.1];  % 0 to max + 10% padding

%% Plot Row 1: Example rasters from 7th LA neuron
% Get the 7th neuron from LA heatmap (after sorting)
if ~isempty(results_all{1})  % LA exists
    res_LA = results_all{1};
    if length(res_LA.leafOrder) >= 7
        % leafOrder(7) gives the 7th position in sorted list
        % Map it back to global index using light_inhibited_idx
        sorted_position = res_LA.leafOrder(7);
        example_neuron_idx = res_LA.light_inhibited_idx(sorted_position);

        % Create nested tiledlayout spanning all 4 columns of row 1
        t_nested = tiledlayout(t, 1, 4, 'TileSpacing', 'compact', 'Padding', 'tight');
        t_nested.Layout.Tile = 1;  % Start at tile 1
        t_nested.Layout.TileSpan = [1 4];  % Span 1 row, 4 columns

        % Add title to nested layout
        title(t_nested, 'Example light inhibited neuron', 'FontSize', 12, 'FontWeight', 'bold');

        % Get trial data for all 4 conditions
        ttl_list = {ttl_nolight{1}, ttl_light{1}, ttl_nolight{2}, ttl_light{2}};  % CS no-light, CS light, US no-light, US light
        raster_titles = {'CS - No light', 'CS - Light', 'US - No light', 'US - Light'};

        % Plot raster for each condition
        for col = 1:4
            ax = nexttile(t_nested, col);  % Use nested layout instead of main layout

            % Get spike times for this neuron and TTL
            ttl_current = ttl_list{col};
            [psth_temp, preAP_norm, postAP_norm] = BAfc_psth_spx('cell_metrics', cell_metrics, ...
                'ttl', ttl_current, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

            % Plot raster
            hold on;

            % Add red shaded area for light conditions (cols 2 and 4)
            if col == 2 || col == 4
                % Create patch for shaded area
                patch([-0.5, 0.5, 0.5, -0.5], [0, 0, 1000, 1000], [1 0.8 0.8], ...
                    'EdgeColor', 'none', 'FaceAlpha', 0.6);
            end

            % Get pre and post spike times for this neuron
            preAP = preAP_norm{example_neuron_idx};
            postAP = postAP_norm{example_neuron_idx};
            n_trials = length(preAP);

            for trial = 1:n_trials
                % Pre-stimulus spikes (already in seconds, negative times)
                if ~isempty(preAP{trial})
                    spk_times_pre = preAP{trial};
                    % Filter to plot window (use g.plotwin for consistency with heatmaps)
                    spk_times_pre = spk_times_pre(spk_times_pre >= -g.plotwin(1) & spk_times_pre <= 0);
                    if ~isempty(spk_times_pre)
                        plot(spk_times_pre, trial * ones(size(spk_times_pre)), 'k.', 'MarkerSize', 6);
                    end
                end

                % Post-stimulus spikes (already in seconds, positive times)
                if ~isempty(postAP{trial})
                    spk_times_post = postAP{trial};
                    % Filter to plot window (use g.plotwin for consistency with heatmaps)
                    spk_times_post = spk_times_post(spk_times_post >= 0 & spk_times_post <= g.plotwin(2));
                    if ~isempty(spk_times_post)
                        plot(spk_times_post, trial * ones(size(spk_times_post)), 'k.', 'MarkerSize', 6);
                    end
                end
            end

            % Vertical line at stimulus onset
            plot([0 0], [0 n_trials+1], 'r-', 'LineWidth', 2);

            % Red line for light illumination on light conditions (cols 2 and 4)
            if col == 2 || col == 4
                y_pos = n_trials + 2;  % Above all trials
                plot([-0.5, 0.5], [y_pos, y_pos], '-', 'Color', [1 0 0], 'LineWidth', 3);
            end

            hold off;

            % Formatting
            xlim([-g.plotwin(1), g.plotwin(2)]);
            ylim([0, n_trials + 1]);

            title(raster_titles{col}, 'FontSize', g.fontSize1);
            if col == 1
                ylabel('Trial', 'FontSize', 12);
                yticks([0 40 80]);
            else
                set(gca, 'YTickLabel', []);
            end
            set(gca, 'XTickLabel', []);  % No x-labels on row 1
            set(gca, 'FontSize', g.fontSize2);
        end
    end
end

%% Plot: Row 2=LA heatmaps, Row 3=LA lineplots, Row 4=Astria heatmaps, Row 5=Astria lineplots
% 4 columns per row: CS no-light, CS light, US no-light, US light
for br = 1:2  % LA, Astria
    if isempty(results_all{br})
        continue;
    end

    res = results_all{br};

    % Get cluster boundaries
    Clusters_sorted = res.Clusters(res.leafOrder);
    n_clu = find(diff(Clusters_sorted) ~= 0);

    % Get all PSTHs sorted
    CS_nolight_z = res.psth_CS_nolight(res.leafOrder, :);
    CS_light_z = res.psth_CS_light(res.leafOrder, :);
    US_nolight_z = res.psth_US_nolight(res.leafOrder, :);
    US_light_z = res.psth_US_light(res.leafOrder, :);
    CS_nolight_hz = res.psth_CS_Hz_nolight(res.leafOrder, :);
    CS_light_hz = res.psth_CS_Hz_light(res.leafOrder, :);
    US_nolight_hz = res.psth_US_Hz_nolight(res.leafOrder, :);
    US_light_hz = res.psth_US_Hz_light(res.leafOrder, :);

    % Calculate row indices: LA gets rows 2-3, Astria gets rows 4-5 (row 1 is for example rasters)
    heatmap_row = (br-1)*2 + 2;
    lineplot_row = (br-1)*2 + 3;

    % Create nested tiledlayout spanning 2 rows and 4 columns for this brain region
    t_nested_region = tiledlayout(t, 2, 4, 'TileSpacing', 'compact', 'Padding', 'tight');
    first_tile = (heatmap_row - 1) * 4 + 1;  % First tile of heatmap row
    t_nested_region.Layout.Tile = first_tile;
    t_nested_region.Layout.TileSpan = [2 4];  % Span 2 rows, 4 columns

    % Add title to nested layout
    if br == 1
        title(t_nested_region, 'All light-inhibited neurons from LA', 'FontSize', 12, 'FontWeight', 'bold');
    else
        title(t_nested_region, 'All light-inhibited neurons from AStria', 'FontSize', 12, 'FontWeight', 'bold');
    end

    % Prepare all PSTH data in order: CS no-light, CS light, US no-light, US light
    psth_z_all = {CS_nolight_z, CS_light_z, US_nolight_z, US_light_z};
    psth_hz_all = {CS_nolight_hz, CS_light_hz, US_nolight_hz, US_light_hz};
    onset_lat_all = {res.CS_onset_lat_nolight, res.CS_onset_lat_light, res.US_onset_lat_nolight, res.US_onset_lat_light};
    offset_lat_all = {res.CS_offset_lat_nolight, res.CS_offset_lat_light, res.US_offset_lat_nolight, res.US_offset_lat_light};
    titles_all = {'CS - No light', 'CS - Light', 'US - No light', 'US - Light'};
    line_colors = {[0.2 0.6 0.2], [0.8 0.2 0.2], [0.2 0.6 0.2], [0.8 0.2 0.2]};

    plot_idx = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    time_vec = g.timeaxis_hmp;
    if length(time_vec) ~= length(plot_idx)
        time_vec = linspace(-g.plotwin(1), g.plotwin(2), length(plot_idx));
    end

    % Plot all 4 conditions (columns)
    for col = 1:4
        %% Row 1 of nested layout: Heatmap
        tile_idx = col;  % Just use column index within nested layout
        ax = nexttile(t_nested_region, tile_idx);
        matrix = psth_z_all{col}(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        imagesc(g.timeaxis_hmp, 1:size(matrix, 1), matrix);
        clim([clim_min clim_max]);
        colormap(ax, g.colors.Heatmap);
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        % Labels (no titles on heatmap rows - titles are on row 1 rasters)
        if col == 1  % First column gets ylabel
            if br == 1
                ylabel('LA neurons', 'FontSize', 12);
            else
                ylabel('AStria neurons', 'FontSize', 12);
            end
        end

        yticks([1, size(matrix, 1)]);
        if col > 1
            set(gca, 'YTickLabel', []);
        end
        set(gca, 'XTickLabel', []);  % Heatmap rows never show x-labels
        set(gca, 'FontSize', g.fontSize2);

        % Add cluster lines
        hold on;
        for i = 1:length(n_clu)
            yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
        end

        % Add onset/offset markers
        for n = 1:length(res.leafOrder)
            idx_n = res.leafOrder(n);
            if ~isnan(onset_lat_all{col}(idx_n))
                plot(onset_lat_all{col}(idx_n), n, 'k.', 'MarkerSize', 4);
            end
            if ~isnan(offset_lat_all{col}(idx_n))
                plot(offset_lat_all{col}(idx_n), n, 'k.', 'MarkerSize', 4);
            end
        end

        % Add red line for light illumination period on light conditions (cols 2 and 4)
        if col == 2 || col == 4
            y_pos = 0.5;  % At the top of heatmap (y-axis is flipped, so lowest value is at top)
            plot([-0.5, 0.5], [y_pos, y_pos], '-', 'Color', [1 0 0], 'LineWidth', 3);
        end

        hold off;

        %% Row 2 of nested layout: Lineplot
        tile_idx = 4 + col;  % Second row starts at tile 5
        ax = nexttile(t_nested_region, tile_idx);

        psth_mean = mean(psth_hz_all{col}(:, plot_idx), 1);

        % Set limits first
        xlim([time_vec(1) time_vec(end)]);
        ylim([0, 30]);

        hold on;

        % Add red shaded area for light conditions (cols 2 and 4)
        if col == 2 || col == 4
            patch([-0.5, 0.5, 0.5, -0.5], [0, 0, 30, 30], [1 0.8 0.8], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end

        xline(0, '--k', 'LineWidth', g.xlinewidth);
        plot(time_vec, psth_mean, '-', 'Color', line_colors{col}, 'LineWidth', 2.5);

        hold off;

        % Y-axis labels and ticks
        if col == 1  % First column shows ylabel and ticks
            ylabel('FR (Hz)', 'FontSize', 12);
            yticks([0 15 30]);
        else
            set(gca, 'YTickLabel', []);
        end

        % X-axis labels
        if br == 2  % Astria lineplot (row 5) - bottom row gets x-labels
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
        end
        set(gca, 'FontSize', g.fontSize2);
    end
end

% Add colorbar
drawnow;
cb_width = 0.008;
cb_left = 0.95;
cb_bottom = 0.605;
cb_height = 0.08;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'right');
%ylabel(cb_ax, 'Z-score', 'FontSize', 12);
cb_ax.FontSize = g.fontSize2;

%% Add panel labels
% A: Raster plots (row 1)
annotation(fig, 'textbox', [0.01, 0.95, 0.03, 0.02], ...
    'String', 'A', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% B: LA heatmaps (row 2)
annotation(fig, 'textbox', [0.01, 0.77, 0.03, 0.02], ...
    'String', 'B', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% C: LA lineplots (row 3)
annotation(fig, 'textbox', [0.01, 0.59, 0.03, 0.02], ...
    'String', 'C', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% D: AStria heatmaps (row 4)
annotation(fig, 'textbox', [0.01, 0.41, 0.03, 0.02], ...
    'String', 'D', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

% E: AStria lineplots (row 5)
annotation(fig, 'textbox', [0.01, 0.23, 0.03, 0.02], ...
    'String', 'E', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');

%% Export data to Excel
export_supp_fig3_to_excel(results_all, brain_regions, g, cell_metrics);
exportgraphics(gcf, 'supplementary_figure_3.png', 'Resolution', 300);
%% Helper function
function [onset_lat, offset_lat] = compute_onset_offset_latency(z_trace, event_inds, threshold, min_consec, bin_time)
    seg = z_trace(event_inds);
    if any(isnan(seg))
        onset_lat = NaN;
        offset_lat = NaN;
        return
    end

    isAbove = seg >= threshold;
    winsum = conv(double(isAbove), ones(1, min_consec), 'same');
    onset_idx = find(winsum >= min_consec, 1, 'first');

    if isempty(onset_idx)
        onset_lat = NaN;
        offset_lat = NaN;
    else
        onset_lat = (onset_idx - 1) * bin_time;

        seg_after_onset = seg(onset_idx:end);
        isBelow = seg_after_onset < threshold;
        winsum_below = conv(double(isBelow), ones(1, min_consec), 'same');
        offset_idx_relative = find(winsum_below >= min_consec, 1, 'first');

        if isempty(offset_idx_relative)
            offset_lat = NaN;
        else
            offset_idx = onset_idx + offset_idx_relative - 1;
            offset_lat = (offset_idx - 1) * bin_time;
        end
    end
end

function export_supp_fig3_to_excel(results_all, brain_regions, g, cell_metrics)
    output_filename = 'supplementary_figure_3_data.xlsx';

    if exist(output_filename, 'file')
        delete(output_filename);
    end

    cluster_names = {'CS-selective', 'US-selective', 'Multisensory', 'Non-responsive', 'Inhibited'};

    %% Sheet 1: Overall Sample Sizes
    sheet_data = {};
    row = 1;

    sheet_data{row, 1} = 'SUPPLEMENTARY FIGURE 3 DATA';
    row = row + 1;
    sheet_data{row, 1} = 'Light-Inhibited Neurons Analysis';
    row = row + 1;
    sheet_data{row, 1} = sprintf('Generated: %s', datestr(now));
    row = row + 2;

    sheet_data{row, 1} = 'OVERALL SAMPLE SIZES';
    row = row + 2;

    % Collect all animals
    all_animals = {};
    for br = 1:2
        if ~isempty(results_all{br})
            neuron_indices = results_all{br}.light_inhibited_idx;
            for n = 1:length(neuron_indices)
                ii = neuron_indices(n);
                all_animals = [all_animals; {cell_metrics.animal{ii}}];
            end
        end
    end
    unique_animals = unique(all_animals);

    sheet_data{row, 1} = 'Number of animals (N):';
    sheet_data{row, 2} = length(unique_animals);
    row = row + 1;
    sheet_data{row, 1} = 'Animal IDs:';
    sheet_data{row, 2} = strjoin(unique_animals, ', ');
    row = row + 2;

    total_neurons = 0;
    for br = 1:2
        if ~isempty(results_all{br})
            n_neurons = results_all{br}.n_neurons;
            total_neurons = total_neurons + n_neurons;
            if strcmp(brain_regions{br}, 'Astria')
                region_name = 'AStria';
            else
                region_name = brain_regions{br};
            end
            sheet_data{row, 1} = sprintf('%s:', region_name);
            sheet_data{row, 2} = sprintf('n = %d light-inhibited neurons', n_neurons);
            row = row + 1;
        end
    end
    sheet_data{row, 1} = 'Total (LA + AStria):';
    sheet_data{row, 2} = total_neurons;
    row = row + 2;

    % Detection criteria
    sheet_data{row, 1} = 'DETECTION CRITERIA';
    row = row + 1;
    sheet_data{row, 1} = 'Baseline window:';
    sheet_data{row, 2} = sprintf('-%.1fs to -%.1fs before stimulus', ...
        g.light_inhib_window_baseline + g.light_inhib_window_recent, g.light_inhib_window_recent);
    row = row + 1;
    sheet_data{row, 1} = 'Recent window:';
    sheet_data{row, 2} = sprintf('-%.1fs to 0s (light illumination)', g.light_inhib_window_recent);
    row = row + 1;
    sheet_data{row, 1} = 'Statistical test:';
    sheet_data{row, 2} = 'Wilcoxon ranksum (right-tailed)';
    row = row + 1;
    sheet_data{row, 1} = 'Significance threshold:';
    sheet_data{row, 2} = sprintf('p < %.2f', g.light_inhib_p_threshold);
    row = row + 1;
    sheet_data{row, 1} = 'FR drop threshold:';
    sheet_data{row, 2} = sprintf('%.0f%%', g.light_inhib_fr_drop * 100);
    row = row + 2;

    % Cluster distributions
    sheet_data{row, 1} = 'CLUSTER DISTRIBUTIONS';
    row = row + 2;

    for br = 1:2
        if ~isempty(results_all{br})
            res = results_all{br};
            if strcmp(brain_regions{br}, 'Astria')
                region_name = 'AStria';
            else
                region_name = brain_regions{br};
            end

            sheet_data{row, 1} = region_name;
            row = row + 1;
            for c = 1:5
                n_cluster = sum(res.Clusters == c);
                pct = 100 * n_cluster / res.n_neurons;
                sheet_data{row, 1} = cluster_names{c};
                sheet_data{row, 2} = n_cluster;
                sheet_data{row, 3} = sprintf('(%.1f%%)', pct);
                row = row + 1;
            end
            row = row + 1;
        end
    end

    writecell(sheet_data, output_filename, 'Sheet', 'Summary');

    %% Sheet 2: PV Silencing Effect (pooled CS+US, -0.5-0s window)
    baseline_idx = round((g.pre_time - 0.5)/g.bin_time + 1 : g.pre_time/g.bin_time);

    sheet_data = {};
    row = 1;
    sheet_data{row, 1} = 'PV SILENCING EFFECT (Pre-stimulus period: -0.5s to 0s)';
    row = row + 1;
    sheet_data{row, 1} = 'Pooled across CS and US trials';
    row = row + 3;

    for br = 1:2
        if ~isempty(results_all{br})
            res = results_all{br};
            if strcmp(brain_regions{br}, 'Astria')
                region_name = 'AStria';
            else
                region_name = brain_regions{br};
            end

            % Pool CS + US trials
            psth_nolight_pooled = [res.psth_CS_Hz_nolight; res.psth_US_Hz_nolight];
            psth_light_pooled = [res.psth_CS_Hz_light; res.psth_US_Hz_light];

            fr_nolight = mean(psth_nolight_pooled(:, baseline_idx), 2);
            fr_light = mean(psth_light_pooled(:, baseline_idx), 2);

            sheet_data{row, 1} = region_name;
            sheet_data{row, 2} = sprintf('(n = %d)', res.n_neurons);
            row = row + 1;

            % Add column headers
            sheet_data{row, 1} = '';
            sheet_data{row, 2} = 'Mean (Hz)';
            sheet_data{row, 3} = 'SEM (Hz)';
            sheet_data{row, 4} = 'Median (Hz)';
            row = row + 1;

            sheet_data{row, 1} = 'No-light FR:';
            sheet_data{row, 2} = mean(fr_nolight);
            sheet_data{row, 3} = std(fr_nolight)/sqrt(length(fr_nolight));
            sheet_data{row, 4} = median(fr_nolight);
            row = row + 1;

            sheet_data{row, 1} = 'Light FR:';
            sheet_data{row, 2} = mean(fr_light);
            sheet_data{row, 3} = std(fr_light)/sqrt(length(fr_light));
            sheet_data{row, 4} = median(fr_light);
            row = row + 2;

            [p_val, ~] = signrank(fr_nolight, fr_light);
            sheet_data{row, 1} = 'Statistical Test:';
            sheet_data{row, 2} = 'Wilcoxon signed-rank (paired)';
            row = row + 1;
            sheet_data{row, 1} = 'P-value:';
            sheet_data{row, 2} = p_val;
            sheet_data{row, 3} = format_significance(p_val);
            row = row + 2;

            fr_change = mean(fr_light) - mean(fr_nolight);
            fr_change_pct = (fr_change / mean(fr_nolight)) * 100;
            sheet_data{row, 1} = 'Mean change:';
            sheet_data{row, 2} = sprintf('%.2f Hz (%.1f%%)', fr_change, fr_change_pct);
            row = row + 3;
        end
    end

    writecell(sheet_data, output_filename, 'Sheet', 'PV_Silencing_Effect');

    %% Sheet 3: Firing Rate Analysis (3 time windows, all conditions)
    % Define time windows
    baseline_early_idx = round((g.pre_time - 1)/g.bin_time + 1 : (g.pre_time - 0.5)/g.bin_time);  % -1 to -0.5s
    baseline_late_idx = round((g.pre_time - 0.5)/g.bin_time + 1 : g.pre_time/g.bin_time);  % -0.5 to 0s (illumination)
    response_idx = round(g.pre_time/g.bin_time + 1 : (g.pre_time + 0.5)/g.bin_time);  % 0 to 0.5s (response)

    % Pool data from both regions (LA + AStria)
    all_data_CS_nolight = [];
    all_data_CS_light = [];
    all_data_US_nolight = [];
    all_data_US_light = [];

    for br = 1:2
        if ~isempty(results_all{br})
            res = results_all{br};
            all_data_CS_nolight = [all_data_CS_nolight; res.psth_CS_Hz_nolight];
            all_data_CS_light = [all_data_CS_light; res.psth_CS_Hz_light];
            all_data_US_nolight = [all_data_US_nolight; res.psth_US_Hz_nolight];
            all_data_US_light = [all_data_US_light; res.psth_US_Hz_light];
        end
    end

    n_total = size(all_data_CS_nolight, 1);

    % Calculate FR for each neuron in each window
    % CS no-light
    fr_CS_nolight_baseline = mean(all_data_CS_nolight(:, baseline_early_idx), 2);
    fr_CS_nolight_illum = mean(all_data_CS_nolight(:, baseline_late_idx), 2);
    fr_CS_nolight_resp = mean(all_data_CS_nolight(:, response_idx), 2);

    % CS light
    fr_CS_light_baseline = mean(all_data_CS_light(:, baseline_early_idx), 2);
    fr_CS_light_illum = mean(all_data_CS_light(:, baseline_late_idx), 2);
    fr_CS_light_resp = mean(all_data_CS_light(:, response_idx), 2);

    % US no-light
    fr_US_nolight_baseline = mean(all_data_US_nolight(:, baseline_early_idx), 2);
    fr_US_nolight_illum = mean(all_data_US_nolight(:, baseline_late_idx), 2);
    fr_US_nolight_resp = mean(all_data_US_nolight(:, response_idx), 2);

    % US light
    fr_US_light_baseline = mean(all_data_US_light(:, baseline_early_idx), 2);
    fr_US_light_illum = mean(all_data_US_light(:, baseline_late_idx), 2);
    fr_US_light_resp = mean(all_data_US_light(:, response_idx), 2);

    sheet_data = {};
    row = 1;
    sheet_data{row, 1} = 'FIRING RATE ANALYSIS (3 Time Windows)';
    row = row + 1;
    sheet_data{row, 1} = sprintf('Pooled LA + AStria (n = %d neurons)', n_total);
    row = row + 1;
    sheet_data{row, 1} = 'Time windows: -1 to -0.5s (baseline), -0.5 to 0s (illumination), 0 to 0.5s (response)';
    row = row + 3;

    % CS No-light
    sheet_data{row, 1} = 'CS - No light';
    row = row + 1;
    sheet_data{row, 1} = 'Time Window';
    sheet_data{row, 2} = 'Mean (Hz)';
    sheet_data{row, 3} = 'SEM (Hz)';
    sheet_data{row, 4} = 'Median (Hz)';
    sheet_data{row, 5} = 'SD (Hz)';
    row = row + 1;
    sheet_data{row, 1} = '-1 to -0.5s (baseline)';
    sheet_data{row, 2} = mean(fr_CS_nolight_baseline);
    sheet_data{row, 3} = std(fr_CS_nolight_baseline)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_CS_nolight_baseline);
    sheet_data{row, 5} = std(fr_CS_nolight_baseline);
    row = row + 1;
    sheet_data{row, 1} = '-0.5 to 0s (illumination)';
    sheet_data{row, 2} = mean(fr_CS_nolight_illum);
    sheet_data{row, 3} = std(fr_CS_nolight_illum)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_CS_nolight_illum);
    sheet_data{row, 5} = std(fr_CS_nolight_illum);
    row = row + 1;
    sheet_data{row, 1} = '0 to 0.5s (response)';
    sheet_data{row, 2} = mean(fr_CS_nolight_resp);
    sheet_data{row, 3} = std(fr_CS_nolight_resp)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_CS_nolight_resp);
    sheet_data{row, 5} = std(fr_CS_nolight_resp);
    row = row + 2;

    % CS Light
    sheet_data{row, 1} = 'CS - Light';
    row = row + 1;
    sheet_data{row, 1} = 'Time Window';
    sheet_data{row, 2} = 'Mean (Hz)';
    sheet_data{row, 3} = 'SEM (Hz)';
    sheet_data{row, 4} = 'Median (Hz)';
    sheet_data{row, 5} = 'SD (Hz)';
    row = row + 1;
    sheet_data{row, 1} = '-1 to -0.5s (baseline)';
    sheet_data{row, 2} = mean(fr_CS_light_baseline);
    sheet_data{row, 3} = std(fr_CS_light_baseline)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_CS_light_baseline);
    sheet_data{row, 5} = std(fr_CS_light_baseline);
    row = row + 1;
    sheet_data{row, 1} = '-0.5 to 0s (illumination)';
    sheet_data{row, 2} = mean(fr_CS_light_illum);
    sheet_data{row, 3} = std(fr_CS_light_illum)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_CS_light_illum);
    sheet_data{row, 5} = std(fr_CS_light_illum);
    row = row + 1;
    sheet_data{row, 1} = '0 to 0.5s (response)';
    sheet_data{row, 2} = mean(fr_CS_light_resp);
    sheet_data{row, 3} = std(fr_CS_light_resp)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_CS_light_resp);
    sheet_data{row, 5} = std(fr_CS_light_resp);
    row = row + 2;

    % US No-light
    sheet_data{row, 1} = 'US - No light';
    row = row + 1;
    sheet_data{row, 1} = 'Time Window';
    sheet_data{row, 2} = 'Mean (Hz)';
    sheet_data{row, 3} = 'SEM (Hz)';
    sheet_data{row, 4} = 'Median (Hz)';
    sheet_data{row, 5} = 'SD (Hz)';
    row = row + 1;
    sheet_data{row, 1} = '-1 to -0.5s (baseline)';
    sheet_data{row, 2} = mean(fr_US_nolight_baseline);
    sheet_data{row, 3} = std(fr_US_nolight_baseline)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_US_nolight_baseline);
    sheet_data{row, 5} = std(fr_US_nolight_baseline);
    row = row + 1;
    sheet_data{row, 1} = '-0.5 to 0s (illumination)';
    sheet_data{row, 2} = mean(fr_US_nolight_illum);
    sheet_data{row, 3} = std(fr_US_nolight_illum)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_US_nolight_illum);
    sheet_data{row, 5} = std(fr_US_nolight_illum);
    row = row + 1;
    sheet_data{row, 1} = '0 to 0.5s (response)';
    sheet_data{row, 2} = mean(fr_US_nolight_resp);
    sheet_data{row, 3} = std(fr_US_nolight_resp)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_US_nolight_resp);
    sheet_data{row, 5} = std(fr_US_nolight_resp);
    row = row + 2;

    % US Light
    sheet_data{row, 1} = 'US - Light';
    row = row + 1;
    sheet_data{row, 1} = 'Time Window';
    sheet_data{row, 2} = 'Mean (Hz)';
    sheet_data{row, 3} = 'SEM (Hz)';
    sheet_data{row, 4} = 'Median (Hz)';
    sheet_data{row, 5} = 'SD (Hz)';
    row = row + 1;
    sheet_data{row, 1} = '-1 to -0.5s (baseline)';
    sheet_data{row, 2} = mean(fr_US_light_baseline);
    sheet_data{row, 3} = std(fr_US_light_baseline)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_US_light_baseline);
    sheet_data{row, 5} = std(fr_US_light_baseline);
    row = row + 1;
    sheet_data{row, 1} = '-0.5 to 0s (illumination)';
    sheet_data{row, 2} = mean(fr_US_light_illum);
    sheet_data{row, 3} = std(fr_US_light_illum)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_US_light_illum);
    sheet_data{row, 5} = std(fr_US_light_illum);
    row = row + 1;
    sheet_data{row, 1} = '0 to 0.5s (response)';
    sheet_data{row, 2} = mean(fr_US_light_resp);
    sheet_data{row, 3} = std(fr_US_light_resp)/sqrt(n_total);
    sheet_data{row, 4} = median(fr_US_light_resp);
    sheet_data{row, 5} = std(fr_US_light_resp);
    row = row + 3;

    % Statistical comparisons
    sheet_data{row, 1} = 'STATISTICAL COMPARISONS';
    row = row + 2;

    % 1. Baseline vs Illumination (baseline early vs baseline late) - pooled CS+US no-light
    sheet_data{row, 1} = '1. Illumination Effect (No-light trials)';
    row = row + 1;
    sheet_data{row, 1} = 'Comparison: Baseline (-1 to -0.5s) vs Illumination (-0.5 to 0s)';
    row = row + 1;
    sheet_data{row, 1} = 'Pooled CS + US, no-light trials';
    row = row + 2;

    fr_pooled_nolight_baseline = [fr_CS_nolight_baseline; fr_US_nolight_baseline];
    fr_pooled_nolight_illum = [fr_CS_nolight_illum; fr_US_nolight_illum];
    [p_baseline_vs_illum, ~] = signrank(fr_pooled_nolight_baseline, fr_pooled_nolight_illum);

    sheet_data{row, 1} = 'Baseline FR:';
    sheet_data{row, 2} = mean(fr_pooled_nolight_baseline);
    sheet_data{row, 3} = std(fr_pooled_nolight_baseline)/sqrt(length(fr_pooled_nolight_baseline));
    row = row + 1;
    sheet_data{row, 1} = 'Illumination FR:';
    sheet_data{row, 2} = mean(fr_pooled_nolight_illum);
    sheet_data{row, 3} = std(fr_pooled_nolight_illum)/sqrt(length(fr_pooled_nolight_illum));
    row = row + 1;
    sheet_data{row, 1} = 'Wilcoxon signed-rank test p-value:';
    sheet_data{row, 2} = p_baseline_vs_illum;
    sheet_data{row, 3} = format_significance(p_baseline_vs_illum);
    row = row + 3;

    % 2. Response comparison: No-light vs Light (0 to 0.5s) - CS trials
    sheet_data{row, 1} = '2. PV Silencing Effect on CS Response (0 to 0.5s)';
    row = row + 1;
    sheet_data{row, 1} = 'Comparison: CS No-light vs CS Light response';
    row = row + 2;

    [p_CS_resp, ~] = signrank(fr_CS_nolight_resp, fr_CS_light_resp);

    sheet_data{row, 1} = 'CS No-light response FR:';
    sheet_data{row, 2} = mean(fr_CS_nolight_resp);
    sheet_data{row, 3} = std(fr_CS_nolight_resp)/sqrt(n_total);
    row = row + 1;
    sheet_data{row, 1} = 'CS Light response FR:';
    sheet_data{row, 2} = mean(fr_CS_light_resp);
    sheet_data{row, 3} = std(fr_CS_light_resp)/sqrt(n_total);
    row = row + 1;
    sheet_data{row, 1} = 'Wilcoxon signed-rank test p-value:';
    sheet_data{row, 2} = p_CS_resp;
    sheet_data{row, 3} = format_significance(p_CS_resp);
    row = row + 3;

    % 3. Response comparison: No-light vs Light (0 to 0.5s) - US trials
    sheet_data{row, 1} = '3. PV Silencing Effect on US Response (0 to 0.5s)';
    row = row + 1;
    sheet_data{row, 1} = 'Comparison: US No-light vs US Light response';
    row = row + 2;

    [p_US_resp, ~] = signrank(fr_US_nolight_resp, fr_US_light_resp);

    sheet_data{row, 1} = 'US No-light response FR:';
    sheet_data{row, 2} = mean(fr_US_nolight_resp);
    sheet_data{row, 3} = std(fr_US_nolight_resp)/sqrt(n_total);
    row = row + 1;
    sheet_data{row, 1} = 'US Light response FR:';
    sheet_data{row, 2} = mean(fr_US_light_resp);
    sheet_data{row, 3} = std(fr_US_light_resp)/sqrt(n_total);
    row = row + 1;
    sheet_data{row, 1} = 'Wilcoxon signed-rank test p-value:';
    sheet_data{row, 2} = p_US_resp;
    sheet_data{row, 3} = format_significance(p_US_resp);
    row = row + 1;

    writecell(sheet_data, output_filename, 'Sheet', 'Firing_Rate_Analysis');

    %% Sheet 4: Complete Neuron Lists with Firing Rate Values
    sheet_data = {};
    row = 1;
    sheet_data{row, 1} = 'COMPLETE NEURON LISTS (sorted order as in heatmaps)';
    row = row + 1;
    sheet_data{row, 1} = 'Firing rate values calculated in three time windows for all conditions';
    row = row + 3;

    % Get the pooled firing rate data calculated earlier (need to recalculate per region)
    for br = 1:2
        if ~isempty(results_all{br})
            res = results_all{br};
            if strcmp(brain_regions{br}, 'Astria')
                region_name = 'AStria';
            else
                region_name = brain_regions{br};
            end

            % Calculate FR for each neuron in each window
            baseline_early_idx = round((g.pre_time - 1)/g.bin_time + 1 : (g.pre_time - 0.5)/g.bin_time);
            baseline_late_idx = round((g.pre_time - 0.5)/g.bin_time + 1 : g.pre_time/g.bin_time);
            response_idx = round(g.pre_time/g.bin_time + 1 : (g.pre_time + 0.5)/g.bin_time);

            % CS no-light
            fr_CS_nolight_baseline = mean(res.psth_CS_Hz_nolight(:, baseline_early_idx), 2);
            fr_CS_nolight_illum = mean(res.psth_CS_Hz_nolight(:, baseline_late_idx), 2);
            fr_CS_nolight_resp = mean(res.psth_CS_Hz_nolight(:, response_idx), 2);

            % CS light
            fr_CS_light_baseline = mean(res.psth_CS_Hz_light(:, baseline_early_idx), 2);
            fr_CS_light_illum = mean(res.psth_CS_Hz_light(:, baseline_late_idx), 2);
            fr_CS_light_resp = mean(res.psth_CS_Hz_light(:, response_idx), 2);

            % US no-light
            fr_US_nolight_baseline = mean(res.psth_US_Hz_nolight(:, baseline_early_idx), 2);
            fr_US_nolight_illum = mean(res.psth_US_Hz_nolight(:, baseline_late_idx), 2);
            fr_US_nolight_resp = mean(res.psth_US_Hz_nolight(:, response_idx), 2);

            % US light
            fr_US_light_baseline = mean(res.psth_US_Hz_light(:, baseline_early_idx), 2);
            fr_US_light_illum = mean(res.psth_US_Hz_light(:, baseline_late_idx), 2);
            fr_US_light_resp = mean(res.psth_US_Hz_light(:, response_idx), 2);

            sheet_data{row, 1} = region_name;
            sheet_data{row, 2} = sprintf('(n = %d)', res.n_neurons);
            row = row + 1;

            % Header row
            sheet_data{row, 1} = 'Position';
            sheet_data{row, 2} = 'Animal';
            sheet_data{row, 3} = 'Cell ID';
            sheet_data{row, 4} = 'Region';
            sheet_data{row, 5} = 'Cell Type';
            sheet_data{row, 6} = 'Cluster';
            sheet_data{row, 7} = 'P-value';
            sheet_data{row, 8} = 'FR drop (%)';
            sheet_data{row, 9} = '';  % Separator
            sheet_data{row, 10} = 'CS NoLight -1 to -0.5s (Hz)';
            sheet_data{row, 11} = 'CS NoLight -0.5 to 0s (Hz)';
            sheet_data{row, 12} = 'CS NoLight 0 to 0.5s (Hz)';
            sheet_data{row, 13} = '';  % Separator
            sheet_data{row, 14} = 'CS Light -1 to -0.5s (Hz)';
            sheet_data{row, 15} = 'CS Light -0.5 to 0s (Hz)';
            sheet_data{row, 16} = 'CS Light 0 to 0.5s (Hz)';
            sheet_data{row, 17} = '';  % Separator
            sheet_data{row, 18} = 'US NoLight -1 to -0.5s (Hz)';
            sheet_data{row, 19} = 'US NoLight -0.5 to 0s (Hz)';
            sheet_data{row, 20} = 'US NoLight 0 to 0.5s (Hz)';
            sheet_data{row, 21} = '';  % Separator
            sheet_data{row, 22} = 'US Light -1 to -0.5s (Hz)';
            sheet_data{row, 23} = 'US Light -0.5 to 0s (Hz)';
            sheet_data{row, 24} = 'US Light 0 to 0.5s (Hz)';
            row = row + 1;

            % Data rows
            for pos = 1:length(res.leafOrder)
                sorted_idx = res.leafOrder(pos);
                global_idx = res.light_inhibited_idx(sorted_idx);

                cell_type_str = 'N/A';
                if isfield(cell_metrics, 'putativeCellType')
                    cell_type_str = cell_metrics.putativeCellType{global_idx};
                end

                sheet_data{row, 1} = pos;
                sheet_data{row, 2} = cell_metrics.animal{global_idx};
                sheet_data{row, 3} = cell_metrics.cellID(global_idx);
                sheet_data{row, 4} = cell_metrics.brainRegion{global_idx};
                sheet_data{row, 5} = cell_type_str;
                sheet_data{row, 6} = cluster_names{res.Clusters(sorted_idx)};
                sheet_data{row, 7} = res.light_inhib_pvalues(sorted_idx);
                sheet_data{row, 8} = res.light_inhib_fr_drops(sorted_idx) * 100;
                sheet_data{row, 9} = '';  % Separator

                % CS no-light
                sheet_data{row, 10} = fr_CS_nolight_baseline(sorted_idx);
                sheet_data{row, 11} = fr_CS_nolight_illum(sorted_idx);
                sheet_data{row, 12} = fr_CS_nolight_resp(sorted_idx);
                sheet_data{row, 13} = '';  % Separator

                % CS light
                sheet_data{row, 14} = fr_CS_light_baseline(sorted_idx);
                sheet_data{row, 15} = fr_CS_light_illum(sorted_idx);
                sheet_data{row, 16} = fr_CS_light_resp(sorted_idx);
                sheet_data{row, 17} = '';  % Separator

                % US no-light
                sheet_data{row, 18} = fr_US_nolight_baseline(sorted_idx);
                sheet_data{row, 19} = fr_US_nolight_illum(sorted_idx);
                sheet_data{row, 20} = fr_US_nolight_resp(sorted_idx);
                sheet_data{row, 21} = '';  % Separator

                % US light
                sheet_data{row, 22} = fr_US_light_baseline(sorted_idx);
                sheet_data{row, 23} = fr_US_light_illum(sorted_idx);
                sheet_data{row, 24} = fr_US_light_resp(sorted_idx);

                row = row + 1;
            end

            % Add mean and SEM rows for verification
            row = row + 1;
            sheet_data{row, 1} = 'Mean:';
            sheet_data{row, 10} = mean(fr_CS_nolight_baseline);
            sheet_data{row, 11} = mean(fr_CS_nolight_illum);
            sheet_data{row, 12} = mean(fr_CS_nolight_resp);
            sheet_data{row, 14} = mean(fr_CS_light_baseline);
            sheet_data{row, 15} = mean(fr_CS_light_illum);
            sheet_data{row, 16} = mean(fr_CS_light_resp);
            sheet_data{row, 18} = mean(fr_US_nolight_baseline);
            sheet_data{row, 19} = mean(fr_US_nolight_illum);
            sheet_data{row, 20} = mean(fr_US_nolight_resp);
            sheet_data{row, 22} = mean(fr_US_light_baseline);
            sheet_data{row, 23} = mean(fr_US_light_illum);
            sheet_data{row, 24} = mean(fr_US_light_resp);
            row = row + 1;

            sheet_data{row, 1} = 'SEM:';
            sheet_data{row, 10} = std(fr_CS_nolight_baseline)/sqrt(res.n_neurons);
            sheet_data{row, 11} = std(fr_CS_nolight_illum)/sqrt(res.n_neurons);
            sheet_data{row, 12} = std(fr_CS_nolight_resp)/sqrt(res.n_neurons);
            sheet_data{row, 14} = std(fr_CS_light_baseline)/sqrt(res.n_neurons);
            sheet_data{row, 15} = std(fr_CS_light_illum)/sqrt(res.n_neurons);
            sheet_data{row, 16} = std(fr_CS_light_resp)/sqrt(res.n_neurons);
            sheet_data{row, 18} = std(fr_US_nolight_baseline)/sqrt(res.n_neurons);
            sheet_data{row, 19} = std(fr_US_nolight_illum)/sqrt(res.n_neurons);
            sheet_data{row, 20} = std(fr_US_nolight_resp)/sqrt(res.n_neurons);
            sheet_data{row, 22} = std(fr_US_light_baseline)/sqrt(res.n_neurons);
            sheet_data{row, 23} = std(fr_US_light_illum)/sqrt(res.n_neurons);
            sheet_data{row, 24} = std(fr_US_light_resp)/sqrt(res.n_neurons);

            row = row + 3;
        end
    end

    writecell(sheet_data, output_filename, 'Sheet', 'Complete_Neuron_Lists');
end

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
