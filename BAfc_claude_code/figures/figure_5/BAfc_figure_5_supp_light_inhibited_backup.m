% BAfc_figure_5_supp_light_inhibited.m
% Supplementary figure: Light-inhibited neurons during optogenetic manipulation
% Identifies and visualizes neurons inhibited by pre-stimulus light
% Similar layout to BAfc_figure_3.m with heatmaps for LA and Astria

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

% Light inhibition detection parameters (matching BAfc_figure_2.m)
g.light_inhib_window_recent = 0.5;  % 0-0.5s before stimulus (recent)
g.light_inhib_window_baseline = 4.5;  % 4.5s baseline period
g.light_inhib_p_threshold = 0.05;  % p-value threshold for Wilcoxon
g.light_inhib_fr_drop = 0.5;  % 50% FR drop threshold (same as BAfc_figure_2)

% Brain regions to analyze
brain_regions = {'LA', 'Astria'};

fprintf('Identifying light-inhibited neurons...\n');

%% Calculate PSTHs for all conditions
fprintf('Calculating PSTHs...\n');
psthZ_nolight = cell(1, 2);  % CS, US (no-light)
psthHz_nolight = cell(1, 2);
psthZ_light = cell(1, 2);  % CS, US (light)
psthHz_light = cell(1, 2);

% Process no-light conditions
for hmp = 1:2  % CS, US
    psth_spx_og = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl_nolight{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Convert to Hz
    num_trials = size(cell_metrics.general.(ttl_nolight{hmp}){1}, 1);
    psth_hz = psth_spx_og / (num_trials * g.bin_time);
    psth_hz_smooth = smoothdata(psth_hz, 2, 'sgolay', g.smoothvalue);
    psthHz_nolight{hmp} = psth_hz_smooth;

    % Z-score using baseline period only
    baseline_idx = 1:(g.pre_time / g.bin_time);
    baseline_mean = mean(psth_spx_og(:, baseline_idx), 2);
    baseline_std = std(psth_spx_og(:, baseline_idx), 0, 2);
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;
    psth_spx = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);
    psthZ_nolight{hmp} = psth_spx;
end

% Process light conditions
for hmp = 1:2  % CS, US
    psth_spx_og = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl_light{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Convert to Hz
    num_trials = size(cell_metrics.general.(ttl_light{hmp}){1}, 1);
    psth_hz = psth_spx_og / (num_trials * g.bin_time);
    psth_hz_smooth = smoothdata(psth_hz, 2, 'sgolay', g.smoothvalue);
    psthHz_light{hmp} = psth_hz_smooth;

    % Z-score using baseline period only
    baseline_idx = 1:(g.pre_time / g.bin_time);
    baseline_mean = mean(psth_spx_og(:, baseline_idx), 2);
    baseline_std = std(psth_spx_og(:, baseline_idx), 0, 2);
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;
    psth_spx = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);
    psthZ_light{hmp} = psth_spx;
end

%% Identify light-inhibited neurons for each brain region
results_all = cell(1, 2);

for br = 1:2
    fprintf('\nProcessing %s...\n', brain_regions{br});

    % Get ALL neuron indices for this brain region (no cell type filter)
    idx_neurons = strcmp(cell_metrics.brainRegion, brain_regions{br});
    n_neurons = sum(idx_neurons);
    fprintf('  %d total neurons in %s\n', n_neurons, brain_regions{br});

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
        fprintf('  Testing %s for light inhibition...\n', ttl_light{hmp});

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

    fprintf('  Found %d light-inhibited neurons (%.1f%%)\n', n_light_inhibited, n_light_inhibited/n_neurons*100);

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
    fprintf('\nNo light-inhibited neurons found. Exiting.\n');
    return;
end

%% Create figure (4-row layout: LA-CS, LA-US, Astria-CS, Astria-US)
fig = figure('Units', 'pixels', 'Position', [-100, -100, 1200, 1400], 'Visible', 'on');
t = tiledlayout(fig, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

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

%% Plot: 4 rows (LA-CS, LA-US, Astria-CS, Astria-US) × 4 columns (nolight heatmap, light heatmap, nolight lineplot, light lineplot)
for br = 1:2
    if isempty(results_all{br})
        continue;
    end

    res = results_all{br};

    % Get cluster boundaries
    Clusters_sorted = res.Clusters(res.leafOrder);
    n_clu = find(diff(Clusters_sorted) ~= 0);

    % Plot CS and US separately (each gets its own row)
    for stim = 1:2  % CS, US
        % Calculate row index: LA-CS=1, LA-US=2, Astria-CS=3, Astria-US=4
        row_idx = (br-1)*2 + stim;

        % Get PSTHs for this stimulus
        if stim == 1  % CS
            psth_nolight_z = res.psth_CS_nolight(res.leafOrder, :);
            psth_light_z = res.psth_CS_light(res.leafOrder, :);
            psth_nolight_hz = res.psth_CS_Hz_nolight(res.leafOrder, :);
            psth_light_hz = res.psth_CS_Hz_light(res.leafOrder, :);
            onset_lat_nolight = res.CS_onset_lat_nolight;
            offset_lat_nolight = res.CS_offset_lat_nolight;
            onset_lat_light = res.CS_onset_lat_light;
            offset_lat_light = res.CS_offset_lat_light;
            stim_title = 'CS';
        else  % US
            psth_nolight_z = res.psth_US_nolight(res.leafOrder, :);
            psth_light_z = res.psth_US_light(res.leafOrder, :);
            psth_nolight_hz = res.psth_US_Hz_nolight(res.leafOrder, :);
            psth_light_hz = res.psth_US_Hz_light(res.leafOrder, :);
            onset_lat_nolight = res.US_onset_lat_nolight;
            offset_lat_nolight = res.US_offset_lat_nolight;
            onset_lat_light = res.US_onset_lat_light;
            offset_lat_light = res.US_offset_lat_light;
            stim_title = 'US';
        end

        %% Column 1: No-light heatmap
        tile_idx = (row_idx-1)*4 + 1;
        ax = nexttile(t, tile_idx);
        matrix = psth_nolight_z(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        imagesc(g.timeaxis_hmp, 1:size(matrix, 1), matrix);
        clim([clim_min clim_max]);
        colormap(ax, g.colors.Heatmap);
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        % Labels and title
        if row_idx == 1
            title(sprintf('%s - No light', stim_title), 'FontSize', g.fontSize1);
        end
        if stim == 1
            ylabel(sprintf('%s\nNeuron #', brain_regions{br}), 'FontSize', g.fontSize2);
        end

        yticks([1, size(matrix, 1)]);
        if row_idx == 4
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
        end
        set(gca, 'FontSize', g.fontSize2);

        % Add cluster lines
        hold on;
        for i = 1:length(n_clu)
            yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
        end

        % Add onset/offset markers (no-light)
        for n = 1:length(res.leafOrder)
            idx_n = res.leafOrder(n);
            if ~isnan(onset_lat_nolight(idx_n))
                plot(onset_lat_nolight(idx_n), n, 'k.', 'MarkerSize', 4);
            end
            if ~isnan(offset_lat_nolight(idx_n))
                plot(offset_lat_nolight(idx_n), n, 'k.', 'MarkerSize', 4);
            end
        end
        hold off;

        %% Column 2: Light heatmap
        tile_idx = (row_idx-1)*4 + 2;
        ax = nexttile(t, tile_idx);
        matrix = psth_light_z(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        imagesc(g.timeaxis_hmp, 1:size(matrix, 1), matrix);
        clim([clim_min clim_max]);
        colormap(ax, g.colors.Heatmap);
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        % Labels and title
        if row_idx == 1
            title(sprintf('%s - Light', stim_title), 'FontSize', g.fontSize1);
        end

        yticks([1, size(matrix, 1)]);
        set(gca, 'YTickLabel', []);
        if row_idx == 4
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
        end
        set(gca, 'FontSize', g.fontSize2);

        % Add cluster lines
        hold on;
        for i = 1:length(n_clu)
            yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
        end

        % Add onset/offset markers (light)
        for n = 1:length(res.leafOrder)
            idx_n = res.leafOrder(n);
            if ~isnan(onset_lat_light(idx_n))
                plot(onset_lat_light(idx_n), n, 'k.', 'MarkerSize', 4);
            end
            if ~isnan(offset_lat_light(idx_n))
                plot(offset_lat_light(idx_n), n, 'k.', 'MarkerSize', 4);
            end
        end
        hold off;

        %% Column 3: No-light lineplot
        tile_idx = (row_idx-1)*4 + 3;
        ax = nexttile(t, tile_idx);
        plot_idx = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        time_vec = g.timeaxis_hmp;
        if length(time_vec) ~= length(plot_idx)
            time_vec = linspace(-g.plotwin(1), g.plotwin(2), length(plot_idx));
        end

        psth_mean = mean(psth_nolight_hz(:, plot_idx), 1);
        hold on;
        xline(0, '--k', 'LineWidth', g.xlinewidth);
        plot(time_vec, psth_mean, '-', 'Color', [0.2 0.6 0.2], 'LineWidth', 2.5);
        hold off;

        xlim([time_vec(1) time_vec(end)]);
        set(gca, 'YTickLabel', []);
        if row_idx == 4
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
        end
        set(gca, 'FontSize', g.fontSize2);

        %% Column 4: Light lineplot
        tile_idx = (row_idx-1)*4 + 4;
        ax = nexttile(t, tile_idx);

        psth_mean = mean(psth_light_hz(:, plot_idx), 1);
        hold on;
        xline(0, '--k', 'LineWidth', g.xlinewidth);
        plot(time_vec, psth_mean, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2.5);
        hold off;

        xlim([time_vec(1) time_vec(end)]);
        set(gca, 'YTickLabel', []);
        if row_idx == 4
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
cb_width = 0.010;
cb_left = 0.025;
cb_bottom = 0.40;
cb_height = 0.20;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'left');
ylabel(cb_ax, 'Z-score', 'FontSize', g.fontSize2);
cb_ax.FontSize = g.fontSize2;

fprintf('\nDone. Light-inhibited neurons visualized.\n');

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
