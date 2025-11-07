% BAfc_figure_4.m
% Monosynaptic responses for LA and Astria
% Same layout as figure 3 but with monosynaptic detection (0-25ms window)

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
g.fontSize1 = 10;
g.fontSize2 = 10;
g.bin_time = 0.001;
g.pre_time = 5;
g.post_time = 0.5;
g.monosyn_window = 0.025;  % 0-25ms
g.smoothvalue = 7;
g.plotwin = [0.05 0.05];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.xlinewidth = 2;
g.clim_percentile = 95;

g.onset_threshold = 3;  % Lower threshold for monosynaptic detection
g.min_consec_bins = max(1, round(0.001 / g.bin_time));  % 3ms for brief monosynaptic responses
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

cluster_colors = [0.8 0.2 0.2; 0.2 0.4 0.8; 0.6 0.2 0.6];  % CS-sel, US-sel, Multi

%% Calculate PSTHs once
fprintf('Calculating PSTHs...\n');
psthZ_full = cell(1,3);
psthHz_full = cell(1,3);
baseline_bins = round(g.pre_time / g.bin_time);
for hmp = 1:3
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

    n_neurons = sum(idx_neurons);

    if n_neurons == 0
        continue;
    end

    % Get number of trials
    num_trials_CS = size(cell_metrics.general.(ttl{1}){1}, 1);
    num_trials_US = size(cell_metrics.general.(ttl{2}){1}, 1);
    num_trials_Both = size(cell_metrics.general.(ttl{3}){1}, 1);

    % Extract monosynaptic window responses for CS, US, CS+US
    % Pool all neurons together
    psth_CS_all = psthZ_full{1}(idx_neurons, :);
    psth_US_all = psthZ_full{2}(idx_neurons, :);
    psth_Both_all = psthZ_full{3}(idx_neurons, :);
    psth_CS_Hz_all = psthHz_full{1}(idx_neurons, :);
    psth_US_Hz_all = psthHz_full{2}(idx_neurons, :);
    psth_Both_Hz_all = psthHz_full{3}(idx_neurons, :);

    % Calculate spike probabilities trial-by-trial
    [~, ~, postAP_norm_CS] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{1}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    [~, ~, postAP_norm_US] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{2}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Process all neurons (pooled)
    fprintf('  Processing all neurons (%d neurons)...\n', n_neurons);
    CS_peak_all = max(psth_CS_all(:, monosyn_window_bins), [], 2);
    US_peak_all = max(psth_US_all(:, monosyn_window_bins), [], 2);

    neuron_indices_all = find(idx_neurons);
    CS_prob_all = zeros(n_neurons, 1);
    US_prob_all = zeros(n_neurons, 1);

    for n = 1:n_neurons
        global_idx = neuron_indices_all(n);

        % CS probability - count trials with at least 1 spike in 0-25ms window
        if ~isempty(postAP_norm_CS{global_idx})
            responsive_trials_CS = 0;
            for trial = 1:length(postAP_norm_CS{global_idx})
                trial_spikes = postAP_norm_CS{global_idx}{trial};
                if any(trial_spikes > 0 & trial_spikes <= g.monosyn_window)
                    responsive_trials_CS = responsive_trials_CS + 1;
                end
            end
            CS_prob_all(n) = responsive_trials_CS / num_trials_CS;
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
            US_prob_all(n) = responsive_trials_US / num_trials_US;
        end
    end

    % Responsiveness for all neurons
    if g.use_two_rule
        % Two-rule: (Rule 1 OR Rule 2)
        CS_responsive_all = (CS_peak_all >= g.zscore_threshold_rule1 & CS_prob_all >= g.prob_threshold_rule1) | ...
                            (CS_peak_all >= g.zscore_threshold_rule2 & CS_prob_all >= g.prob_threshold_rule2);
        US_responsive_all = (US_peak_all >= g.zscore_threshold_rule1 & US_prob_all >= g.prob_threshold_rule1) | ...
                            (US_peak_all >= g.zscore_threshold_rule2 & US_prob_all >= g.prob_threshold_rule2);
    else
        % One-rule: z-score only
        CS_responsive_all = CS_peak_all >= g.zscore_threshold_one_rule;
        US_responsive_all = US_peak_all >= g.zscore_threshold_one_rule;
    end

    fprintf('  All neurons - CS responsive: %d, US responsive: %d\n', sum(CS_responsive_all), sum(US_responsive_all));

    % Categorize all neurons
    Clusters_all = zeros(n_neurons, 1);
    Clusters_all(CS_responsive_all & ~US_responsive_all) = 1;  % CS-selective
    Clusters_all(US_responsive_all & ~CS_responsive_all) = 2;  % US-selective
    Clusters_all(CS_responsive_all & US_responsive_all) = 3;   % Multisensory

    % Calculate onset/offset latencies for sorting (using helper function)
    CS_onset_lat_all = nan(n_neurons, 1);
    CS_offset_lat_all = nan(n_neurons, 1);
    US_onset_lat_all = nan(n_neurons, 1);
    US_offset_lat_all = nan(n_neurons, 1);
    Both_onset_lat_all = nan(n_neurons, 1);
    Both_offset_lat_all = nan(n_neurons, 1);

    for n = 1:n_neurons
        [CS_onset_lat_all(n), CS_offset_lat_all(n)] = compute_onset_offset_latency(psth_CS_all(n, :), monosyn_window_bins, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [US_onset_lat_all(n), US_offset_lat_all(n)] = compute_onset_offset_latency(psth_US_all(n, :), monosyn_window_bins, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [Both_onset_lat_all(n), Both_offset_lat_all(n)] = compute_onset_offset_latency(psth_Both_all(n, :), monosyn_window_bins, g.onset_threshold, g.min_consec_bins, g.bin_time);
    end

    fprintf('  All neurons - CS onset latencies (ms): valid=%d\n', sum(~isnan(CS_onset_lat_all)));
    fprintf('  All neurons - US onset latencies (ms): valid=%d\n', sum(~isnan(US_onset_lat_all)));

    % Sort all neurons by category using rank score (onset + alpha * duration)
    leafOrder_all = [];
    for c = [1 2 3]
        clust_idx = find(Clusters_all == c);
        if isempty(clust_idx)
            continue;
        end

        if c == 1  % CS-selective
            onset_c = CS_onset_lat_all(clust_idx);
            offset_c = CS_offset_lat_all(clust_idx);
        elseif c == 2  % US-selective
            onset_c = US_onset_lat_all(clust_idx);
            offset_c = US_offset_lat_all(clust_idx);
        else  % Multisensory - use CS rank score
            onset_c = CS_onset_lat_all(clust_idx);
            offset_c = CS_offset_lat_all(clust_idx);
        end

        duration_c = offset_c - onset_c;
        rank_score = onset_c + g.alpha * duration_c;
        sort_matrix = [isnan(rank_score), rank_score];
        [~, sort_idx] = sortrows(sort_matrix, [1 2]);
        leafOrder_all = [leafOrder_all; clust_idx(sort_idx)];
    end

    % Store results
    results_all{br}.Clusters_all = Clusters_all;
    results_all{br}.leafOrder_all = leafOrder_all;
    results_all{br}.psth_CS_all = psth_CS_all;
    results_all{br}.psth_US_all = psth_US_all;
    results_all{br}.psth_Both_all = psth_Both_all;
    results_all{br}.psth_CS_Hz_all = psth_CS_Hz_all;
    results_all{br}.psth_US_Hz_all = psth_US_Hz_all;
    results_all{br}.psth_Both_Hz_all = psth_Both_Hz_all;
    results_all{br}.CS_onset_lat_all = CS_onset_lat_all;
    results_all{br}.CS_offset_lat_all = CS_offset_lat_all;
    results_all{br}.US_onset_lat_all = US_onset_lat_all;
    results_all{br}.US_offset_lat_all = US_offset_lat_all;
    results_all{br}.Both_onset_lat_all = Both_onset_lat_all;
    results_all{br}.Both_offset_lat_all = Both_offset_lat_all;
    results_all{br}.n_neurons = n_neurons;

    fprintf('  All neurons: %d | CS-sel: %d | US-sel: %d | Multi: %d\n', ...
        n_neurons, sum(Clusters_all==1), sum(Clusters_all==2), sum(Clusters_all==3));
end

%% Create figure - heatmaps + bar charts + latency boxplots
fig = figure('Units', 'pixels', 'Position', [100, 100, 1000, 600], 'Visible', 'on');
t = tiledlayout(fig, 4, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

% Determine global color limits across all heatmaps (CS, US, and CS+US)
all_values = [];
for br = 1:2
    if ~isempty(results_all{br})
        res = results_all{br};

        % All neurons pooled
        Clusters_sorted_all = res.Clusters_all(res.leafOrder_all);
        responsive_mask_all = ismember(Clusters_sorted_all, [1, 2, 3]);
        responsive_leafOrder_all = res.leafOrder_all(responsive_mask_all);

        for stim = 1:3
            if stim == 1
                psth_sorted_all = res.psth_CS_all(responsive_leafOrder_all, :);
            elseif stim == 2
                psth_sorted_all = res.psth_US_all(responsive_leafOrder_all, :);
            else
                psth_sorted_all = res.psth_Both_all(responsive_leafOrder_all, :);
            end
            plot_bins = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
            matrix_all = psth_sorted_all(:, plot_bins);
            all_values = [all_values; matrix_all(:)];
        end
    end
end

clim_min = prctile(all_values, (100 - g.clim_percentile) / 2);
clim_max = prctile(all_values, 100 - (100 - g.clim_percentile) / 2);
fprintf('Global color limits: [%.2f, %.2f]\n', clim_min, clim_max);

% Create nested tiledlayout for all 6 heatmaps (2 brain regions × 3 stimuli)
% This spans rows 1-4, columns 1-3 of the main layout
t_heatmaps = tiledlayout(t, 2, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
t_heatmaps.Layout.Tile = 1;  % Start at tile 1 (row 1, column 1)
t_heatmaps.Layout.TileSpan = [4 3];  % Span 4 rows, 3 columns

for br = 1:2
    if isempty(results_all{br})
        continue;
    end

    res = results_all{br};

    % All neurons pooled
    Clusters_sorted_all = res.Clusters_all(res.leafOrder_all);
    responsive_mask_all = ismember(Clusters_sorted_all, [1, 2, 3]);
    responsive_leafOrder_all = res.leafOrder_all(responsive_mask_all);
    Clusters_sorted_responsive_all = Clusters_sorted_all(responsive_mask_all);
    n_clu_all = find(diff(Clusters_sorted_responsive_all) ~= 0);

    % Heatmaps for CS, US, CS+US (columns 1-3)
    for stim = 1:3
        if stim == 1
            psth_sorted_all = res.psth_CS_all(responsive_leafOrder_all, :);
            stim_title = 'CS';
        elseif stim == 2
            psth_sorted_all = res.psth_US_all(responsive_leafOrder_all, :);
            stim_title = 'US';
        else
            psth_sorted_all = res.psth_Both_all(responsive_leafOrder_all, :);
            stim_title = 'CS+US';
        end

        n_responsive = sum(responsive_mask_all);

        % Heatmaps in nested tiledlayout (2 rows × 3 columns)
        % Row 1 (br=1, LA): tiles 1, 2, 3
        % Row 2 (br=2, Astria): tiles 4, 5, 6
        tile_num = (br-1)*3 + stim;
        ax = nexttile(t_heatmaps, tile_num);
        plot_bins = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        matrix = psth_sorted_all(:, plot_bins);
        imagesc(g.timeaxis_hmp, 1:size(matrix,1), matrix);
        clim([clim_min clim_max]);
        colormap(ax, g.colors.Heatmap);
        hold on;

        % Indicate artifact period (0-10ms) with semi-transparent gray overlay
        if stim >= 1  % All stimuli (CS, US, CS+US)
            patch([0 0.010 0.010 0], [0.5 0.5 size(matrix,1)+0.5 size(matrix,1)+0.5], ...
                'k', 'FaceAlpha', 0.8, 'EdgeColor', 'none');

            % Add lightning bolt symbol at top of artifact region (only for US and CS+US)
            if stim == 2 || stim == 3
                text(0.005, 2, '⚡', 'Color', 'y', 'FontSize', g.fontSize1, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            end
        end

        xline(0, '--k', 'LineWidth', g.xlinewidth);
        xline(g.monosyn_window, '-k', 'LineWidth', 1);  % End of monosyn window

        % Labels
        if br == 1
            title(stim_title, 'FontSize', g.fontSize1, 'FontWeight', 'bold');
        end

        if stim == 1
            % Display AStria instead of Astria
            if strcmp(brain_regions{br}, 'Astria')
                ylabel('AStria', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
            else
                ylabel(brain_regions{br}, 'FontSize', g.fontSize2, 'FontWeight', 'bold');
            end
        else
            set(gca, 'YTickLabel', []);
        end

        if br == 2
            xlabel('Time (s)', 'FontSize', g.fontSize2);
        else
            set(gca, 'XTickLabel', []);
        end

        if size(matrix,1) > 1
            yticks([1, size(matrix,1)]);
        else
            yticks(1);
        end

        % Cluster boundaries (black lines)
        for i = 1:length(n_clu_all)
            yline(n_clu_all(i) + 0.5, 'k-', 'LineWidth', 1);
        end

        % Add onset and offset markers for all neurons
        for n = 1:n_responsive
            idx_n = responsive_leafOrder_all(n);

            % Get appropriate onset/offset latencies based on stimulus
            if stim == 1  % CS
                onset_lat = res.CS_onset_lat_all(idx_n);
                offset_lat = res.CS_offset_lat_all(idx_n);
            elseif stim == 2  % US
                onset_lat = res.US_onset_lat_all(idx_n);
                offset_lat = res.US_offset_lat_all(idx_n);
            else  % CS+US
                onset_lat = res.Both_onset_lat_all(idx_n);
                offset_lat = res.Both_offset_lat_all(idx_n);
            end

            % Plot onset marker (black dot)
            if ~isnan(onset_lat)
                plot(onset_lat, n, 'k.', 'MarkerSize', 4);
            end

            % Plot offset marker (black dot)
            if ~isnan(offset_lat)
                plot(offset_lat, n, 'k.', 'MarkerSize', 4);
            end
        end

        hold off;

        set(gca, 'FontSize', g.fontSize2);
    end
end

% Add colorbar - manually create with full control, positioned at bottom heatmap (Astria)
drawnow;  % Ensure all positions are updated
cb_width = 0.010;
cb_left = 0.033;
cb_bottom = 0.04;  % Aligned with Astria heatmap (row 2 in 4-row layout)
cb_height = 0.15;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'left');
ylabel(cb_ax, 'Z-score', 'FontSize', g.fontSize2);
cb_ax.FontSize = g.fontSize2;

%% Add metric bar charts (3 clusters × 2 metrics in columns 4-5)
cluster_names = {'CS-sel', 'US-sel', 'Multi'};
metric_names = {'\DeltaMean FR (Hz)', '\DeltaPeak FR (Hz)'};

% Create nested tiledlayout for bar charts spanning rows 1-3, columns 4-5
t_bars = tiledlayout(t, 3, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
t_bars.Layout.Tile = 4;  % Start at tile 4 (row 1, column 4)
t_bars.Layout.TileSpan = [3 2];  % Span 3 rows, 2 columns

for c = [1 2 3]  % CS-selective, US-selective, Multisensory
    for metric = 1:2  % Mean FR, Peak FR only
        % Bar charts in nested layout: 3 rows × 2 columns
        ax_metric = nexttile(t_bars, (c-1)*2 + metric);

        % Collect data for LA and Astria (all neurons pooled)
        data_LA = [];
        data_Astria = [];

        for br = 1:2
            if isempty(results_all{br})
                continue;
            end

            res = results_all{br};
            clust_idx_all = find(res.Clusters_all == c);

            % Process all neurons
            if ~isempty(clust_idx_all)
                if metric == 1  % Mean firing rate change (Hz) in monosyn window
                    CS_metric_all = zeros(length(clust_idx_all), 1);
                    US_metric_all = zeros(length(clust_idx_all), 1);
                    Both_metric_all = zeros(length(clust_idx_all), 1);

                    baseline_idx = 1:baseline_bins;
                    % Define response window (12ms to g.monosyn_window)
                    response_start_bin = round((g.pre_time + 0.012) / g.bin_time);
                    response_end_bin = round((g.pre_time + g.monosyn_window) / g.bin_time);

                    for n = 1:length(clust_idx_all)
                        idx_n = clust_idx_all(n);
                        baseline_fr = mean(res.psth_CS_Hz_all(idx_n, baseline_idx));

                        % Calculate mean FR in 12-25ms window for all neurons
                        response_fr_CS = mean(res.psth_CS_Hz_all(idx_n, response_start_bin:response_end_bin));
                        CS_metric_all(n) = response_fr_CS - baseline_fr;

                        response_fr_US = mean(res.psth_US_Hz_all(idx_n, response_start_bin:response_end_bin));
                        US_metric_all(n) = response_fr_US - baseline_fr;

                        response_fr_Both = mean(res.psth_Both_Hz_all(idx_n, response_start_bin:response_end_bin));
                        Both_metric_all(n) = response_fr_Both - baseline_fr;
                    end

                else  % metric == 2: Peak firing rate change (Hz) in monosyn window
                    CS_metric_all = zeros(length(clust_idx_all), 1);
                    US_metric_all = zeros(length(clust_idx_all), 1);
                    Both_metric_all = zeros(length(clust_idx_all), 1);

                    baseline_idx = 1:baseline_bins;
                    % Define response window (12ms to g.monosyn_window)
                    response_start_bin = round((g.pre_time + 0.012) / g.bin_time);
                    response_end_bin = round((g.pre_time + g.monosyn_window) / g.bin_time);

                    for n = 1:length(clust_idx_all)
                        idx_n = clust_idx_all(n);
                        baseline_fr = mean(res.psth_CS_Hz_all(idx_n, baseline_idx));

                        % Calculate peak FR in 12-25ms window for all neurons
                        peak_fr_CS = max(res.psth_CS_Hz_all(idx_n, response_start_bin:response_end_bin));
                        CS_metric_all(n) = peak_fr_CS - baseline_fr;

                        peak_fr_US = max(res.psth_US_Hz_all(idx_n, response_start_bin:response_end_bin));
                        US_metric_all(n) = peak_fr_US - baseline_fr;

                        peak_fr_Both = max(res.psth_Both_Hz_all(idx_n, response_start_bin:response_end_bin));
                        Both_metric_all(n) = peak_fr_Both - baseline_fr;
                    end
                end

                if br == 1
                    data_LA = [CS_metric_all, US_metric_all, Both_metric_all];
                else
                    data_Astria = [CS_metric_all, US_metric_all, Both_metric_all];
                end
            end
        end

        % Plot grouped bars - LA [2 3 4], Astria [7 8 9]
        hold on;

        % LA bars at positions 2, 3, 4 - darker cluster color
        if ~isempty(data_LA)
            means_LA = mean(data_LA, 1, 'omitnan');
            sems_LA = std(data_LA, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_LA(:,1))));
            bar_color_LA = cluster_colors(c, :) * 0.7;  % Darker
            bar([2 3 4], means_LA, 0.6, 'FaceColor', bar_color_LA, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([2 3 4], means_LA, sems_LA, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
        end

        % Astria bars at positions 7, 8, 9 - lighter cluster color
        if ~isempty(data_Astria)
            means_Astria = mean(data_Astria, 1, 'omitnan');
            sems_Astria = std(data_Astria, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_Astria(:,1))));
            bar_color_Astria = cluster_colors(c, :) + (1 - cluster_colors(c, :)) * 0.5;  % Lighter
            bar([7 8 9], means_Astria, 0.6, 'FaceColor', bar_color_Astria, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([7 8 9], means_Astria, sems_Astria, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
        end

        hold off;

        % Statistical comparisons - Wilcoxon signed rank test (paired data)
        % LA comparisons
        if ~isempty(data_LA) && sum(~isnan(data_LA(:,1))) > 0 && sum(~isnan(data_LA(:,2))) > 0 && sum(~isnan(data_LA(:,3))) > 0
            y_max_LA = max(means_LA + sems_LA);
            hold on;

            % Level 1: Short comparisons
            y_pos_level1 = y_max_LA + y_max_LA * 0.08;

            % CS vs US
            try
                [p_LA_CS_US, ~] = signrank(data_LA(:,1), data_LA(:,2));
                if p_LA_CS_US < 0.05
                    sig_text = get_sig_stars(p_LA_CS_US);
                    plot([2 3], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                    text(2.5, y_pos_level1, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
                end
            catch
            end

            % US vs Both
            try
                [p_LA_US_Both, ~] = signrank(data_LA(:,2), data_LA(:,3));
                if p_LA_US_Both < 0.05
                    sig_text = get_sig_stars(p_LA_US_Both);
                    plot([3 4], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                    text(3.5, y_pos_level1, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
                end
            catch
            end

            % Level 2: Long comparison
            y_pos_level2 = y_max_LA + y_max_LA * 0.25;

            % CS vs Both
            try
                [p_LA_CS_Both, ~] = signrank(data_LA(:,1), data_LA(:,3));
                if p_LA_CS_Both < 0.05
                    sig_text = get_sig_stars(p_LA_CS_Both);
                    plot([2 4], [y_pos_level2 y_pos_level2], 'k-', 'LineWidth', 1.5);
                    text(3, y_pos_level2, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
                end
            catch
            end

            hold off;
        end

        % Astria comparisons
        if ~isempty(data_Astria) && sum(~isnan(data_Astria(:,1))) > 0 && sum(~isnan(data_Astria(:,2))) > 0 && sum(~isnan(data_Astria(:,3))) > 0
            y_max_Astria = max(means_Astria + sems_Astria);
            hold on;

            % Level 1: Short comparisons
            y_pos_level1 = y_max_Astria + y_max_Astria * 0.08;

            % CS vs US
            try
                [p_Astria_CS_US, ~] = signrank(data_Astria(:,1), data_Astria(:,2));
                if p_Astria_CS_US < 0.05
                    sig_text = get_sig_stars(p_Astria_CS_US);
                    plot([7 8], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                    text(7.5, y_pos_level1, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
                end
            catch
            end

            % US vs Both
            try
                [p_Astria_US_Both, ~] = signrank(data_Astria(:,2), data_Astria(:,3));
                if p_Astria_US_Both < 0.05
                    sig_text = get_sig_stars(p_Astria_US_Both);
                    plot([8 9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                    text(8.5, y_pos_level1, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
                end
            catch
            end

            % Level 2: Long comparison
            y_pos_level2 = y_max_Astria + y_max_Astria * 0.25;

            % CS vs Both
            try
                [p_Astria_CS_Both, ~] = signrank(data_Astria(:,1), data_Astria(:,3));
                if p_Astria_CS_Both < 0.05
                    sig_text = get_sig_stars(p_Astria_CS_Both);
                    plot([7 9], [y_pos_level2 y_pos_level2], 'k-', 'LineWidth', 1.5);
                    text(8, y_pos_level2, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
                end
            catch
            end

            hold off;
        end

        % Formatting
        xlim([0.5 10.5]);
        xticks([2 3 4 7 8 9]);

        % Set fixed y-limits based on metric
        if metric == 1  % Delta Mean FR
            ylim([0 200]);
            y_ticks = [0, 100, 200];
        else  % Delta Peak FR
            ylim([0 300]);
            y_ticks = [0, 150, 300];
        end

        if c == 1
            title(metric_names{metric}, 'FontSize', g.fontSize1, 'FontWeight', 'bold');
        end

        % X-tick labels only on bottom row
        if c == 3
            xticklabels({'CS', 'US', 'CS+US', 'CS', 'US', 'CS+US'});
        else
            set(gca, 'XTickLabel', []);
        end

        % Set y-ticks
        yticks(y_ticks);

        % Add LA/Astria labels on top row
        if c == 1
            curr_ylim_extended = ylim;
            y_range = curr_ylim_extended(2) - curr_ylim_extended(1);
            if metric == 1
                y_pos = curr_ylim_extended(2) - 0.05*y_range;
            else
                y_pos = curr_ylim_extended(2) - 0.10*y_range;
            end
            text(3, y_pos, 'LA', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            text(8, y_pos, 'AStria', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end

        % Add y-axis label
        if metric == 1
            ylabel(cluster_names{c}, 'FontSize', g.fontSize2, 'FontWeight', 'bold');
        end

        set(gca, 'FontSize', g.fontSize2);
        box off;
    end
end

%% Add latency boxplot panel (spanning columns 4-5, row 4) - Brain region comparison
ax_latency = nexttile(t, 19, [1 2]);  % Row 4, columns 4-5

% Collect latency data - POOLING selective and multisensory together for each region
% LA CS latency (CS-selective + Multisensory, using CS latency)
LA_CS_pooled_lat = [results_all{1}.CS_onset_lat_all(results_all{1}.Clusters_all == 1); ...
                    results_all{1}.CS_onset_lat_all(results_all{1}.Clusters_all == 3)] * 1000;  % Convert to ms

% Astria CS latency (CS-selective + Multisensory, using CS latency)
Astria_CS_pooled_lat = [results_all{2}.CS_onset_lat_all(results_all{2}.Clusters_all == 1); ...
                        results_all{2}.CS_onset_lat_all(results_all{2}.Clusters_all == 3)] * 1000;

% LA US latency (US-selective + Multisensory, using US latency)
LA_US_pooled_lat = [results_all{1}.US_onset_lat_all(results_all{1}.Clusters_all == 2); ...
                    results_all{1}.US_onset_lat_all(results_all{1}.Clusters_all == 3)] * 1000;

% Astria US latency (US-selective + Multisensory, using US latency)
Astria_US_pooled_lat = [results_all{2}.US_onset_lat_all(results_all{2}.Clusters_all == 2); ...
                        results_all{2}.US_onset_lat_all(results_all{2}.Clusters_all == 3)] * 1000;

% Remove NaN values
LA_CS_pooled_lat = LA_CS_pooled_lat(~isnan(LA_CS_pooled_lat));
Astria_CS_pooled_lat = Astria_CS_pooled_lat(~isnan(Astria_CS_pooled_lat));
LA_US_pooled_lat = LA_US_pooled_lat(~isnan(LA_US_pooled_lat));
Astria_US_pooled_lat = Astria_US_pooled_lat(~isnan(Astria_US_pooled_lat));

hold on;

% Positions for 2 pairs of boxplots (LA vs Astria for CS and US)
positions = [1 2, 5 6];
all_data = {LA_CS_pooled_lat, Astria_CS_pooled_lat, LA_US_pooled_lat, Astria_US_pooled_lat};

% Plot boxplots
boxplot(ax_latency, [all_data{1}; all_data{2}; all_data{3}; all_data{4}], ...
        [ones(length(all_data{1}),1)*1; ones(length(all_data{2}),1)*2; ...
         ones(length(all_data{3}),1)*5; ones(length(all_data{4}),1)*6], ...
        'Positions', positions, 'Widths', 0.6, 'Colors', 'k');

% Overlay scatter points with jitter - color by brain region
jitter_amount = 0.15;
region_colors = {[0.8 0.2 0.2], [0.2 0.4 0.8]};  % Red for LA, Blue for Astria
for i = 1:4
    if ~isempty(all_data{i})
        x_jitter = positions(i) + (rand(length(all_data{i}), 1) - 0.5) * jitter_amount;

        % Color based on region: LA (odd) vs Astria (even)
        if mod(i, 2) == 1  % LA
            color = region_colors{1};
        else  % Astria
            color = region_colors{2};
        end

        scatter(x_jitter, all_data{i}, 20, color, 'filled', 'MarkerFaceAlpha', 0.4);
    end
end

hold off;

% Formatting
xlim([0 7]);
ylim([0 50]);
xticks([1.5 5.5]);
xticklabels({'CS', 'US'});
ylabel('Onset Latency (ms)', 'FontSize', g.fontSize2);
title('Onset Latencies', 'FontSize', g.fontSize1, 'FontWeight', 'bold');
set(gca, 'FontSize', g.fontSize2);
box off;

% Add statistics - ranksum test for all comparisons
hold on;

% Level 1 (y = 45): LA vs Astria within each stimulus (short lines)
y_pos_level1 = 45;
pairs_level1 = {
    {LA_CS_pooled_lat, Astria_CS_pooled_lat, 1, 2};  % CS: LA vs Astria
    {LA_US_pooled_lat, Astria_US_pooled_lat, 5, 6};  % US: LA vs Astria
};

for p = 1:2
    data1 = pairs_level1{p}{1};
    data2 = pairs_level1{p}{2};
    x_left = pairs_level1{p}{3};
    x_right = pairs_level1{p}{4};

    if ~isempty(data1) && ~isempty(data2) && length(data1) > 0 && length(data2) > 0
        [p_val, ~] = ranksum(data1, data2);

        % Draw line connecting the two boxplots
        plot([x_left x_right], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);

        if p_val < 0.05
            sig_text = get_sig_stars(p_val);
        else
            sig_text = 'n.s.';
        end

        text((x_left + x_right) / 2, y_pos_level1, sig_text, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
    end
end

% Level 2 (y = 48): CS vs US within each region (longer lines)
y_pos_level2 = 48;
pairs_level2 = {
    {LA_CS_pooled_lat, LA_US_pooled_lat, 1, 5};      % LA: CS vs US
    {Astria_CS_pooled_lat, Astria_US_pooled_lat, 2, 6};  % Astria: CS vs US
};

for p = 1:2
    data1 = pairs_level2{p}{1};
    data2 = pairs_level2{p}{2};
    x_left = pairs_level2{p}{3};
    x_right = pairs_level2{p}{4};

    if ~isempty(data1) && ~isempty(data2) && length(data1) > 0 && length(data2) > 0
        [p_val, ~] = ranksum(data1, data2);

        % Draw line connecting the two boxplots
        plot([x_left x_right], [y_pos_level2 y_pos_level2], 'k-', 'LineWidth', 1.5);

        if p_val < 0.05
            sig_text = get_sig_stars(p_val);
        else
            sig_text = 'n.s.';
        end

        text((x_left + x_right) / 2, y_pos_level2, sig_text, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
    end
end

hold off;

% Add legend for brain regions
legend_ax = axes('Position', [0.85 0.15 0.1 0.05], 'Visible', 'off');
hold(legend_ax, 'on');
scatter(legend_ax, 1, 1, 50, region_colors{1}, 'filled');
scatter(legend_ax, 1, 2, 50, region_colors{2}, 'filled');
legend(legend_ax, {'LA', 'Astria'}, 'Location', 'best', 'FontSize', g.fontSize2);
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
