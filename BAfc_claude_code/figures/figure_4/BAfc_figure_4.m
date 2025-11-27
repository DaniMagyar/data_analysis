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
g.smoothvalue = 5;
g.plotwin = [0.05 0.05];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.xlinewidth = 2;
g.clim_percentile = 95;

g.onset_threshold = 5;  % Lower threshold for monosynaptic detection
g.min_consec_bins = max(1, round(0.001 / g.bin_time));  % 3ms for brief monosynaptic responses
g.alpha = 0.0;

% Responsiveness detection method
g.use_two_rule = false;  % true: two-rule (Rule 1 OR Rule 2), false: one-rule (z-score only)

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

% Calculate Savitzky-Golay filter delay correction
filter_delay = floor(g.smoothvalue / 2);  % Half the filter width (symmetric kernel)
fprintf('Savitzky-Golay filter width: %d bins, delay correction: %d bins (%.1f ms)\n', ...
    g.smoothvalue, filter_delay, filter_delay * g.bin_time * 1000);

for hmp = 1:3
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
    baseline_std(baseline_std == 0) = 1;  % Avoid division by zero
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;
    psth_spx_smooth = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);

    % Correct for filter delay by shifting forward (move data earlier in time)
    psth_spx_corrected = zeros(size(psth_spx_smooth));
    psth_spx_corrected(:, filter_delay+1:end) = psth_spx_smooth(:, 1:end-filter_delay);
    psth_spx_corrected(:, 1:filter_delay) = repmat(psth_spx_smooth(:, 1), 1, filter_delay);
    psthZ_full{hmp} = psth_spx_corrected;
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
    results_all{br}.animals = cell_metrics.animal(idx_neurons);

    fprintf('  All neurons: %d | CS-sel: %d | US-sel: %d | Multi: %d\n', ...
        n_neurons, sum(Clusters_all==1), sum(Clusters_all==2), sum(Clusters_all==3));
end

%% Create figure - 4x6 layout: heatmaps (cols 1-3) + Delta FR bars (col 4) + remaining (cols 5-6)
fig = figure('Units', 'pixels', 'Position', [100, 100, 1200, 800], 'Visible', 'on');
t = tiledlayout(fig, 4, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add panel labels for main figure
% A: LA heatmaps (row 1, col 1)
annotation(fig, 'textbox', [0.01 0.94 0.05 0.05], 'String', 'A', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% B: LA delta FR plots (row 1, col 4)
annotation(fig, 'textbox', [0.48 0.94 0.05 0.05], 'String', 'B', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% C: AStria heatmaps (row 3, col 1)
annotation(fig, 'textbox', [0.01 0.47 0.05 0.05], 'String', 'C', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% D: AStria delta FR plots (row 3, col 4)
annotation(fig, 'textbox', [0.48 0.47 0.05 0.05], 'String', 'D', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% E: Pie charts (row 1, col 5)
annotation(fig, 'textbox', [0.68 0.94 0.05 0.05], 'String', 'E', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% F: Across region comparison (row 2, col 5)
annotation(fig, 'textbox', [0.68 0.62 0.05 0.05], 'String', 'F', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% G: Onset latency comparison (row 3, col 5)
annotation(fig, 'textbox', [0.68 0.28 0.05 0.05], 'String', 'G', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

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

% Create nested tiledlayout for heatmaps + bar graphs (2 brain regions × 4 columns)
% This spans rows 1-4, columns 1-4 of the main layout
t_heatmaps = tiledlayout(t, 2, 4, 'TileSpacing', 'tight', 'Padding', 'tight');
t_heatmaps.Layout.Tile = 1;  % Start at tile 1 (row 1, column 1)
t_heatmaps.Layout.TileSpan = [4 4];  % Span 4 rows, 4 columns

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

        % Heatmaps in nested tiledlayout (2 rows × 4 columns)
        % Row 1 (br=1, LA): tiles 1, 2, 3 (column 4 is bars)
        % Row 2 (br=2, Astria): tiles 5, 6, 7 (column 8 is bars)
        tile_num = (br-1)*4 + stim;
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
                ylabel('AStria neurons', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
            else
                ylabel([brain_regions{br} ' neurons'], 'FontSize', g.fontSize2, 'FontWeight', 'bold');
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
cb_width = 0.008;
cb_left = 0.185;
cb_bottom = 0.06;  % Aligned with Astria heatmap (row 2 in 4-row layout)
cb_height = 0.08;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'right');
%ylabel(cb_ax, 'Z-score', 'FontSize', g.fontSize2);
cb_ax.FontSize = g.fontSize2;

%% Add Delta Peak FR bar charts in column 4 of the same nested layout
cluster_names = {'CS-sel', 'US-sel', 'Multi'};

% Storage for Kruskal-Wallis test data
% Structure: kw_data_storage{region, cluster, stimulus}
kw_data_storage = cell(2, 3, 3);  % 2 regions × 3 clusters × 3 stimuli (CS, US, Both)

for br = 1:2
    if isempty(results_all{br})
        continue;
    end

    res = results_all{br};

    % Create nested 3×1 tiledlayout for bars in column 4
    % Row 1 (br=1, LA): tile 4
    % Row 2 (br=2, Astria): tile 8
    t_bars = tiledlayout(t_heatmaps, 3, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    t_bars.Layout.Tile = (br-1)*4 + 4;  % Tile 4 for LA, Tile 8 for Astria

    % Plot each cluster in separate row
    for c = [1 2 3]  % CS-selective, US-selective, Multisensory
        ax_bar = nexttile(t_bars, c);

        clust_idx = find(res.Clusters_all == c);

        if ~isempty(clust_idx)
            % Calculate Delta Peak FR in monosyn window (12ms to g.monosyn_window)
            CS_metric = zeros(length(clust_idx), 1);
            US_metric = zeros(length(clust_idx), 1);
            Both_metric = zeros(length(clust_idx), 1);

            baseline_idx = 1:baseline_bins;
            response_start_bin = round((g.pre_time + 0.012) / g.bin_time);
            response_end_bin = round((g.pre_time + g.monosyn_window) / g.bin_time);

            for n = 1:length(clust_idx)
                idx_n = clust_idx(n);
                baseline_fr = mean(res.psth_CS_Hz_all(idx_n, baseline_idx));

                % CS peak firing rate change
                peak_fr_CS = max(res.psth_CS_Hz_all(idx_n, response_start_bin:response_end_bin));
                CS_metric(n) = peak_fr_CS - baseline_fr;

                % US peak firing rate change
                peak_fr_US = max(res.psth_US_Hz_all(idx_n, response_start_bin:response_end_bin));
                US_metric(n) = peak_fr_US - baseline_fr;

                % CS+US peak firing rate change
                peak_fr_Both = max(res.psth_Both_Hz_all(idx_n, response_start_bin:response_end_bin));
                Both_metric(n) = peak_fr_Both - baseline_fr;
            end

            % Store data for Kruskal-Wallis test (region, cluster, stimulus)
            kw_data_storage{br, c, 1} = CS_metric;
            kw_data_storage{br, c, 2} = US_metric;
            kw_data_storage{br, c, 3} = Both_metric;

            % Calculate means and SEMs
            means_data = mean([CS_metric, US_metric, Both_metric], 1, 'omitnan');
            sems_data = std([CS_metric, US_metric, Both_metric], 0, 1, 'omitnan') ./ sqrt(sum(~isnan(CS_metric)));

            % Plot bars
            hold on;
            bar_color = cluster_colors(c, :);  % Same color as pie charts
            bar([1 2 3], means_data, 0.4, 'FaceColor', bar_color, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([1 2 3], means_data, sems_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

            % First perform Kruskal-Wallis test
            all_values = [CS_metric; US_metric; Both_metric];
            group_labels = [ones(length(CS_metric), 1); 2*ones(length(US_metric), 1); 3*ones(length(Both_metric), 1)];
            [p_kw, ~, ~] = kruskalwallis(all_values, group_labels, 'off');

            % Only perform post-hoc tests if Kruskal-Wallis is significant
            if p_kw < 0.05
                [p_CS_US, ~] = signrank(CS_metric, US_metric);
                [p_CS_Both, ~] = signrank(CS_metric, Both_metric);
                [p_US_Both, ~] = signrank(US_metric, Both_metric);
            else
                % Set p-values to 1 (non-significant) if KW is not significant
                p_CS_US = 1;
                p_CS_Both = 1;
                p_US_Both = 1;
            end

            hold off;
        end

        % Formatting
        xlim([0.5 3.5]);
        xticks([1 2 3]);
        ylim([0 320]);
        yticks([0 160 320]);

        % Add significance markers after setting ylim
        if ~isempty(clust_idx)
            hold on;
            y_max_data = max(means_data + sems_data);
            curr_ylim = ylim;
            y_range = curr_ylim(2) - curr_ylim(1);

            % Level 1: Short comparisons (CS vs US, US vs CS+US)
            y_pos_level1 = y_max_data + 0.08 * y_range;

            % CS vs US
            if p_CS_US < 0.001
                sig_text = '***';
            elseif p_CS_US < 0.01
                sig_text = '**';
            elseif p_CS_US < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([1.1 1.9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                text(1.5, y_pos_level1 + 0.02 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            % US vs CS+US
            if p_US_Both < 0.001
                sig_text = '***';
            elseif p_US_Both < 0.01
                sig_text = '**';
            elseif p_US_Both < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([2.1 2.9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                text(2.5, y_pos_level1 + 0.02 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            % Level 2: Long comparison (CS vs CS+US)
            y_pos_level2 = y_pos_level1 + 25;  % Increased spacing between levels

            % CS vs CS+US
            if p_CS_Both < 0.001
                sig_text = '***';
            elseif p_CS_Both < 0.01
                sig_text = '**';
            elseif p_CS_Both < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([1 3], [y_pos_level2 y_pos_level2], 'k-', 'LineWidth', 1.5);
                text(2, y_pos_level2 + 0.02 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            hold off;
        end

        % Add title on top cluster of LA
        if br == 1 && c == 1
            title('ΔFR (Hz)', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end

        % X-tick labels only on bottom cluster of Astria
        if br == 2 && c == 3
            xticklabels({'CS', 'US', 'CS+US'});
        else
            set(gca, 'XTickLabel', []);
        end

        set(gca, 'FontSize', g.fontSize2);
        box off;
    end
end

%% Add right panel with pie charts, across-region comparison, and latency comparison
% Create main nested 3x2 tiledlayout spanning rows 1-4, columns 5-6
t_right = tiledlayout(t, 3, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
t_right.Layout.Tile = 5;  % Row 1, columns 5-6
t_right.Layout.TileSpan = [4 2];  % Span 4 rows, 2 columns

%% Row 1: Pie charts (1×2 nested tiledlayout in tiles 1-2)
% Create nested 1×2 tiledlayout for pie charts
t_pie = tiledlayout(t_right, 1, 2, 'TileSpacing', 'none', 'Padding', 'tight');
t_pie.Layout.Tile = 1;  % Start at tile 1
t_pie.Layout.TileSpan = [1 2];  % Span tiles 1-2

% LA pie chart
ax_pie_LA = nexttile(t_pie, 1);

if ~isempty(results_all{1})
    res = results_all{1};

    % Count neurons in each cluster (only clusters 1, 2, 3)
    n_CS_sel = sum(res.Clusters_all == 1);
    n_US_sel = sum(res.Clusters_all == 2);
    n_Multi = sum(res.Clusters_all == 3);
    total_n = n_CS_sel + n_US_sel + n_Multi;

    pie_data = [n_CS_sel, n_US_sel, n_Multi];
    % Calculate percentages
    pct_CS = (n_CS_sel / total_n) * 100;
    pct_US = (n_US_sel / total_n) * 100;
    pct_Multi = (n_Multi / total_n) * 100;

    pie_labels = {sprintf('%.0f%%', pct_CS), ...
                  sprintf('%.0f%%', pct_US), ...
                  sprintf('%.0f%%', pct_Multi)};

    p = pie(ax_pie_LA, pie_data, pie_labels);

    % Set colors and text properties
    for i = 1:2:length(p)  % Every other element is a patch
        patch_idx = (i+1)/2;
        set(p(i), 'FaceColor', cluster_colors(patch_idx, :));
        set(p(i), 'EdgeColor', 'k');
        set(p(i), 'LineWidth', 1);
    end

    % Set text properties and move inside pie - FontSize 14
    for i = 2:2:length(p)  % Text elements
        set(p(i), 'FontSize', 14);
        set(p(i), 'FontWeight', 'bold');
        set(p(i), 'Color', 'w');  % White font color
        % Move text closer to center (inside the pie)
        pos = get(p(i), 'Position');
        set(p(i), 'Position', pos * 0.3);  % Move 30% toward center (more inside)
    end

    title('LA', 'FontSize', 12, 'FontWeight', 'bold');
end

% Astria pie chart
ax_pie_Astria = nexttile(t_pie, 2);

if ~isempty(results_all{2})
    res = results_all{2};

    % Count neurons in each cluster (only clusters 1, 2, 3)
    n_CS_sel = sum(res.Clusters_all == 1);
    n_US_sel = sum(res.Clusters_all == 2);
    n_Multi = sum(res.Clusters_all == 3);
    total_n = n_CS_sel + n_US_sel + n_Multi;

    pie_data = [n_CS_sel, n_US_sel, n_Multi];
    % Calculate percentages
    pct_CS = (n_CS_sel / total_n) * 100;
    pct_US = (n_US_sel / total_n) * 100;
    pct_Multi = (n_Multi / total_n) * 100;

    pie_labels = {sprintf('%.0f%%', pct_CS), ...
                  sprintf('%.0f%%', pct_US), ...
                  sprintf('%.0f%%', pct_Multi)};

    p = pie(ax_pie_Astria, pie_data, pie_labels);

    % Set colors and text properties
    for i = 1:2:length(p)  % Every other element is a patch
        patch_idx = (i+1)/2;
        set(p(i), 'FaceColor', cluster_colors(patch_idx, :));
        set(p(i), 'EdgeColor', 'k');
        set(p(i), 'LineWidth', 1);
    end

    % Set text properties and move inside pie - FontSize 14
    for i = 2:2:length(p)  % Text elements
        set(p(i), 'FontSize', 14);
        set(p(i), 'FontWeight', 'bold');
        set(p(i), 'Color', 'w');  % White font color
        % Move text closer to center (inside the pie)
        pos = get(p(i), 'Position');
        set(p(i), 'Position', pos * 0.33);  % Move 30% toward center (more inside)
    end

    title('AStria', 'FontSize', 12, 'FontWeight', 'bold');
end

% Add legend between the two pie charts
% Create a dummy invisible axes positioned between pie charts
ax_legend = axes('Position', [0.65 0.85 0.05 0.05], 'Visible', 'off');
hold(ax_legend, 'on');
% Create invisible plot objects with square markers for legend
h1 = plot(ax_legend, NaN, NaN, 's', 'MarkerSize', 12, 'MarkerFaceColor', cluster_colors(1, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
h2 = plot(ax_legend, NaN, NaN, 's', 'MarkerSize', 12, 'MarkerFaceColor', cluster_colors(2, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
h3 = plot(ax_legend, NaN, NaN, 's', 'MarkerSize', 12, 'MarkerFaceColor', cluster_colors(3, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
lgd = legend(ax_legend, [h1, h2, h3], {'CS-sel', 'US-sel', 'Multi'}, ...
             'Location', 'east', 'FontSize', 10, 'Box', 'off');
lgd.ItemTokenSize = [30, 30];  % Increase marker size in legend

% Add title to pie chart section - FontSize 12
title(t_pie, 'Response categories', 'FontSize', 12, 'FontWeight', 'bold');

%% Row 2: Across-region comparison bar plots (3 nested tiles spanning both columns)
% Create nested 1×3 tiledlayout for comparison bars in row 2
t_comp = tiledlayout(t_right, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
t_comp.Layout.Tile = 3;  % Row 2, spanning tiles 3-4
t_comp.Layout.TileSpan = [1 2];  % Span 1 row, 2 columns

stim_names = {'CS', 'US', 'CS+US'};

for stim = 1:3  % CS, US, CS+US
    ax_comp = nexttile(t_comp, stim);

    % Collect data from LA and Astria
    data_regions = [];
    region_labels = {};

    for br = 1:2
        if isempty(results_all{br})
            continue;
        end

        res = results_all{br};

        % Select appropriate PSTH
        if stim == 1  % CS
            psth_Hz = res.psth_CS_Hz_all;
        elseif stim == 2  % US
            psth_Hz = res.psth_US_Hz_all;
        else  % CS+US
            psth_Hz = res.psth_Both_Hz_all;
        end

        % Find responsive neurons (clusters 1, 2, 3)
        responsive_idx = find(ismember(res.Clusters_all, [1 2 3]));

        if ~isempty(responsive_idx)
            % Calculate Delta Peak FR in monosyn window (12ms to g.monosyn_window)
            delta_peak_fr = zeros(length(responsive_idx), 1);
            baseline_idx = 1:baseline_bins;
            response_start_bin = round((g.pre_time + 0.012) / g.bin_time);
            response_end_bin = round((g.pre_time + g.monosyn_window) / g.bin_time);

            for n = 1:length(responsive_idx)
                idx_n = responsive_idx(n);

                % Calculate baseline firing rate
                baseline_fr = mean(psth_Hz(idx_n, baseline_idx));

                % Calculate peak firing rate in monosyn window
                peak_fr = max(psth_Hz(idx_n, response_start_bin:response_end_bin));

                % Delta Peak FR
                delta_peak_fr(n) = peak_fr - baseline_fr;
            end

            data_regions = [data_regions; delta_peak_fr];

            if br == 1
                region_labels = [region_labels; repmat({'LA'}, length(responsive_idx), 1)];
            else
                region_labels = [region_labels; repmat({'AStria'}, length(responsive_idx), 1)];
            end
        end
    end

    % Plot bars
    if ~isempty(data_regions)
        % Separate LA and Astria data
        LA_data = data_regions(strcmp(region_labels, 'LA'));
        Astria_data = data_regions(strcmp(region_labels, 'AStria'));

        means_data = [mean(LA_data, 'omitnan'), mean(Astria_data, 'omitnan')];
        sems_data = [std(LA_data, 0, 'omitnan') / sqrt(length(LA_data)), ...
                     std(Astria_data, 0, 'omitnan') / sqrt(length(Astria_data))];

        hold on;

        % Bar colors
        bar_colors = [0.7 0.2 0.2; 0.2 0.4 0.7];  % LA dark red, Astria blue

        b = bar([1 2], means_data, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1);
        b.CData = bar_colors;

        % Plot individual data points as grey circles with jitter
        n_LA = length(LA_data);
        n_Astria = length(Astria_data);
        x_jitter_LA = 1 + (rand(n_LA, 1) - 0.5) * 0.15;
        x_jitter_Astria = 2 + (rand(n_Astria, 1) - 0.5) * 0.15;
        scatter(x_jitter_LA, LA_data, 16, [0.5 0.5 0.5], 'LineWidth', 0.5);
        scatter(x_jitter_Astria, Astria_data, 16, [0.5 0.5 0.5], 'LineWidth', 0.5);

        errorbar([1 2], means_data, sems_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

        % Statistical comparison - Wilcoxon rank-sum test (unpaired)
        [p_val, ~] = ranksum(LA_data, Astria_data);

        % Use 95th percentile to avoid compression from extreme outliers
        all_data_values = [LA_data; Astria_data];
        data_95th = prctile(all_data_values, 95);
        y_max = max(max(means_data + sems_data), data_95th);
        y_range = y_max * 1.5;
        y_pos = y_max + 0.1 * y_range;

        plot([1 2], [y_pos y_pos], 'k-', 'LineWidth', 1.5);

        if p_val < 0.001
            sig_text = '***';
        elseif p_val < 0.01
            sig_text = '**';
        elseif p_val < 0.05
            sig_text = '*';
        else
            sig_text = 'n.s.';
        end

        text(1.5, y_pos, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', g.fontSize2);

        hold off;

        % Y-axis scaled to 95th percentile to avoid compression from outliers
        ylim([0 y_range]);
    else
        % No data - use default limits
        ylim([0 250]);
    end

    % Formatting
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'LA', 'AStria'});

    % Y-label and Y-tick labels only on first panel
    if stim == 1
        ylabel('ΔFR (Hz)', 'FontSize', g.fontSize2);
    else
        set(gca, 'YTickLabel', []);
    end

    title(stim_names{stim}, 'FontSize', g.fontSize1, 'FontWeight', 'bold');
    set(gca, 'FontSize', g.fontSize2);
    box off;
end

% Add title to comparison section - FontSize 12
title(t_comp, 'Across region comparison', 'FontSize', 12, 'FontWeight', 'bold');

%% Row 3: Latency comparison boxplots (3 nested tiles spanning both columns)
% Create nested 1×3 tiledlayout for latency comparisons in row 3
t_latency = tiledlayout(t_right, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
t_latency.Layout.Tile = 5;  % Row 3, spanning tiles 5-6
t_latency.Layout.TileSpan = [1 2];  % Span 1 row, 2 columns

stim_names_lat = {'CS', 'US', 'CS+US'};

for stim = 1:3  % CS, US, CS+US
    ax_lat = nexttile(t_latency, stim);

    % Collect latency data for this stimulus
    data_regions = [];
    region_labels = {};

    for br = 1:2
        if isempty(results_all{br})
            continue;
        end

        res = results_all{br};

        % Select appropriate latencies based on stimulus and cluster
        if stim == 1  % CS - from CS-selective and Multisensory
            responsive_idx = find(ismember(res.Clusters_all, [1 3]));
            onset_lat = res.CS_onset_lat_all(responsive_idx);
        elseif stim == 2  % US - from US-selective and Multisensory
            responsive_idx = find(ismember(res.Clusters_all, [2 3]));
            onset_lat = res.US_onset_lat_all(responsive_idx);
        else  % CS+US - from all responsive neurons
            responsive_idx = find(ismember(res.Clusters_all, [1 2 3]));
            onset_lat = res.Both_onset_lat_all(responsive_idx);
        end

        % Remove NaN values and convert to ms
        onset_lat = onset_lat(~isnan(onset_lat)) * 1000;

        if ~isempty(onset_lat)
            data_regions = [data_regions; onset_lat];

            if br == 1
                region_labels = [region_labels; repmat({'LA'}, length(onset_lat), 1)];
            else
                region_labels = [region_labels; repmat({'AStria'}, length(onset_lat), 1)];
            end
        end
    end

    % Plot boxplots and scatter
    if ~isempty(data_regions)
        % Separate LA and Astria data
        LA_data = data_regions(strcmp(region_labels, 'LA'));
        Astria_data = data_regions(strcmp(region_labels, 'AStria'));

        hold on;

        % Plot boxplots
        boxplot(ax_lat, [LA_data; Astria_data], ...
                [ones(length(LA_data),1)*1; ones(length(Astria_data),1)*2], ...
                'Positions', [1 2], 'Widths', 0.6, 'Colors', 'k');

        % Overlay scatter points with jitter
        jitter_amount = 0.15;
        region_colors_lat = {[0.8 0.2 0.2], [0.2 0.4 0.8]};  % Red for LA, Blue for Astria

        if ~isempty(LA_data)
            x_jitter_LA = 1 + (rand(length(LA_data), 1) - 0.5) * jitter_amount;
            scatter(x_jitter_LA, LA_data, 20, region_colors_lat{1}, 'filled', 'MarkerFaceAlpha', 0.4);
        end

        if ~isempty(Astria_data)
            x_jitter_Astria = 2 + (rand(length(Astria_data), 1) - 0.5) * jitter_amount;
            scatter(x_jitter_Astria, Astria_data, 20, region_colors_lat{2}, 'filled', 'MarkerFaceAlpha', 0.4);
        end

        % Statistical comparison - Wilcoxon rank-sum test (unpaired)
        if ~isempty(LA_data) && ~isempty(Astria_data)
            [p_val, ~] = ranksum(LA_data, Astria_data);

            % Add significance marker
            y_max = max([LA_data; Astria_data]);
            y_pos = y_max * 1.15;

            plot([1 2], [y_pos y_pos], 'k-', 'LineWidth', 1.5);

            if p_val < 0.001
                sig_text = '***';
            elseif p_val < 0.01
                sig_text = '**';
            elseif p_val < 0.05
                sig_text = '*';
            else
                sig_text = 'n.s.';
            end

            text(1.5, y_pos, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', g.fontSize2);
        end

        hold off;
    end

    % Formatting
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'LA', 'AStria'});
    ylim([0 50]);
    yticks([0 25 50]);

    % Y-label and Y-tick labels only on first panel
    if stim == 1
        ylabel('Onset Latency (ms)', 'FontSize', g.fontSize2);
    else
        set(gca, 'YTickLabel', []);
    end

    title(stim_names_lat{stim}, 'FontSize', g.fontSize1, 'FontWeight', 'bold');
    set(gca, 'FontSize', g.fontSize2);
    box off;
end

% Add title to latency section - FontSize 12
title(t_latency, 'Onset latency comparison', 'FontSize', 12, 'FontWeight', 'bold');

fprintf('\nMain figure complete.\n');

%% Supplementary figure - selective vs multisensory latency comparison + chi-square test
fprintf('\nGenerating supplementary figure...\n');

fig_supp = figure('Units', 'pixels', 'Position', [200, 200, 1000, 900], 'Visible', 'on');
t_supp = tiledlayout(fig_supp, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add panel labels for supplementary figure
% A: Onset latencies (row 1)
annotation(fig_supp, 'textbox', [0.01 0.94 0.05 0.05], 'String', 'A', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Row 1: Latency comparison (original boxplot content)
ax_latency = nexttile(t_supp, 1);

% Collect latency data for each group
% LA CS-selective and Multisensory (CS latency)
LA_CS_sel_lat = results_all{1}.CS_onset_lat_all(results_all{1}.Clusters_all == 1) * 1000;  % Convert to ms
LA_Multi_CS_lat = results_all{1}.CS_onset_lat_all(results_all{1}.Clusters_all == 3) * 1000;

% Astria CS-selective and Multisensory (CS latency)
Astria_CS_sel_lat = results_all{2}.CS_onset_lat_all(results_all{2}.Clusters_all == 1) * 1000;
Astria_Multi_CS_lat = results_all{2}.CS_onset_lat_all(results_all{2}.Clusters_all == 3) * 1000;

% LA US-selective and Multisensory (US latency)
LA_US_sel_lat = results_all{1}.US_onset_lat_all(results_all{1}.Clusters_all == 2) * 1000;
LA_Multi_US_lat = results_all{1}.US_onset_lat_all(results_all{1}.Clusters_all == 3) * 1000;

% Astria US-selective and Multisensory (US latency)
Astria_US_sel_lat = results_all{2}.US_onset_lat_all(results_all{2}.Clusters_all == 2) * 1000;
Astria_Multi_US_lat = results_all{2}.US_onset_lat_all(results_all{2}.Clusters_all == 3) * 1000;

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
xticklabels({'LA CS', 'AStria CS', 'LA US', 'AStria US'});
ylabel('Onset Latency (ms)', 'FontSize', 10);
title('Onset Latencies: Selective vs Multisensory', 'FontSize', 10, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);
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
            'VerticalAlignment', 'bottom', 'FontSize', 10);
    end
end
hold off;

% Add legend
legend_ax = axes('Position', [0.75 0.8 0.15 0.1], 'Visible', 'off');
hold(legend_ax, 'on');
h1 = scatter(legend_ax, NaN, NaN, 50, [0.5 0.5 0.5], 'filled');
h2 = scatter(legend_ax, NaN, NaN, 50, [0.6 0.2 0.6], 'filled');
legend(legend_ax, [h1, h2], {'Selective', 'Multisensory'}, 'Location', 'best', 'FontSize', 10);
hold(legend_ax, 'off');

%% Row 2: Chi-square test analysis (3 panels)
fprintf('\n=== Chi-square Test: LA vs Astria (Monosynaptic) ===\n');

% Build contingency table: 2 regions × 3 clusters (CS-sel, US-sel, Multi)
contingency_table = zeros(2, 3);
for br = 1:2
    if ~isempty(results_all{br})
        for c = 1:3
            contingency_table(br, c) = sum(results_all{br}.Clusters_all == c);
        end
    end
end

fprintf('Contingency Table (regions × clusters):\n');
fprintf('%12s%12s%12s%12s%12s\n', 'Region', 'CS-sel', 'US-sel', 'Multi', 'Total');
for br = 1:2
    if strcmp(brain_regions{br}, 'Astria')
        fprintf('%12s', 'AStria');
    else
        fprintf('%12s', brain_regions{br});
    end
    for c = 1:3
        fprintf('%12d', contingency_table(br, c));
    end
    fprintf('%12d\n', sum(contingency_table(br, :)));
end
fprintf('%12s', 'Total');
for c = 1:3
    fprintf('%12d', sum(contingency_table(:, c)));
end
fprintf('%12d\n', sum(contingency_table(:)));

% Calculate chi-square statistic
[chi2_obs, p_parametric] = calculate_chi_square(contingency_table);
fprintf('\nObserved chi-square statistic: %.4f\n', chi2_obs);
fprintf('Parametric p-value: %.4f\n', p_parametric);

% Cramér's V for effect size
n_total = sum(contingency_table(:));
df_cramer = min(size(contingency_table, 1) - 1, size(contingency_table, 2) - 1);
cramers_v = sqrt(chi2_obs / (n_total * df_cramer));
fprintf('Cramér''s V (effect size): %.4f\n', cramers_v);

% Check expected counts
row_totals = sum(contingency_table, 2);
col_totals = sum(contingency_table, 1);
expected = (row_totals * col_totals) / n_total;
min_expected = min(expected(:));
if min_expected < 5
    fprintf('WARNING: Minimum expected count = %.2f (< 5). Chi-square approximation may be unreliable.\n', min_expected);
else
    fprintf('Minimum expected count = %.2f (>= 5). Chi-square assumptions satisfied.\n', min_expected);
end

% Permutation test
n_permutations = 10000;
chi2_perm = zeros(n_permutations, 1);

% Pool all cluster assignments
all_clusters = [];
region_indices = [];
for br = 1:2
    if ~isempty(results_all{br})
        n_neurons = results_all{br}.n_neurons;
        all_clusters = [all_clusters; results_all{br}.Clusters_all];
        region_indices = [region_indices; br * ones(n_neurons, 1)];
    end
end

fprintf('Running %d permutations...\n', n_permutations);
for perm = 1:n_permutations
    % Shuffle cluster assignments
    shuffled_clusters = all_clusters(randperm(length(all_clusters)));

    % Build permuted contingency table
    perm_table = zeros(2, 3);
    for br = 1:2
        br_mask = region_indices == br;
        br_clusters = shuffled_clusters(br_mask);
        for c = 1:3
            perm_table(br, c) = sum(br_clusters == c);
        end
    end

    % Calculate chi-square for permuted data
    chi2_perm(perm) = calculate_chi_square(perm_table);

    if mod(perm, 1000) == 0
        fprintf('  Completed %d/%d permutations\n', perm, n_permutations);
    end
end

% Calculate permutation p-value
p_perm = sum(chi2_perm >= chi2_obs) / n_permutations;
fprintf('\nPermutation p-value: %.4f\n', p_perm);
if p_perm < 0.05
    fprintf('Significance (p < 0.05): YES\n');
else
    fprintf('Significance (p < 0.05): NO\n');
end

% Create nested 1×3 tiledlayout for chi-square panels
t_chi = tiledlayout(t_supp, 1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
t_chi.Layout.Tile = 2;

% B: Permutation distribution (row 2, col 1)
annotation(fig_supp, 'textbox', [0.01 0.63 0.05 0.05], 'String', 'B', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% C: Contingency table (row 2, col 2)
annotation(fig_supp, 'textbox', [0.35 0.63 0.05 0.05], 'String', 'C', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% D: Statistics (row 2, col 3)
annotation(fig_supp, 'textbox', [0.68 0.63 0.05 0.05], 'String', 'D', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% E: Kruskal-Wallis tests (row 3)
annotation(fig_supp, 'textbox', [0.01 0.30 0.05 0.05], 'String', 'E', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Panel 1: Permutation distribution
ax_perm = nexttile(t_chi, 1);
histogram(ax_perm, chi2_perm, 50, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
hold on;
xline(chi2_obs, 'r-', 'LineWidth', 2);
hold off;
xlabel('Chi-square statistic', 'FontSize', 10);
ylabel('Frequency', 'FontSize', 10);
title('Permutation Distribution', 'FontSize', 10, 'FontWeight', 'bold');
legend({'Null', 'Observed'}, 'Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 10);
box off;

% Panel 2: Contingency table
ax_table = nexttile(t_chi, 2);
axis(ax_table, 'off');
% Display contingency table as text
y_start = 0.7;
y_step = 0.2;
text(0.1, y_start, 'Region', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.3, y_start, 'CS-sel', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.5, y_start, 'US-sel', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.7, y_start, 'Multi', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
for br = 1:2
    if strcmp(brain_regions{br}, 'Astria')
        region_name = 'AStria';
    else
        region_name = brain_regions{br};
    end
    text(0.1, y_start - br*y_step, region_name, 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(0.3, y_start - br*y_step, sprintf('%d', contingency_table(br, 1)), 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(0.5, y_start - br*y_step, sprintf('%d', contingency_table(br, 2)), 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(0.7, y_start - br*y_step, sprintf('%d', contingency_table(br, 3)), 'FontSize', 10, 'HorizontalAlignment', 'center');
end
title(ax_table, sprintf('Contingency Table (χ²=%.2f)', chi2_obs), 'FontSize', 10, 'FontWeight', 'bold');

% Panel 3: Statistical summary
ax_stats = nexttile(t_chi, 3);
axis(ax_stats, 'off');
% Compact display of statistics
y_pos = 0.7;
y_step = 0.15;
text(0.5, y_pos, sprintf('χ²: %.2f', chi2_obs), 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.5, y_pos - y_step, sprintf('p: %.4f', p_perm), 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.5, y_pos - 2*y_step, sprintf('V: %.3f', cramers_v), 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Add significance symbol
if p_perm < 0.001
    sig_symbol = '***';
    sig_color = [0.8 0.2 0.2];
elseif p_perm < 0.01
    sig_symbol = '**';
    sig_color = [0.8 0.4 0.2];
elseif p_perm < 0.05
    sig_symbol = '*';
    sig_color = [0.8 0.6 0.2];
else
    sig_symbol = 'n.s.';
    sig_color = [0.5 0.5 0.5];
end
text(0.5, y_pos - 3*y_step, sig_symbol, 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Color', sig_color);
title(ax_stats, 'Statistics', 'FontSize', 10, 'FontWeight', 'bold');

%% Row 3: Kruskal-Wallis test visualization
% Perform Kruskal-Wallis tests for each region × cluster combination
kw_results = struct();

fprintf('\n=== Kruskal-Wallis Tests for Bar Charts (CS vs US vs CS+US within each cluster) ===\n');

cluster_names_kw = {'CS-sel', 'US-sel', 'Multi'};

% Create nested 2×3 tiledlayout for KW test panels (2 regions × 3 clusters)
t_kw = tiledlayout(t_supp, 2, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
t_kw.Layout.Tile = 3;

for br = 1:2
    if isempty(results_all{br})
        continue;
    end

    fprintf('\n--- %s ---\n', brain_regions{br});

    for c = 1:3  % CS-selective, US-selective, Multisensory
        % Check if we have data for this cluster
        if isempty(kw_data_storage{br, c, 1})
            fprintf('  %s: No data\n', cluster_names_kw{c});
            continue;
        end

        % Combine data for KW test: CS, US, Both stimuli
        all_values = [kw_data_storage{br, c, 1}; kw_data_storage{br, c, 2}; kw_data_storage{br, c, 3}];
        group_labels = [ones(length(kw_data_storage{br, c, 1}), 1); ...
                        2*ones(length(kw_data_storage{br, c, 2}), 1); ...
                        3*ones(length(kw_data_storage{br, c, 3}), 1)];

        % Kruskal-Wallis test
        [p_kw, ~, ~] = kruskalwallis(all_values, group_labels, 'off');

        % Only perform post-hoc tests if KW is significant
        if p_kw < 0.05
            % Post-hoc pairwise Wilcoxon signed-rank tests
            p_values = zeros(3, 1);
            p_values(1) = signrank(kw_data_storage{br, c, 1}, kw_data_storage{br, c, 2});  % CS vs US
            p_values(2) = signrank(kw_data_storage{br, c, 1}, kw_data_storage{br, c, 3});  % CS vs Both
            p_values(3) = signrank(kw_data_storage{br, c, 2}, kw_data_storage{br, c, 3});  % US vs Both
        else
            % Set p-values to 1 (non-significant) if KW is not significant
            p_values = ones(3, 1);
        end

        % Store results
        kw_results(br, c).p_kw = p_kw;
        kw_results(br, c).p_values = p_values;

        fprintf('  %s: KW p=%.4f, CS-US p=%.4f, CS-Multi p=%.4f, US-Multi p=%.4f\n', ...
            cluster_names_kw{c}, p_kw, p_values(1), p_values(2), p_values(3));

        % Calculate tile index for this region × cluster combination
        % Row 1 (br=1, LA): tiles 1, 2, 3
        % Row 2 (br=2, AStria): tiles 4, 5, 6
        tile_idx = (br-1)*3 + c;

        % Create panel with KW p-value and post-hoc results combined
        ax = nexttile(t_kw, tile_idx);
        axis off;

        % KW p-value at top
        text(0.5, 0.85, sprintf('KW p = %.4f', p_kw), 'FontSize', 9, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        if p_kw < 0.001
            sig_text = '***';
            text_color = [0.8 0.2 0.2];
        elseif p_kw < 0.01
            sig_text = '**';
            text_color = [0.8 0.2 0.2];
        elseif p_kw < 0.05
            sig_text = '*';
            text_color = [0.8 0.2 0.2];
        else
            sig_text = 'n.s.';
            text_color = [0.4 0.4 0.4];
        end
        text(0.5, 0.65, sig_text, 'FontSize', 11, 'FontWeight', 'bold', 'Color', text_color, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

        % Post-hoc p-values below (only show if KW is significant)
        if p_kw < 0.05
            y_pos = 0.45;
            y_step = 0.15;
            text(0.5, y_pos, 'Post-hoc:', 'FontSize', 8, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            y_pos = y_pos - y_step;
            text(0.5, y_pos, sprintf('CS-US: %.3f', p_values(1)), 'FontSize', 7, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            y_pos = y_pos - y_step;
            text(0.5, y_pos, sprintf('CS-Multi: %.3f', p_values(2)), 'FontSize', 7, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            y_pos = y_pos - y_step;
            text(0.5, y_pos, sprintf('US-Multi: %.3f', p_values(3)), 'FontSize', 7, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end

        % Title showing region and cluster
        if strcmp(brain_regions{br}, 'Astria')
            title_text = sprintf('AStria %s', cluster_names_kw{c});
        else
            title_text = sprintf('%s %s', brain_regions{br}, cluster_names_kw{c});
        end
        title(title_text, 'FontSize', 10, 'FontWeight', 'bold');
    end
end

fprintf('\nSupplementary figure complete.\n');

%% Export statistics
export_figure4_stats(results_all, contingency_table, chi2_obs, p_perm, cramers_v, ...
    kw_results, kw_data_storage, brain_regions, cluster_names, g);

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

function [chi2_stat, p_value] = calculate_chi_square(contingency_table)
    % Calculate chi-square statistic from contingency table
    % Returns chi2 statistic and parametric p-value

    % Calculate expected frequencies
    row_totals = sum(contingency_table, 2);
    col_totals = sum(contingency_table, 1);
    n_total = sum(contingency_table(:));

    expected = (row_totals * col_totals) / n_total;

    % Calculate chi-square statistic
    chi2_stat = sum(((contingency_table(:) - expected(:)).^2) ./ expected(:));

    % Degrees of freedom
    df = (size(contingency_table, 1) - 1) * (size(contingency_table, 2) - 1);

    % Parametric p-value (chi-square distribution)
    p_value = 1 - chi2cdf(chi2_stat, df);
end

function export_figure4_stats(results_all, contingency_table, chi2_obs, p_perm, cramers_v, ...
    kw_results, kw_data_storage, brain_regions, cluster_names, g)

    fid = fopen('figure_4_stats.txt', 'w');

    fprintf(fid, '========================================\n');
    fprintf(fid, 'FIGURE 4 STATISTICS (Monosynaptic Responses)\n');
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, '========================================\n\n');

    fprintf(fid, 'This file contains statistics for:\n');
    fprintf(fid, '  - MAIN FIGURE (Panels A-G): Heatmaps, delta FR bars, pie charts, region comparison, latency comparison\n');
    fprintf(fid, '  - SUPPLEMENTARY FIGURE (Panels A-E): Latency boxplots, chi-square test, Kruskal-Wallis tests\n\n');

    fprintf(fid, 'Analysis window: 0-%d ms (monosynaptic)\n', g.monosyn_window*1000);
    fprintf(fid, 'Rank score alpha: %.1f\n\n', g.alpha);

    %% Sample sizes
    fprintf(fid, '### SAMPLE SIZES BY REGION ###\n\n');

    all_animals = {};
    for br = 1:2
        if ~isempty(results_all{br})
            animals_br = results_all{br}.animals;
            if ~iscell(animals_br)
                animals_br = {animals_br};
            end
            all_animals = [all_animals; animals_br(:)];
        end
    end
    unique_animals = unique(all_animals);

    fprintf(fid, 'Number of animals (N): %d\n', length(unique_animals));
    fprintf(fid, 'Animal IDs: %s\n\n', strjoin(unique_animals, ', '));

    total_neurons = 0;
    for br = 1:2
        if ~isempty(results_all{br})
            n_neurons = results_all{br}.n_neurons;
            total_neurons = total_neurons + n_neurons;
            fprintf(fid, '%s: n = %d neurons\n', brain_regions{br}, n_neurons);
        end
    end
    fprintf(fid, 'Total neurons (LA + AStria): %d\n\n', total_neurons);

    %% Cluster distributions
    fprintf(fid, '========================================\n');
    fprintf(fid, 'MAIN FIGURE - PANEL E: PIE CHARTS\n');
    fprintf(fid, 'SUPPLEMENTARY FIGURE - PANELS B, C, D: CHI-SQUARE TEST\n');
    fprintf(fid, '========================================\n\n');

    fprintf(fid, '### CLUSTER DISTRIBUTIONS BY REGION ###\n\n');

    fprintf(fid, 'Contingency Table (neurons per cluster):\n');
    fprintf(fid, '%-10s', 'Region');
    for c = 1:3
        fprintf(fid, '%-12s', cluster_names{c});
    end
    fprintf(fid, '%-10s\n', 'Total');

    for br = 1:2
        if strcmp(brain_regions{br}, 'Astria')
            region_name = 'AStria';
        else
            region_name = brain_regions{br};
        end
        fprintf(fid, '%-10s', region_name);
        for c = 1:3
            fprintf(fid, '%-12d', contingency_table(br, c));
        end
        fprintf(fid, '%-10d\n', sum(contingency_table(br, :)));
    end
    fprintf(fid, '\n');

    fprintf(fid, 'Cluster Proportions by Region:\n');
    for br = 1:2
        if strcmp(brain_regions{br}, 'Astria')
            region_name = 'AStria';
        else
            region_name = brain_regions{br};
        end
        fprintf(fid, '%s:\n', region_name);
        total_region = sum(contingency_table(br, :));
        if total_region > 0
            for c = 1:3
                n_cluster = contingency_table(br, c);
                pct = 100 * n_cluster / total_region;
                fprintf(fid, '  %s: %d (%.1f%%)\n', cluster_names{c}, n_cluster, pct);
            end
        end
        fprintf(fid, '\n');
    end

    %% Chi-square test
    fprintf(fid, '\n### CHI-SQUARE TEST (LA vs AStria) ###\n\n');

    fprintf(fid, 'Overall Chi-square Test:\n');
    fprintf(fid, '  Chi-square statistic: %.4f\n', chi2_obs);
    fprintf(fid, '  Permutation p-value (10,000 permutations): %.4f\n', p_perm);
    fprintf(fid, '  Cramér''s V (effect size): %.4f\n', cramers_v);
    fprintf(fid, '  Significant (p < 0.05): %s\n\n', char(string(p_perm < 0.05)));

    %% Kruskal-Wallis tests
    fprintf(fid, '========================================\n');
    fprintf(fid, 'MAIN FIGURE - PANELS B, D: DELTA PEAK FR BARS\n');
    fprintf(fid, 'SUPPLEMENTARY FIGURE - PANEL E: KRUSKAL-WALLIS TESTS\n');
    fprintf(fid, '========================================\n\n');

    fprintf(fid, '### KRUSKAL-WALLIS TESTS (CS vs US vs CS+US WITHIN EACH CLUSTER) ###\n\n');

    cluster_names_kw = {'CS-sel', 'US-sel', 'Multi'};

    for br = 1:2
        if strcmp(brain_regions{br}, 'Astria')
            region_name = 'AStria';
        else
            region_name = brain_regions{br};
        end

        fprintf(fid, '--- %s ---\n', region_name);

        for c = 1:3
            fprintf(fid, '\n%s:\n', cluster_names_kw{c});
            fprintf(fid, '  Kruskal-Wallis p-value: %.4f\n', kw_results(br, c).p_kw);
            fprintf(fid, '  Significant (p < 0.05): %s\n', char(string(kw_results(br, c).p_kw < 0.05)));

            if kw_results(br, c).p_kw < 0.05
                fprintf(fid, '  Post-hoc comparisons (Wilcoxon signed-rank):\n');
                fprintf(fid, '    CS vs US: p = %.4f\n', kw_results(br, c).p_values(1));
                fprintf(fid, '    CS vs CS+US: p = %.4f\n', kw_results(br, c).p_values(2));
                fprintf(fid, '    US vs CS+US: p = %.4f\n', kw_results(br, c).p_values(3));
            end
        end
        fprintf(fid, '\n');
    end

    %% Descriptive statistics for Delta Peak FR bars
    fprintf(fid, '### DESCRIPTIVE STATISTICS FOR DELTA PEAK FR BARS (PANELS B, D) ###\n\n');

    cluster_names_stats = {'CS-selective', 'US-selective', 'Multisensory'};
    stim_names = {'CS', 'US', 'CS+US'};

    for br = 1:2
        if strcmp(brain_regions{br}, 'Astria')
            region_name = 'AStria';
        else
            region_name = brain_regions{br};
        end

        fprintf(fid, '--- %s ---\n', region_name);

        for c = 1:3
            % Use the exact data stored during plotting
            if isempty(kw_data_storage{br, c, 1})
                continue;
            end

            n_cluster = length(kw_data_storage{br, c, 1});
            fprintf(fid, '\n%s (n = %d neurons):\n', cluster_names_stats{c}, n_cluster);

            for stim = 1:3
                data = kw_data_storage{br, c, stim};
                fprintf(fid, '  %s: %.2f ± %.2f Hz (mean ± SEM), median = %.2f Hz, SD = %.2f Hz\n', ...
                    stim_names{stim}, mean(data), std(data)/sqrt(length(data)), ...
                    median(data), std(data));
            end
        end
        fprintf(fid, '\n');
    end

    %% Cluster-specific latencies
    fprintf(fid, '========================================\n');
    fprintf(fid, 'MAIN FIGURE - PANEL G: ONSET LATENCY COMPARISON\n');
    fprintf(fid, 'SUPPLEMENTARY FIGURE - PANEL A: LATENCY BOXPLOTS\n');
    fprintf(fid, '========================================\n\n');

    fprintf(fid, '### CLUSTER-SPECIFIC MONOSYNAPTIC LATENCIES ###\n\n');

    for br = 1:2
        if isempty(results_all{br})
            continue;
        end

        if strcmp(brain_regions{br}, 'Astria')
            region_name = 'AStria';
        else
            region_name = brain_regions{br};
        end

        fprintf(fid, '--- %s ---\n', region_name);
        res = results_all{br};

        for c = 1:3
            clust_idx = res.Clusters_all == c;
            n_cluster = sum(clust_idx);

            if n_cluster == 0
                continue;
            end

            fprintf(fid, '\n%s (n = %d):\n', cluster_names{c}, n_cluster);

            % Monosynaptic onset latencies
            if c == 1 || c == 3  % CS-sel or Multi
                cs_onsets = res.CS_onset_lat_all(clust_idx);
                cs_onsets = cs_onsets(~isnan(cs_onsets));
                if ~isempty(cs_onsets)
                    fprintf(fid, '  CS onset (monosynaptic): %.3f ± %.3f ms (mean ± SD), median = %.3f ms, n = %d\n', ...
                        mean(cs_onsets)*1000, std(cs_onsets)*1000, median(cs_onsets)*1000, length(cs_onsets));
                end
            end

            if c == 2 || c == 3  % US-sel or Multi
                us_onsets = res.US_onset_lat_all(clust_idx);
                us_onsets = us_onsets(~isnan(us_onsets));
                if ~isempty(us_onsets)
                    fprintf(fid, '  US onset (monosynaptic): %.3f ± %.3f ms (mean ± SD), median = %.3f ms, n = %d\n', ...
                        mean(us_onsets)*1000, std(us_onsets)*1000, median(us_onsets)*1000, length(us_onsets));
                end
            end
        end
        fprintf(fid, '\n');
    end

    %% Across-region comparison (Panel F)
    fprintf(fid, '========================================\n');
    fprintf(fid, 'MAIN FIGURE - PANEL F: ACROSS-REGION COMPARISON\n');
    fprintf(fid, '========================================\n\n');

    fprintf(fid, '### ACROSS-REGION COMPARISON (ALL RESPONSIVE NEURONS) ###\n\n');

    baseline_bins = round(g.pre_time / g.bin_time);
    stim_names_full = {'CS', 'US', 'CS+US'};

    for stim = 1:3
        fprintf(fid, '--- %s ---\n', stim_names_full{stim});

        for br = 1:2
            if isempty(results_all{br})
                continue;
            end

            res = results_all{br};

            % Select appropriate PSTH
            if stim == 1  % CS
                psth_Hz = res.psth_CS_Hz_all;
            elseif stim == 2  % US
                psth_Hz = res.psth_US_Hz_all;
            else  % CS+US
                psth_Hz = res.psth_Both_Hz_all;
            end

            % Find responsive neurons (clusters 1, 2, 3)
            responsive_idx = find(ismember(res.Clusters_all, [1 2 3]));

            if ~isempty(responsive_idx)
                % Calculate Delta Peak FR in monosyn window (12ms to g.monosyn_window)
                delta_peak_fr = zeros(length(responsive_idx), 1);
                baseline_idx = 1:baseline_bins;
                response_start_bin = round((g.pre_time + 0.012) / g.bin_time);
                response_end_bin = round((g.pre_time + g.monosyn_window) / g.bin_time);

                for n = 1:length(responsive_idx)
                    idx_n = responsive_idx(n);

                    % Calculate baseline firing rate
                    baseline_fr = mean(psth_Hz(idx_n, baseline_idx));

                    % Calculate peak firing rate in monosyn window
                    peak_fr = max(psth_Hz(idx_n, response_start_bin:response_end_bin));

                    % Delta Peak FR
                    delta_peak_fr(n) = peak_fr - baseline_fr;
                end

                if strcmp(brain_regions{br}, 'Astria')
                    region_name = 'AStria';
                else
                    region_name = brain_regions{br};
                end

                fprintf(fid, '%s (n = %d responsive neurons): %.2f ± %.2f Hz (mean ± SEM), median = %.2f Hz, SD = %.2f Hz\n', ...
                    region_name, length(responsive_idx), ...
                    mean(delta_peak_fr), std(delta_peak_fr)/sqrt(length(delta_peak_fr)), ...
                    median(delta_peak_fr), std(delta_peak_fr));
            end
        end

        % Perform Wilcoxon rank-sum test between LA and AStria
        if ~isempty(results_all{1}) && ~isempty(results_all{2})
            % Recalculate for both regions for the test
            LA_data = [];
            AStria_data = [];

            for br = 1:2
                res = results_all{br};

                if stim == 1
                    psth_Hz = res.psth_CS_Hz_all;
                elseif stim == 2
                    psth_Hz = res.psth_US_Hz_all;
                else
                    psth_Hz = res.psth_Both_Hz_all;
                end

                responsive_idx = find(ismember(res.Clusters_all, [1 2 3]));

                if ~isempty(responsive_idx)
                    delta_peak_fr = zeros(length(responsive_idx), 1);
                    baseline_idx = 1:baseline_bins;
                    response_start_bin = round((g.pre_time + 0.012) / g.bin_time);
                    response_end_bin = round((g.pre_time + g.monosyn_window) / g.bin_time);

                    for n = 1:length(responsive_idx)
                        idx_n = responsive_idx(n);
                        baseline_fr = mean(psth_Hz(idx_n, baseline_idx));
                        peak_fr = max(psth_Hz(idx_n, response_start_bin:response_end_bin));
                        delta_peak_fr(n) = peak_fr - baseline_fr;
                    end

                    if br == 1
                        LA_data = delta_peak_fr;
                    else
                        AStria_data = delta_peak_fr;
                    end
                end
            end

            if ~isempty(LA_data) && ~isempty(AStria_data)
                [p_val, ~] = ranksum(LA_data, AStria_data);
                fprintf(fid, 'Wilcoxon rank-sum test (LA vs AStria): p = %.4f\n', p_val);
            end
        end

        fprintf(fid, '\n');
    end

    fprintf(fid, '========================================\n');
    fprintf(fid, 'END OF FIGURE 4 STATISTICS\n');
    fprintf(fid, '========================================\n');

    fclose(fid);
    fprintf('Statistics exported to: figure_4_stats.txt\n');
end
