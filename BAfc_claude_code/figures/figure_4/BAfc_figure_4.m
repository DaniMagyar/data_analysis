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
g.artifact_exclusion = 0.010;  % Artifact exclusion period (12ms)
g.smoothvalue = 5;
g.plotwin = [0.05 0.05];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.xlinewidth = 2;
g.clim_percentile = 95;

g.onset_threshold = 5;  % Lower threshold for monosynaptic detection
g.min_consec_bins = max(1, round(0.001 / g.bin_time));  % 3ms for brief monosynaptic responses
g.alpha = 0.0;

% Responsiveness detection method
g.use_two_rule = true;  % true: two-rule (Rule 1 OR Rule 2), false: one-rule (z-score only)

% Two-rule responsiveness parameters (used if g.use_two_rule = true)
g.zscore_threshold_rule1 = 5;   % Rule 1: z-score threshold
g.prob_threshold_rule1 = 0.05;  % Rule 1: probability threshold
g.zscore_threshold_rule2 = 5;  % Rule 2: z-score threshold
g.prob_threshold_rule2 = 0.05;   % Rule 2: probability threshold

% One-rule responsiveness parameter (used if g.use_two_rule = false)
g.zscore_threshold_one_rule = 5;  % One-rule: z-score threshold only

cluster_colors = [0.8 0.2 0.2; 0.2 0.4 0.8; 0.6 0.2 0.6];  % CS-sel, US-sel, Multi

%% Calculate PSTHs once
psthZ_full = cell(1,3);
psthHz_full = cell(1,3);
baseline_bins = round(g.pre_time / g.bin_time);
filter_delay = floor(g.smoothvalue / 2);

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

%% Calculate trial probabilities once
postAP_norm_all = cell(1,3);
for hmp = 1:3
    [~, ~, postAP_norm_all{hmp}] = BAfc_psth_spx('cell_metrics', cell_metrics, 'ttl', ttl{hmp}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
end

%% Monosynaptic detection for each brain region
results_all = cell(1,2);
monosyn_window_bins = round((g.pre_time)/g.bin_time+1 : (g.pre_time+g.monosyn_window)/g.bin_time);

for br = 1:2
    idx_neurons = strcmp(cell_metrics.brainRegion, brain_regions{br});
    n_neurons = sum(idx_neurons);
    if n_neurons == 0, continue; end

    % Extract monosynaptic window responses for CS, US, CS+US
    psth_all = cell(2,3);  % (z/Hz, CS/US/Both)
    for hmp = 1:3
        psth_all{1,hmp} = psthZ_full{hmp}(idx_neurons, :);
        psth_all{2,hmp} = psthHz_full{hmp}(idx_neurons, :);
    end

    neuron_indices_all = find(idx_neurons);

    % Calculate spike probabilities
    prob_all = zeros(n_neurons, 2);  % CS, US
    for stim = 1:2
        for n = 1:n_neurons
            global_idx = neuron_indices_all(n);
            if ~isempty(postAP_norm_all{stim}{global_idx})
                responsive_trials = 0;
                for trial = 1:length(postAP_norm_all{stim}{global_idx})
                    trial_spikes = postAP_norm_all{stim}{global_idx}{trial};
                    if any(trial_spikes > 0 & trial_spikes <= g.monosyn_window)
                        responsive_trials = responsive_trials + 1;
                    end
                end
                num_trials = size(cell_metrics.general.(ttl{stim}){global_idx}, 1);
                prob_all(n, stim) = responsive_trials / num_trials;
            end
        end
    end

    % Responsiveness for all neurons
    peak_all = [max(psth_all{1,1}(:, monosyn_window_bins), [], 2), ...
                max(psth_all{1,2}(:, monosyn_window_bins), [], 2)];

    if g.use_two_rule
        responsive_all = (peak_all >= g.zscore_threshold_rule1 & prob_all >= g.prob_threshold_rule1) | ...
                         (peak_all >= g.zscore_threshold_rule2 & prob_all >= g.prob_threshold_rule2);
    else
        responsive_all = peak_all >= g.zscore_threshold_one_rule;
    end

    % Categorize neurons
    Clusters_all = zeros(n_neurons, 1);
    Clusters_all(responsive_all(:,1) & ~responsive_all(:,2)) = 1;  % CS-selective
    Clusters_all(responsive_all(:,2) & ~responsive_all(:,1)) = 2;  % US-selective
    Clusters_all(responsive_all(:,1) & responsive_all(:,2)) = 3;   % Multisensory

    % Calculate onset/offset latencies
    lat_all = nan(n_neurons, 6);  % CS_onset, CS_offset, US_onset, US_offset, Both_onset, Both_offset
    for n = 1:n_neurons
        for stim = 1:3
            [lat_all(n, stim*2-1), lat_all(n, stim*2)] = compute_onset_offset_latency(...
                psth_all{1,stim}(n, :), monosyn_window_bins, g.onset_threshold, g.min_consec_bins, g.bin_time, g.artifact_exclusion);
        end
    end

    % Sort neurons by category using rank score
    leafOrder_all = [];
    for c = [1 2 3]
        clust_idx = find(Clusters_all == c);
        if isempty(clust_idx), continue; end

        % Select onset/offset based on cluster
        lat_col = (c == 1) * 1 + (c == 2) * 3 + (c == 3) * 1;  % CS for c=1,3; US for c=2
        onset_c = lat_all(clust_idx, lat_col*2-1);
        offset_c = lat_all(clust_idx, lat_col*2);

        rank_score = onset_c + g.alpha * (offset_c - onset_c);
        [~, sort_idx] = sortrows([isnan(rank_score), rank_score], [1 2]);
        leafOrder_all = [leafOrder_all; clust_idx(sort_idx)];
    end

    % Store results
    results_all{br}.Clusters_all = Clusters_all;
    results_all{br}.leafOrder_all = leafOrder_all;
    results_all{br}.psth_CS_all = psth_all{1,1};
    results_all{br}.psth_US_all = psth_all{1,2};
    results_all{br}.psth_Both_all = psth_all{1,3};
    results_all{br}.psth_CS_Hz_all = psth_all{2,1};
    results_all{br}.psth_US_Hz_all = psth_all{2,2};
    results_all{br}.psth_Both_Hz_all = psth_all{2,3};
    results_all{br}.postAP_norm_CS = postAP_norm_all{1};
    results_all{br}.postAP_norm_US = postAP_norm_all{2};
    results_all{br}.postAP_norm_Both = postAP_norm_all{3};
    results_all{br}.CS_onset_lat_all = lat_all(:,1);
    results_all{br}.CS_offset_lat_all = lat_all(:,2);
    results_all{br}.US_onset_lat_all = lat_all(:,3);
    results_all{br}.US_offset_lat_all = lat_all(:,4);
    results_all{br}.Both_onset_lat_all = lat_all(:,5);
    results_all{br}.Both_offset_lat_all = lat_all(:,6);
    results_all{br}.n_neurons = n_neurons;
    results_all{br}.neuron_indices_all = neuron_indices_all;
    results_all{br}.animals = cell_metrics.animal(idx_neurons);
end

%% Create figure - 4x6 layout
fig = figure('Units', 'pixels', 'Position', [100, 100, 1200, 800], 'Visible', 'on');
t = tiledlayout(fig, 4, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add panel labels
label_positions = {[0.01 0.94], [0.48 0.94], [0.01 0.47], [0.48 0.47], ...
                   [0.68 0.94], [0.68 0.62], [0.68 0.28]};
label_texts = {'A', 'C', 'B', 'D', 'E', 'F', 'G'};
for i = 1:length(label_positions)
    annotation(fig, 'textbox', [label_positions{i} 0.05 0.05], 'String', label_texts{i}, ...
        'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

% Determine global color limits
all_values = [];
for br = 1:2
    if ~isempty(results_all{br})
        res = results_all{br};
        Clusters_sorted = res.Clusters_all(res.leafOrder_all);
        responsive_mask = ismember(Clusters_sorted, [1, 2, 3]);
        responsive_leafOrder = res.leafOrder_all(responsive_mask);
        plot_bins = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);

        psth_data = {res.psth_CS_all, res.psth_US_all, res.psth_Both_all};
        for stim = 1:3
            matrix_data = psth_data{stim}(responsive_leafOrder, plot_bins);
            all_values = [all_values; matrix_data(:)];
        end
    end
end
clim_min = prctile(all_values, (100 - g.clim_percentile) / 2);
clim_max = prctile(all_values, 100 - (100 - g.clim_percentile) / 2);

%% Generate heatmaps + bar charts
t_heatmaps = tiledlayout(t, 2, 4, 'TileSpacing', 'tight', 'Padding', 'tight');
t_heatmaps.Layout.Tile = 1;
t_heatmaps.Layout.TileSpan = [4 4];

kw_data_storage = cell(2, 3, 3);  % region × cluster × stimulus

for br = 1:2
    if isempty(results_all{br}), continue; end
    res = results_all{br};

    % Get responsive neurons
    Clusters_sorted = res.Clusters_all(res.leafOrder_all);
    responsive_mask = ismember(Clusters_sorted, [1, 2, 3]);
    responsive_leafOrder = res.leafOrder_all(responsive_mask);
    Clusters_sorted_resp = Clusters_sorted(responsive_mask);
    n_clu = find(diff(Clusters_sorted_resp) ~= 0);

    % Heatmaps (columns 1-3)
    psth_data = {res.psth_CS_all, res.psth_US_all, res.psth_Both_all};
    lat_data = {[res.CS_onset_lat_all, res.CS_offset_lat_all], ...
                [res.US_onset_lat_all, res.US_offset_lat_all], ...
                [res.Both_onset_lat_all, res.Both_offset_lat_all]};
    stim_titles = {'CS', 'US', 'CS+US'};

    for stim = 1:3
        ax = nexttile(t_heatmaps, (br-1)*4 + stim);
        plot_bins = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        matrix = psth_data{stim}(responsive_leafOrder, plot_bins);

        imagesc(g.timeaxis_hmp, 1:size(matrix,1), matrix);
        clim([clim_min clim_max]);
        colormap(ax, g.colors.Heatmap);
        hold on;

        % Artifact overlay
        patch([0 g.artifact_exclusion g.artifact_exclusion 0], ...
              [0.5 0.5 size(matrix,1)+0.5 size(matrix,1)+0.5], ...
              'k', 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        if stim >= 2
            text(0.005, 2, '⚡', 'Color', 'y', 'FontSize', g.fontSize1, ...
                 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        end

        xline(0, '--k', 'LineWidth', g.xlinewidth);
        xline(g.monosyn_window, '-k', 'LineWidth', 1);

        % Labels
        if br == 1, title(stim_titles{stim}, 'FontSize', g.fontSize1, 'FontWeight', 'bold'); end
        if stim == 1
            if strcmp(brain_regions{br}, 'Astria')
                ylabel_text = 'AStria neurons';
            else
                ylabel_text = [brain_regions{br} ' neurons'];
            end
            ylabel(ylabel_text, 'FontSize', g.fontSize2, 'FontWeight', 'bold');
        else
            set(gca, 'YTickLabel', []);
        end
        if br == 2, xlabel('Time (s)', 'FontSize', g.fontSize2); else, set(gca, 'XTickLabel', []); end
        if size(matrix,1) > 1
            yticks([1, size(matrix,1)]);
        else
            yticks(1);
        end

        % Cluster boundaries
        for i = 1:length(n_clu), yline(n_clu(i) + 0.5, 'k-', 'LineWidth', 1); end

        % Onset/offset markers
        for n = 1:length(responsive_leafOrder)
            idx_n = responsive_leafOrder(n);
            onset_lat = lat_data{stim}(idx_n, 1);
            offset_lat = lat_data{stim}(idx_n, 2);
            if ~isnan(onset_lat), plot(onset_lat, n, 'k.', 'MarkerSize', 4); end
            if ~isnan(offset_lat), plot(offset_lat, n, 'k.', 'MarkerSize', 4); end
        end
        hold off;
        set(gca, 'FontSize', g.fontSize2);
    end

    % Bar charts (column 4)
    [kw_data_storage, ~] = plot_delta_fr_bars(t_heatmaps, (br-1)*4 + 4, res, g, cluster_colors, kw_data_storage, br);
end

% Colorbar
drawnow;
cb_ax = axes('Position', [0.185, 0.06, 0.008, 0.08]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal', 'XTick', [], 'YAxisLocation', 'right', 'FontSize', g.fontSize2);

%% Right panels
t_right = tiledlayout(t, 3, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
t_right.Layout.Tile = 5;
t_right.Layout.TileSpan = [4 2];

% Pie charts
[t_pie, contingency_table] = plot_pie_charts(t_right, results_all, cluster_colors, brain_regions, g);

% Region comparison
plot_region_comparison(t_right, results_all, g, 'ΔFR (Hz)', 3);

% Latency comparison
plot_latency_comparison(t_right, results_all, g);

%% Export statistics
kw_results = compute_friedman_tests(kw_data_storage);
[chi2_obs, p_perm, cramers_v] = compute_chi_square_stats(contingency_table, results_all);

export_figure4_to_excel_simple(results_all, kw_data_storage, kw_results, contingency_table, ...
    chi2_obs, p_perm, cramers_v, brain_regions, {'CS-sel', 'US-sel', 'Multi'}, cell_metrics, g, 'figure_4_data.xlsx');
exportgraphics(gcf, 'figure_4.png', 'Resolution', 300);

%% Helper functions
function [onset_lat, offset_lat] = compute_onset_offset_latency(z_trace, event_inds, threshold, min_consec, bin_time, artifact_exclusion)
    seg = z_trace(event_inds);
    if any(isnan(seg)) || length(seg) <= round(artifact_exclusion / bin_time)
        onset_lat = NaN; offset_lat = NaN; return
    end

    artifact_bins = round(artifact_exclusion / bin_time);
    seg_post = seg(artifact_bins+1:end);
    isAbove = seg_post >= threshold;
    winsum = conv(double(isAbove), ones(1, min_consec), 'same');
    onset_idx_rel = find(winsum >= min_consec, 1, 'first');

    if isempty(onset_idx_rel)
        onset_lat = NaN; offset_lat = NaN;
    else
        onset_idx = artifact_bins + onset_idx_rel;
        onset_lat = (onset_idx - 1) * bin_time;

        isBelow = seg(onset_idx:end) < threshold;
        winsum_below = conv(double(isBelow), ones(1, min_consec), 'same');
        offset_idx_rel = find(winsum_below >= min_consec, 1, 'first');

        if isempty(offset_idx_rel)
            offset_lat = (length(event_inds) - 1) * bin_time;
        else
            offset_lat = (onset_idx + offset_idx_rel - 2) * bin_time;
        end
    end
end

function [kw_data, kw_results] = plot_delta_fr_bars(t_parent, tile_num, res, g, cluster_colors, kw_data, br)
    t_bars = tiledlayout(t_parent, 3, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    t_bars.Layout.Tile = tile_num;

    window_duration = g.monosyn_window - g.artifact_exclusion;
    postAP_data = {res.postAP_norm_CS, res.postAP_norm_US, res.postAP_norm_Both};

    for c = 1:3
        ax_bar = nexttile(t_bars, c);
        clust_idx = find(res.Clusters_all == c);

        if ~isempty(clust_idx)
            metrics = zeros(length(clust_idx), 3);
            for stim = 1:3
                metrics(:, stim) = calculate_delta_fr(clust_idx, res.neuron_indices_all, ...
                    postAP_data{stim}, g.artifact_exclusion, g.monosyn_window, window_duration);
            end

            for stim = 1:3, kw_data{br, c, stim} = metrics(:, stim); end

            means_data = mean(metrics, 1, 'omitnan');
            sems_data = std(metrics, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(metrics(:,1))));

            hold on;
            bar([1 2 3], means_data, 0.4, 'FaceColor', cluster_colors(c, :), 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([1 2 3], means_data, sems_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

            % Statistics
            p_vals = ones(3, 1);
            if length(clust_idx) >= 2
                [p_friedman, ~, ~] = friedman(metrics, 1, 'off');
                if p_friedman < 0.05
                    p_vals(1) = signrank(metrics(:,1), metrics(:,2));
                    p_vals(2) = signrank(metrics(:,1), metrics(:,3));
                    p_vals(3) = signrank(metrics(:,2), metrics(:,3));
                end
            end

            % Significance markers
            y_max = max(means_data + sems_data);
            y_range = 80;
            y_level1 = y_max + 0.08 * y_range;
            y_level2 = y_level1 + 10;

            plot_sig_marker([1.1 1.9], y_level1, p_vals(1), g.fontSize2, y_range);
            plot_sig_marker([2.1 2.9], y_level1, p_vals(3), g.fontSize2, y_range);
            plot_sig_marker([1 3], y_level2, p_vals(2), g.fontSize2, y_range);

            hold off;
        end

        xlim([0.5 3.5]); xticks([1 2 3]); ylim([0 80]); yticks([0 40 80]);
        if br == 1 && c == 1, title('FR (Hz)', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Interpreter', 'tex'); end
        if br == 2 && c == 3, xticklabels({'CS', 'US', 'CS+US'}); else, set(gca, 'XTickLabel', []); end
        set(gca, 'FontSize', g.fontSize2); box off;
    end
    kw_results = [];
end

function delta_fr = calculate_delta_fr(clust_idx, neuron_indices, postAP_norm, artifact_excl, monosyn_win, win_dur)
    delta_fr = zeros(length(clust_idx), 1);
    for n = 1:length(clust_idx)
        global_idx = neuron_indices(clust_idx(n));
        if ~isempty(postAP_norm{global_idx})
            spike_count = 0;
            for trial = 1:length(postAP_norm{global_idx})
                spikes = postAP_norm{global_idx}{trial};
                spike_count = spike_count + sum(spikes >= artifact_excl & spikes <= monosyn_win);
            end
            delta_fr(n) = (spike_count / length(postAP_norm{global_idx})) / win_dur;
        end
    end
end

function plot_sig_marker(x_range, y_pos, p_val, font_size, y_range)
    sig_text = '';
    if p_val < 0.001, sig_text = '***';
    elseif p_val < 0.01, sig_text = '**';
    elseif p_val < 0.05, sig_text = '*';
    end
    if ~isempty(sig_text)
        plot(x_range, [y_pos y_pos], 'k-', 'LineWidth', 1.5);
        text(mean(x_range), y_pos + 0.02 * y_range, sig_text, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', font_size);
    end
end

function [t_pie, contingency_table] = plot_pie_charts(t_parent, results_all, cluster_colors, brain_regions, g)
    t_pie = tiledlayout(t_parent, 1, 2, 'TileSpacing', 'none', 'Padding', 'tight');
    t_pie.Layout.Tile = 1;
    t_pie.Layout.TileSpan = [1 2];

    contingency_table = zeros(2, 3);
    region_names = {'LA', 'AStria'};

    for br = 1:2
        if isempty(results_all{br}), continue; end

        nexttile(t_pie, br);
        res = results_all{br};
        counts = [sum(res.Clusters_all == 1), sum(res.Clusters_all == 2), sum(res.Clusters_all == 3)];
        contingency_table(br, :) = counts;
        total_n = sum(counts);

        pie_labels = arrayfun(@(x) sprintf('%.0f%%', x/total_n*100), counts, 'UniformOutput', false);
        p = pie(counts, pie_labels);

        for i = 1:2:length(p)
            set(p(i), 'FaceColor', cluster_colors((i+1)/2, :), 'EdgeColor', 'k', 'LineWidth', 1);
        end
        for i = 2:2:length(p)
            set(p(i), 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'w');
            pos = get(p(i), 'Position');
            set(p(i), 'Position', pos * (0.3 + 0.03*br));
        end

        title(region_names{br}, 'FontSize', 12, 'FontWeight', 'bold');
    end

    % Legend
    ax_legend = axes('Position', [0.65 0.85 0.05 0.05], 'Visible', 'off');
    hold(ax_legend, 'on');
    h = arrayfun(@(i) plot(ax_legend, NaN, NaN, 's', 'MarkerSize', 12, ...
        'MarkerFaceColor', cluster_colors(i, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1), 1:3);
    lgd = legend(ax_legend, h, {'CS-sel', 'US-sel', 'Multi'}, 'Location', 'east', 'FontSize', 10, 'Box', 'off');
    lgd.ItemTokenSize = [30, 30];

    title(t_pie, 'Response categories', 'FontSize', 12, 'FontWeight', 'bold');
end

function plot_region_comparison(t_parent, results_all, g, ylabel_text, tile_start)
    t_comp = tiledlayout(t_parent, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
    t_comp.Layout.Tile = tile_start;
    t_comp.Layout.TileSpan = [1 2];

    stim_names = {'CS', 'US', 'CS+US'};
    window_duration = g.monosyn_window - g.artifact_exclusion;

    for stim = 1:3
        nexttile(t_comp, stim);

        data_all = [];
        labels_all = {};

        for br = 1:2
            if isempty(results_all{br}), continue; end
            res = results_all{br};

            % Select responsive neurons based on stimulus
            % CS: CS-sel + Multi (clusters 1, 3)
            % US: US-sel + Multi (clusters 2, 3)
            % CS+US: All responsive (clusters 1, 2, 3)
            if stim == 1
                responsive_idx = find(ismember(res.Clusters_all, [1 3]));
                postAP_data = res.postAP_norm_CS;
            elseif stim == 2
                responsive_idx = find(ismember(res.Clusters_all, [2 3]));
                postAP_data = res.postAP_norm_US;
            else
                responsive_idx = find(ismember(res.Clusters_all, [1 2 3]));
                postAP_data = res.postAP_norm_Both;
            end

            if ~isempty(responsive_idx)
                delta_fr = calculate_delta_fr(responsive_idx, res.neuron_indices_all, ...
                    postAP_data, g.artifact_exclusion, g.monosyn_window, window_duration);
                data_all = [data_all; delta_fr];
                if br == 1
                    labels_all = [labels_all; repmat({'LA'}, length(responsive_idx), 1)];
                else
                    labels_all = [labels_all; repmat({'AStria'}, length(responsive_idx), 1)];
                end
            end
        end

        if ~isempty(data_all)
            LA_data = data_all(strcmp(labels_all, 'LA'));
            Astria_data = data_all(strcmp(labels_all, 'AStria'));

            means_data = [mean(LA_data, 'omitnan'), mean(Astria_data, 'omitnan')];
            sems_data = [std(LA_data, 0, 'omitnan') / sqrt(length(LA_data)), ...
                         std(Astria_data, 0, 'omitnan') / sqrt(length(Astria_data))];

            hold on;
            bar_colors = [0.7 0.2 0.2; 0.2 0.4 0.7];
            b = bar([1 2], means_data, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1);
            b.CData = bar_colors;

            scatter(1 + (rand(length(LA_data), 1) - 0.5) * 0.15, LA_data, 16, [0.5 0.5 0.5], 'LineWidth', 0.5);
            scatter(2 + (rand(length(Astria_data), 1) - 0.5) * 0.15, Astria_data, 16, [0.5 0.5 0.5], 'LineWidth', 0.5);
            errorbar([1 2], means_data, sems_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

            [p_val, ~] = ranksum(LA_data, Astria_data);
            y_max = max(max(means_data + sems_data), prctile([LA_data; Astria_data], 95));
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
            text(1.5, y_pos, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
            hold off;

            ylim([0 y_range]);
        else
            ylim([0 250]);
        end

        xlim([0.5 2.5]); xticks([1 2]); xticklabels({'LA', 'AStria'});
        if stim == 1, ylabel(ylabel_text, 'FontSize', g.fontSize2); else, set(gca, 'YTickLabel', []); end
        title(stim_names{stim}, 'FontSize', g.fontSize1, 'FontWeight', 'bold');
        set(gca, 'FontSize', g.fontSize2); box off;
    end

    title(t_comp, 'Across region comparison', 'FontSize', 12, 'FontWeight', 'bold');
end

function plot_latency_comparison(t_parent, results_all, g)
    t_latency = tiledlayout(t_parent, 1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');
    t_latency.Layout.Tile = 5;
    t_latency.Layout.TileSpan = [1 2];

    stim_names = {'CS', 'US', 'CS+US'};
    region_colors = {[0.8 0.2 0.2], [0.2 0.4 0.8]};

    for stim = 1:3
        nexttile(t_latency, stim);

        data_all = [];
        labels_all = {};

        for br = 1:2
            if isempty(results_all{br}), continue; end
            res = results_all{br};

            % Select latencies based on stimulus
            if stim == 1
                responsive_idx = find(ismember(res.Clusters_all, [1 3]));
                onset_lat = res.CS_onset_lat_all(responsive_idx);
            elseif stim == 2
                responsive_idx = find(ismember(res.Clusters_all, [2 3]));
                onset_lat = res.US_onset_lat_all(responsive_idx);
            else
                responsive_idx = find(ismember(res.Clusters_all, [1 2 3]));
                onset_lat = res.Both_onset_lat_all(responsive_idx);
            end

            onset_lat = onset_lat(~isnan(onset_lat)) * 1000;
            if ~isempty(onset_lat)
                data_all = [data_all; onset_lat];
                if br == 1
                    labels_all = [labels_all; repmat({'LA'}, length(onset_lat), 1)];
                else
                    labels_all = [labels_all; repmat({'AStria'}, length(onset_lat), 1)];
                end
            end
        end

        if ~isempty(data_all)
            LA_data = data_all(strcmp(labels_all, 'LA'));
            Astria_data = data_all(strcmp(labels_all, 'AStria'));

            hold on;
            boxplot([LA_data; Astria_data], [ones(length(LA_data),1); 2*ones(length(Astria_data),1)], ...
                'Positions', [1 2], 'Widths', 0.6, 'Colors', 'k');

            if ~isempty(LA_data)
                scatter(1 + (rand(length(LA_data), 1) - 0.5) * 0.15, LA_data, 20, region_colors{1}, 'filled', 'MarkerFaceAlpha', 0.4);
            end
            if ~isempty(Astria_data)
                scatter(2 + (rand(length(Astria_data), 1) - 0.5) * 0.15, Astria_data, 20, region_colors{2}, 'filled', 'MarkerFaceAlpha', 0.4);
            end

            if ~isempty(LA_data) && ~isempty(Astria_data)
                [p_val, ~] = ranksum(LA_data, Astria_data);
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
                text(1.5, y_pos, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', g.fontSize2);
            end
            hold off;
        end

        xlim([0.5 2.5]); xticks([1 2]); xticklabels({'LA', 'AStria'});
        ylim([0 50]); yticks([0 25 50]);
        if stim == 1, ylabel('Onset Latency (ms)', 'FontSize', g.fontSize2); else, set(gca, 'YTickLabel', []); end
        title(stim_names{stim}, 'FontSize', g.fontSize1, 'FontWeight', 'bold');
        set(gca, 'FontSize', g.fontSize2); box off;
    end

    title(t_latency, 'Onset latency comparison', 'FontSize', 12, 'FontWeight', 'bold');
end

function kw_results = compute_friedman_tests(kw_data_storage)
    kw_results = struct();
    for br = 1:2
        for c = 1:3
            if isempty(kw_data_storage{br, c, 1}) || length(kw_data_storage{br, c, 1}) < 2
                kw_results(br, c).p_friedman = 1;
                kw_results(br, c).p_values = ones(3, 1);
                continue;
            end

            data_matrix = [kw_data_storage{br, c, 1}(:), kw_data_storage{br, c, 2}(:), kw_data_storage{br, c, 3}(:)];
            [p_friedman, ~, ~] = friedman(data_matrix, 1, 'off');

            if p_friedman < 0.05
                p_values = zeros(3, 1);
                p_values(1) = signrank(kw_data_storage{br, c, 1}, kw_data_storage{br, c, 2});
                p_values(2) = signrank(kw_data_storage{br, c, 1}, kw_data_storage{br, c, 3});
                p_values(3) = signrank(kw_data_storage{br, c, 2}, kw_data_storage{br, c, 3});
            else
                p_values = ones(3, 1);
            end

            kw_results(br, c).p_friedman = p_friedman;
            kw_results(br, c).p_values = p_values;
        end
    end
end

function [chi2_obs, p_perm, cramers_v] = compute_chi_square_stats(contingency_table, results_all)
    [chi2_obs, ~] = calculate_chi_square(contingency_table);
    n_total = sum(contingency_table(:));
    df_cramer = min(size(contingency_table, 1) - 1, size(contingency_table, 2) - 1);
    cramers_v = sqrt(chi2_obs / (n_total * df_cramer));

    % Permutation test
    all_clusters = [];
    region_indices = [];
    for br = 1:2
        if ~isempty(results_all{br})
            all_clusters = [all_clusters; results_all{br}.Clusters_all];
            region_indices = [region_indices; br * ones(results_all{br}.n_neurons, 1)];
        end
    end

    n_permutations = 10000;
    chi2_perm = zeros(n_permutations, 1);
    for perm = 1:n_permutations
        shuffled = all_clusters(randperm(length(all_clusters)));
        perm_table = zeros(2, 3);
        for br = 1:2
            br_clusters = shuffled(region_indices == br);
            for c = 1:3
                perm_table(br, c) = sum(br_clusters == c);
            end
        end
        chi2_perm(perm) = calculate_chi_square(perm_table);
    end
    p_perm = sum(chi2_perm >= chi2_obs) / n_permutations;
end

function [chi2_stat, p_value] = calculate_chi_square(contingency_table)
    row_totals = sum(contingency_table, 2);
    col_totals = sum(contingency_table, 1);
    n_total = sum(contingency_table(:));
    expected = (row_totals * col_totals) / n_total;
    chi2_stat = sum(((contingency_table(:) - expected(:)).^2) ./ expected(:));
    df = (size(contingency_table, 1) - 1) * (size(contingency_table, 2) - 1);
    p_value = 1 - chi2cdf(chi2_stat, df);
end
