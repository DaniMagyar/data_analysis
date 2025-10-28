% BAfc_figure_3
% Compare CS-selective, US-selective, and Multisensory neurons identified from CS and US
% Show their responses to CS+US (both) stimulus
% 4 brain regions (LA, BA, Astria, CeA) with heatmaps showing CS+US responses

clear all; close all

%% Setup
recordings = {...
    'MD292_002_kilosort',...
    'MD293_001_kilosort',...
    'MD294_001_kilosort',...
    'MD295_001_kilosort',...
    'MD296_001_kilosort',...
    'MD297_001_kilosort',...
    'MD298_001_kilosort',...
    'MD299_001_kilosort',...
    'MD300_001_kilosort',...
    'MD304_001_kilosort',...
    'MD305_001_kilosort',...
    'MD307_001_kilosort',...
    'MD309_001_kilosort',...
    'MD310_001_kilosort',...
    'MD311_002_kilosort',...
    'MD312_001_kilosort',...
    'MD313_001_kilosort',...
    'MD314_001_kilosort',...
    'MD315_001_kilosort',...
    'MD316_002_kilosort',...
    'MD317_001_kilosort',...
    'MD318_001_kilosort',...
    'MD318_002_kilosort',...
    'MD319_003_kilosort'};

ttl = {'triptest_sound_only','triptest_shocks_only','triptest_both'};
hmptitles = {'CS', 'US', 'CS+US'};

g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
g.cell_metrics = BAfc_match_celltypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 14;
g.fontSize2 = 12;
g.xlinewidth = 2;
g.pre_time = 5;
g.post_time = 4;
g.test_time = 1;

% Clustering Parameters
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

% Brain regions and cell types to analyze
brain_regions = {'LA', 'Astria'};
cell_type_filter = {'all', 'all'};

cluster_colors = [
    0.8 0.2 0.2;    % CS-selective
    0.2 0.4 0.8;    % US-selective
    0.6 0.2 0.6;    % Multisensory
];

cluster_names = {'CS-sel', 'US-sel', 'Multi'};

%% Calculate PSTHs once
fprintf('Calculating PSTHs...\n');
psthZ_full = cell(1, numel(ttl));
psthHz_full = cell(1, numel(ttl));

for hmp = 1:numel(ttl)
    psth_spx_og = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Convert to Hz: divide by number of trials and bin time
    num_trials = size(g.cell_metrics.general.(ttl{hmp}){1}, 1);
    psth_hz = psth_spx_og / (num_trials * g.bin_time);

    % Smooth Hz data
    psth_hz_smooth = smoothdata(psth_hz, 2, 'sgolay', g.smoothvalue);
    psthHz_full{hmp} = psth_hz_smooth;

    baseline_idx = 1:(g.pre_time / g.bin_time);
    baseline_mean = mean(psth_spx_og(:, baseline_idx), 2);
    baseline_std = std(psth_spx_og(:, baseline_idx), 0, 2);
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;
    psth_spx = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);
    psthZ_full{hmp} = psth_spx;
end

%% Process each brain region
results_all = cell(1, 2);

for br = 1:2
    fprintf('\nProcessing %s...\n', brain_regions{br});

    % Get neuron indices
    if strcmp(cell_type_filter{br}, 'PN')
        idx_neurons = strcmp(g.cell_metrics.brainRegion, brain_regions{br}) & strcmp(g.cell_metrics.putativeCellType, 'PN');
    else
        idx_neurons = strcmp(g.cell_metrics.brainRegion, brain_regions{br});
    end

    n_neurons = sum(idx_neurons);
    fprintf('  %d neurons\n', n_neurons);

    if n_neurons == 0
        continue;
    end

    % Extract PSTHs for CS, US, and CS+US
    psth_CS = psthZ_full{1}(idx_neurons, :);
    psth_US = psthZ_full{2}(idx_neurons, :);
    psth_Both = psthZ_full{3}(idx_neurons, :);
    psth_CS_Hz = psthHz_full{1}(idx_neurons, :);
    psth_US_Hz = psthHz_full{2}(idx_neurons, :);
    psth_Both_Hz = psthHz_full{3}(idx_neurons, :);

    % Calculate responses based on CS and US only (not Both)
    CS_peak = max(psth_CS(:, g.roi), [], 2);
    US_peak = max(psth_US(:, g.roi), [], 2);

    baseline_idx = 1:(g.pre_time / g.bin_time);
    CS_baseline_fr = mean(psth_CS_Hz(:, baseline_idx), 2);
    US_baseline_fr = mean(psth_US_Hz(:, baseline_idx), 2);
    CS_test_fr = mean(psth_CS_Hz(:, g.roi), 2);
    US_test_fr = mean(psth_US_Hz(:, g.roi), 2);
    CS_fr_drop = (CS_baseline_fr - CS_test_fr) ./ (CS_baseline_fr + eps);
    US_fr_drop = (US_baseline_fr - US_test_fr) ./ (US_baseline_fr + eps);

    % Classify neurons based on CS and US responses only (same as figure 2)
    CS_excited = CS_peak >= g.excitation_threshold;
    US_excited = US_peak >= g.excitation_threshold;
    CS_inhibited = CS_fr_drop >= g.inhibition_fr_drop;
    US_inhibited = US_fr_drop >= g.inhibition_fr_drop;

    Clusters = zeros(n_neurons, 1);
    Clusters(CS_excited & ~US_excited) = 1;  % CS-selective
    Clusters(US_excited & ~CS_excited) = 2;  % US-selective
    Clusters(CS_excited & US_excited) = 3;   % Multisensory
    Clusters(~CS_excited & ~US_excited & (CS_inhibited | US_inhibited)) = 5;  % Inhibited
    Clusters(~CS_excited & ~US_excited & ~CS_inhibited & ~US_inhibited) = 4;  % Non-responsive

    % Compute latencies for sorting (using CS and US) - for ALL neurons
    CS_onset_lat = nan(n_neurons, 1);
    CS_offset_lat = nan(n_neurons, 1);
    US_onset_lat = nan(n_neurons, 1);
    US_offset_lat = nan(n_neurons, 1);

    for n = 1:n_neurons
        [CS_onset_lat(n), CS_offset_lat(n)] = compute_onset_offset_latency(psth_CS(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [US_onset_lat(n), US_offset_lat(n)] = compute_onset_offset_latency(psth_US(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
    end

    % Sort neurons (same as figure 2)
    leafOrder = [];
    cluster_order = [1, 2, 3, 4, 5];

    for c = cluster_order
        clust_idx = find(Clusters == c);
        if isempty(clust_idx)
            continue;
        end

        if c == 1  % CS-selective
            onset_c = CS_onset_lat(clust_idx);
            offset_c = CS_offset_lat(clust_idx);
        elseif c == 2  % US-selective
            onset_c = US_onset_lat(clust_idx);
            offset_c = US_offset_lat(clust_idx);
        elseif c == 3  % Multisensory
            CS_duration = CS_offset_lat(clust_idx) - CS_onset_lat(clust_idx);
            US_duration = US_offset_lat(clust_idx) - US_onset_lat(clust_idx);
            CS_rank_score = CS_onset_lat(clust_idx) + g.alpha * CS_duration;
            US_rank_score = US_onset_lat(clust_idx) + g.alpha * US_duration;
            rank_score = mean([CS_rank_score, US_rank_score], 2, 'omitnan');
            sort_matrix = [isnan(rank_score), rank_score];
            [~, sort_idx] = sortrows(sort_matrix, [1 2]);
            leafOrder = [leafOrder; clust_idx(sort_idx)];
            continue;
        elseif c == 4  % Non-responsive
            mean_zscore = mean([mean(psth_CS(clust_idx, g.roi), 2), mean(psth_US(clust_idx, g.roi), 2)], 2);
            [~, sort_idx] = sort(mean_zscore, 'descend');
            leafOrder = [leafOrder; clust_idx(sort_idx)];
            continue;
        else  % Inhibited
            mean_zscore = mean([mean(psth_CS(clust_idx, g.roi), 2), mean(psth_US(clust_idx, g.roi), 2)], 2);
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
    results_all{br}.psth_CS = psth_CS;
    results_all{br}.psth_US = psth_US;
    results_all{br}.psth_Both = psth_Both;
    results_all{br}.psth_CS_Hz = psth_CS_Hz;
    results_all{br}.psth_US_Hz = psth_US_Hz;
    results_all{br}.psth_Both_Hz = psth_Both_Hz;
    results_all{br}.n_neurons = n_neurons;
end

%% Create figure
fig = figure('Units', 'pixels', 'Position', [-100, -100, 1500, 3000], 'Visible', 'on');
t = tiledlayout(fig, 7, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Determine global color limits across all heatmaps (CS, US, and CS+US)
all_values = [];
for br = 1:2
    if ~isempty(results_all{br})
        for stim = 1:3
            if stim == 1
                psth_sorted = results_all{br}.psth_CS(results_all{br}.leafOrder, :);
            elseif stim == 2
                psth_sorted = results_all{br}.psth_US(results_all{br}.leafOrder, :);
            else
                psth_sorted = results_all{br}.psth_Both(results_all{br}.leafOrder, :);
            end
            matrix = psth_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
            all_values = [all_values; matrix(:)];
        end
    end
end

clim_min = prctile(all_values, (100 - g.clim_percentile) / 2);
clim_max = prctile(all_values, 100 - (100 - g.clim_percentile) / 2);

%% Plot each brain region (rows) × heatmaps+lineplots (columns)
for br = 1:2
    if isempty(results_all{br})
        continue;
    end

    res = results_all{br};

    % Get cluster boundaries for responsive clusters only (1, 2, 3)
    Clusters_sorted = res.Clusters(res.leafOrder);

    % Filter to only plot clusters 1, 2, 3
    responsive_mask = ismember(Clusters_sorted, [1, 2, 3]);
    responsive_leafOrder = res.leafOrder(responsive_mask);
    Clusters_sorted_responsive = Clusters_sorted(responsive_mask);
    n_clu = find(diff(Clusters_sorted_responsive) ~= 0);

    % Heatmaps in column 1 (nested layout for 3 stimuli)
    t_heatmaps = tiledlayout(t, 1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
    if br == 1
        t_heatmaps.Layout.Tile = 1;  % LA: row 1, col 1
    else
        t_heatmaps.Layout.Tile = 5;  % Astria: row 3, col 1
    end
    t_heatmaps.Layout.TileSpan = [2 1];

    for stim = 1:3  % CS, US, CS+US
        % Select appropriate PSTH and filter to responsive neurons only
        if stim == 1
            psth_sorted = res.psth_CS(responsive_leafOrder, :);
            stim_title = 'CS';
        elseif stim == 2
            psth_sorted = res.psth_US(responsive_leafOrder, :);
            stim_title = 'US';
        else
            psth_sorted = res.psth_Both(responsive_leafOrder, :);
            stim_title = 'CS+US';
        end

        % Heatmap
        ax = nexttile(t_heatmaps, stim);
        matrix = psth_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        imagesc(g.timeaxis_hmp, 1:size(matrix, 1), matrix);
        clim([clim_min clim_max]);
        colormap(ax, g.colors.Heatmap);
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        % Labels
        if br == 1
            title(stim_title, 'FontSize', g.fontSize1);
        end

        if stim == 1 && br == 1
            ylabel('Neuron #', 'FontSize', g.fontSize2);
        end

        % Set yticks to first and last
        n_neurons_plot = size(matrix, 1);
        yticks([1, n_neurons_plot]);

        if br == 2
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
        end

        % Add cluster lines (black)
        hold on;
        for i = 1:length(n_clu)
            yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
        end
        hold off;

        % Add brain region label (only on first column)
        if stim == 1
            text(-0.9, sum(responsive_mask)/2, brain_regions{br}, 'FontSize', g.fontSize1+2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);
        end
    end

    %% Add lineplots for each stimulus in column 2 (right side)
    % Nested tiledlayout for CS, US, CS+US lineplots
    t_outer = tiledlayout(t, 1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
    if br == 1
        t_outer.Layout.Tile = 2;  % LA: row 1, col 2
    else
        t_outer.Layout.Tile = 6;  % Astria: row 3, col 2
    end
    t_outer.Layout.TileSpan = [2 1];

    % CS lineplots
    t_nested_CS = tiledlayout(t_outer, 3, 1, 'TileSpacing', 'none', 'Padding', 'tight');
    t_nested_CS.Layout.Tile = 1;
    if br == 1
        t_nested_CS.Title.String = 'CS';
        t_nested_CS.Title.FontSize = g.fontSize1;
        t_nested_CS.Title.FontWeight = 'bold';
    end

    % US lineplots
    t_nested_US = tiledlayout(t_outer, 3, 1, 'TileSpacing', 'none', 'Padding', 'tight');
    t_nested_US.Layout.Tile = 2;
    if br == 1
        t_nested_US.Title.String = 'US';
        t_nested_US.Title.FontSize = g.fontSize1;
        t_nested_US.Title.FontWeight = 'bold';
    end

    % CS+US lineplots
    t_nested_Both = tiledlayout(t_outer, 3, 1, 'TileSpacing', 'none', 'Padding', 'tight');
    t_nested_Both.Layout.Tile = 3;
    if br == 1
        t_nested_Both.Title.String = 'CS+US';
        t_nested_Both.Title.FontSize = g.fontSize1;
        t_nested_Both.Title.FontWeight = 'bold';
    end

    plot_idx = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    time_vec = g.timeaxis_hmp;

    % Ensure consistent lengths
    if length(time_vec) ~= length(plot_idx)
        time_vec = linspace(-g.plotwin(1), g.plotwin(2), length(plot_idx));
    end

    % Calculate global y-limits for this brain region (only clusters 1-3)
    y_max = 0;
    y_min = 0;
    for c = [1 2 3]
        clust_idx = find(res.Clusters == c);
        if ~isempty(clust_idx)
            psth_CS_mean = mean(res.psth_CS_Hz(clust_idx, plot_idx), 1);
            psth_US_mean = mean(res.psth_US_Hz(clust_idx, plot_idx), 1);
            psth_Both_mean = mean(res.psth_Both_Hz(clust_idx, plot_idx), 1);
            y_max = max([y_max, max(psth_CS_mean), max(psth_US_mean), max(psth_Both_mean)]);
            y_min = min([y_min, min(psth_CS_mean), min(psth_US_mean), min(psth_Both_mean)]);
        end
    end

    fprintf('  %s lineplots: y_min = %.2f, y_max = %.2f\n', brain_regions{br}, y_min*1.1, y_max*1.1);

    for c = [1 2 3]
        clust_idx = find(res.Clusters == c);

        % CS panel
        ax_CS = nexttile(t_nested_CS, c);
        hold on;
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        if ~isempty(clust_idx)
            psth_CS_mean = mean(res.psth_CS_Hz(clust_idx, plot_idx), 1);
            plot(time_vec, psth_CS_mean, '-', 'Color', cluster_colors(c, :), 'LineWidth', 2.5);
        end

        hold off;
        xlim([time_vec(1) time_vec(end)]);
        ylim([y_min*1.1 y_max*1.1]);

        if br == 2 && c == 3
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            set(gca, 'XColor', 'k');
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
            set(gca, 'XColor', 'none');
        end

        set(gca, 'YTickLabel', []);
        set(gca, 'YColor', 'none');
        set(gca, 'FontSize', g.fontSize2);

        % US panel
        ax_US = nexttile(t_nested_US, c);
        hold on;
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        if ~isempty(clust_idx)
            psth_US_mean = mean(res.psth_US_Hz(clust_idx, plot_idx), 1);
            plot(time_vec, psth_US_mean, '-', 'Color', cluster_colors(c, :), 'LineWidth', 2.5);
        end

        hold off;
        xlim([time_vec(1) time_vec(end)]);
        ylim([y_min*1.1 y_max*1.1]);

        if br == 2 && c == 3
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            set(gca, 'XColor', 'k');
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
            set(gca, 'XColor', 'none');
        end

        set(gca, 'YTickLabel', []);
        set(gca, 'YColor', 'none');
        set(gca, 'FontSize', g.fontSize2);

        % CS+US panel
        ax_Both = nexttile(t_nested_Both, c);
        hold on;
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        if ~isempty(clust_idx)
            psth_Both_mean = mean(res.psth_Both_Hz(clust_idx, plot_idx), 1);
            plot(time_vec, psth_Both_mean, '-', 'Color', cluster_colors(c, :), 'LineWidth', 2.5);
        end

        hold off;
        xlim([time_vec(1) time_vec(end)]);
        ylim([y_min*1.1 y_max*1.1]);

        if br == 2 && c == 3
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            set(gca, 'XColor', 'k');
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
            set(gca, 'XColor', 'none');
        end

        set(gca, 'YTickLabel', []);
        set(gca, 'YColor', 'none');
        set(gca, 'FontSize', g.fontSize2);

        % Add scalebar on bottom cluster (c=3) at right side of CS+US panel
        if c == 3
            curr_ylim = ylim;
            y_range = curr_ylim(2) - curr_ylim(1);
            scalebar_size = 20;
            x_pos = time_vec(end) - 0.3;
            y_pos = curr_ylim(2) - scalebar_size - 0.05*y_range;

            hold on;
            plot([x_pos x_pos], [y_pos y_pos+scalebar_size], 'k-', 'LineWidth', 3);
            text(x_pos+0.2, y_pos+scalebar_size/2, sprintf('%d Hz', scalebar_size), ...
                'FontSize', g.fontSize2-2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            hold off;
        end
    end
end

% Add colorbar - manually create with full control, positioned at bottom heatmap (Astria)
drawnow;  % Ensure all positions are updated
cb_width = 0.010;
cb_left = 0.025;
cb_bottom = 0.48;  % Aligned with Astria heatmap (rows 3-4 in 7-row layout)
cb_height = 0.15;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'left');
ylabel(cb_ax, 'Z-score', 'FontSize', g.fontSize2);
cb_ax.FontSize = g.fontSize2;

%% Add metric bar charts (3 rows × 3 columns in bottom section)
% Rows 5-7: CS-sel, US-sel, Multisensory
% Columns span width of figure, divided into 3 metrics

% Create nested tiledlayout for all metrics
t_metrics = tiledlayout(t, 3, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
t_metrics.Layout.Tile = 9;  % Start at row 5, column 1
t_metrics.Layout.TileSpan = [3 2];  % Span 3 rows × 2 columns (full width)

metric_names = {'Response AUC', 'Response Peak', 'Peak Time'};

for c = [1 2 3]  % CS-selective, US-selective, Multisensory
    for metric = 1:3  % AUC, Peak, Time
        ax_metric = nexttile(t_metrics, (c-1)*3 + metric);

        % Collect data for LA and Astria
        data_LA = [];
        data_Astria = [];

        for br = 1:2
            if isempty(results_all{br})
                continue;
            end

            res = results_all{br};
            clust_idx = find(res.Clusters == c);

            fprintf('Cluster %d, Metric %d, Region %s: %d neurons\n', c, metric, brain_regions{br}, length(clust_idx));

            if ~isempty(clust_idx)
                if metric == 1  % Response AUC
                    CS_metric = sum(max(0, res.psth_CS(clust_idx, g.roi)), 2);
                    US_metric = sum(max(0, res.psth_US(clust_idx, g.roi)), 2);
                    Both_metric = sum(max(0, res.psth_Both(clust_idx, g.roi)), 2);
                elseif metric == 2  % Response Peak
                    CS_metric = max(res.psth_CS(clust_idx, g.roi), [], 2);
                    US_metric = max(res.psth_US(clust_idx, g.roi), [], 2);
                    Both_metric = max(res.psth_Both(clust_idx, g.roi), [], 2);
                else  % Peak Time (latency in ms)
                    CS_metric = nan(length(clust_idx), 1);
                    US_metric = nan(length(clust_idx), 1);
                    Both_metric = nan(length(clust_idx), 1);
                    for n = 1:length(clust_idx)
                        [~, idx] = max(res.psth_CS(clust_idx(n), g.roi));
                        CS_metric(n) = (idx - 1) * g.bin_time * 1000;  % Convert to ms
                        [~, idx] = max(res.psth_US(clust_idx(n), g.roi));
                        US_metric(n) = (idx - 1) * g.bin_time * 1000;
                        [~, idx] = max(res.psth_Both(clust_idx(n), g.roi));
                        Both_metric(n) = (idx - 1) * g.bin_time * 1000;
                    end
                end

                if br == 1
                    data_LA = [CS_metric, US_metric, Both_metric];
                else
                    data_Astria = [CS_metric, US_metric, Both_metric];
                end
            end
        end

        % Plot grouped bars - LA bars [2 3 4], Astria bars [7 8 9] with larger gap
        hold on;

        if ~isempty(data_LA) && ~isempty(data_Astria)
            means_LA = mean(data_LA, 1, 'omitnan');
            means_Astria = mean(data_Astria, 1, 'omitnan');
            sems_LA = std(data_LA, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_LA(:,1))));
            sems_Astria = std(data_Astria, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_Astria(:,1))));

            % LA bars at positions 2, 3, 4 (narrower bars) - use cluster color for dark version
            bar_color_LA = cluster_colors(c, :) * 0.7;  % Darker version
            bar([2 3 4], means_LA, 0.6, 'FaceColor', bar_color_LA, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([2 3 4], means_LA, sems_LA, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

            % Astria bars at positions 7, 8, 9 (gap at positions 5, 6) - use cluster color for light version
            bar_color_Astria = cluster_colors(c, :) + (1 - cluster_colors(c, :)) * 0.5;  % Lighter version
            bar([7 8 9], means_Astria, 0.6, 'FaceColor', bar_color_Astria, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([7 8 9], means_Astria, sems_Astria, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

        elseif ~isempty(data_LA)
            means_LA = mean(data_LA, 1, 'omitnan');
            sems_LA = std(data_LA, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_LA(:,1))));
            bar_color_LA = cluster_colors(c, :) * 0.7;  % Darker version
            bar([2 3 4], means_LA, 0.6, 'FaceColor', bar_color_LA, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([2 3 4], means_LA, sems_LA, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
        elseif ~isempty(data_Astria)
            means_Astria = mean(data_Astria, 1, 'omitnan');
            sems_Astria = std(data_Astria, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_Astria(:,1))));
            bar_color_Astria = cluster_colors(c, :) + (1 - cluster_colors(c, :)) * 0.5;  % Lighter version
            bar([7 8 9], means_Astria, 0.6, 'FaceColor', bar_color_Astria, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([7 8 9], means_Astria, sems_Astria, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
        end

        hold off;

        % Formatting
        xlim([0.5 10.5]);
        xticks([2 3 4 7 8 9]);

        % Calculate y-ticks based on data range BEFORE extending ylim
        curr_ylim_data = ylim;
        y_ticks = round([curr_ylim_data(1), mean(curr_ylim_data), curr_ylim_data(2)], 1);

        % Add titles on top row and adjust ylim AFTER calculating ticks
        if c == 1
            title(metric_names{metric}, 'FontSize', g.fontSize1, 'FontWeight', 'bold');
            % Increase upper ylim to give space for labels
            y_range = curr_ylim_data(2) - curr_ylim_data(1);
            ylim([curr_ylim_data(1), curr_ylim_data(2) + 0.15*y_range]);
        end

        % X-tick labels only on bottom row
        if c == 3
            xticklabels({'CS', 'US', 'CS+US', 'CS', 'US', 'CS+US'});
        else
            set(gca, 'XTickLabel', []);
        end

        % Set y-ticks (based on data range, not extended range)
        yticks(y_ticks);

        % Add LA/Astria labels on top row
        if c == 1
            curr_ylim_extended = ylim;  % Get the extended ylim
            y_range = curr_ylim_extended(2) - curr_ylim_extended(1);
            % Lower position for columns 2 and 3
            if metric == 1
                y_pos = curr_ylim_extended(2) - 0.05*y_range;
            else  % metric 2 or 3
                y_pos = curr_ylim_extended(2) - 0.10*y_range;
            end
            text(3, y_pos, 'LA', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            text(8, y_pos, 'Astria', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end

        % Add y-axis label for metric type
        if metric == 1
            ylabel(sprintf('%s\nZ-score', cluster_names{c}), 'FontSize', g.fontSize2, 'FontWeight', 'bold');
        elseif metric == 2
            ylabel('Z-score', 'FontSize', g.fontSize2);
        else  % metric == 3
            ylabel('Time (ms)', 'FontSize', g.fontSize2);
        end

        set(gca, 'FontSize', g.fontSize2);
        box off;
    end
end

fprintf('\nDone.\n');

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
