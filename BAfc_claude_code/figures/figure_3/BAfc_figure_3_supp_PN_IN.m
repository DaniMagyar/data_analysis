% BAfc_figure_3_supp_PN_IN
% Supplementary figure for BAfc_figure_3
% Compare CS-selective, US-selective, and Multisensory neurons identified from CS and US
% Show their responses to CS+US (both) stimulus
% LA only, separated by cell type (PN vs IN)

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
g.fontSize1 = 10;
g.fontSize2 = 10;
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

% LA only, separated by cell type
cell_types = {'PN', 'IN'};

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

% Calculate Savitzky-Golay filter delay correction
filter_delay = floor(g.smoothvalue / 2);  % Half the filter width (symmetric kernel)
fprintf('Savitzky-Golay filter width: %d bins, delay correction: %d bins (%.1f ms)\n', ...
    g.smoothvalue, filter_delay, filter_delay * g.bin_time * 1000);

for hmp = 1:numel(ttl)
    psth_spx_og = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Convert to Hz
    num_trials = size(g.cell_metrics.general.(ttl{hmp}){1}, 1);
    psth_hz = psth_spx_og / (num_trials * g.bin_time);
    psth_hz_smooth = smoothdata(psth_hz, 2, 'sgolay', g.smoothvalue);

    % Correct for filter delay by shifting forward (move data earlier in time)
    psth_hz_corrected = zeros(size(psth_hz_smooth));
    psth_hz_corrected(:, filter_delay+1:end) = psth_hz_smooth(:, 1:end-filter_delay);
    psth_hz_corrected(:, 1:filter_delay) = repmat(psth_hz_smooth(:, 1), 1, filter_delay);
    psthHz_full{hmp} = psth_hz_corrected;

    % Z-score
    baseline_idx = 1:(g.pre_time / g.bin_time);
    baseline_mean = mean(psth_spx_og(:, baseline_idx), 2);
    baseline_std = std(psth_spx_og(:, baseline_idx), 0, 2);
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;
    psth_spx_smooth = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);

    % Correct for filter delay by shifting forward (move data earlier in time)
    psth_spx_corrected = zeros(size(psth_spx_smooth));
    psth_spx_corrected(:, filter_delay+1:end) = psth_spx_smooth(:, 1:end-filter_delay);
    psth_spx_corrected(:, 1:filter_delay) = repmat(psth_spx_smooth(:, 1), 1, filter_delay);
    psthZ_full{hmp} = psth_spx_corrected;
end

%% Process each cell type (PN and IN)
results_all = cell(1, 2);

for ct = 1:2
    fprintf('\nProcessing LA %s...\n', cell_types{ct});

    % Get neuron indices for LA and specific cell type
    idx_neurons = strcmp(g.cell_metrics.brainRegion, 'LA') & strcmp(g.cell_metrics.putativeCellType, cell_types{ct});

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
    Both_onset_lat = nan(n_neurons, 1);
    Both_offset_lat = nan(n_neurons, 1);

    for n = 1:n_neurons
        [CS_onset_lat(n), CS_offset_lat(n)] = compute_onset_offset_latency(psth_CS(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [US_onset_lat(n), US_offset_lat(n)] = compute_onset_offset_latency(psth_US(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [Both_onset_lat(n), Both_offset_lat(n)] = compute_onset_offset_latency(psth_Both(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
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
    results_all{ct}.Clusters = Clusters;
    results_all{ct}.leafOrder = leafOrder;
    results_all{ct}.psth_CS = psth_CS;
    results_all{ct}.psth_US = psth_US;
    results_all{ct}.psth_Both = psth_Both;
    results_all{ct}.psth_CS_Hz = psth_CS_Hz;
    results_all{ct}.psth_US_Hz = psth_US_Hz;
    results_all{ct}.psth_Both_Hz = psth_Both_Hz;
    results_all{ct}.CS_onset_lat = CS_onset_lat;
    results_all{ct}.CS_offset_lat = CS_offset_lat;
    results_all{ct}.US_onset_lat = US_onset_lat;
    results_all{ct}.US_offset_lat = US_offset_lat;
    results_all{ct}.Both_onset_lat = Both_onset_lat;
    results_all{ct}.Both_offset_lat = Both_offset_lat;
    results_all{ct}.n_neurons = n_neurons;
end

%% Create figure
fig = figure('Units', 'pixels', 'Position', [100, 100, 1000, 1000], 'Visible', 'on');
t = tiledlayout(fig, 6, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add panel labels
% A: PN heatmaps (row 1, col 1)
annotation(fig, 'textbox', [0.01 0.95 0.05 0.05], 'String', 'A', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% B: PN lineplots (row 1, col 2)
annotation(fig, 'textbox', [0.50 0.95 0.05 0.05], 'String', 'B', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% C: IN heatmaps (row 2, col 1)
annotation(fig, 'textbox', [0.01 0.71 0.05 0.05], 'String', 'C', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% D: IN lineplots (row 2, col 2)
annotation(fig, 'textbox', [0.50 0.71 0.05 0.05], 'String', 'D', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% E: Bar graphs (row 4, col 1)
annotation(fig, 'textbox', [0.01 0.46 0.05 0.05], 'String', 'E', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Determine global color limits across all heatmaps (CS, US, and CS+US)
all_values = [];
for ct = 1:2
    if ~isempty(results_all{ct})
        for stim = 1:3
            if stim == 1
                psth_sorted = results_all{ct}.psth_CS(results_all{ct}.leafOrder, :);
            elseif stim == 2
                psth_sorted = results_all{ct}.psth_US(results_all{ct}.leafOrder, :);
            else
                psth_sorted = results_all{ct}.psth_Both(results_all{ct}.leafOrder, :);
            end
            matrix = psth_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
            all_values = [all_values; matrix(:)];
        end
    end
end

clim_min = prctile(all_values, (100 - g.clim_percentile) / 2);
clim_max = prctile(all_values, 100 - (100 - g.clim_percentile) / 2);

%% Create single 2×6 tiledlayout for all heatmaps and lineplots
t_plots = tiledlayout(t, 2, 6, 'TileSpacing', 'compact', 'Padding', 'tight');
t_plots.Layout.Tile = 1;  % Start at row 1, col 1
t_plots.Layout.TileSpan = [3 2];  % Span rows 1-3, cols 1-2

% Column arrangement: CS_hmp, US_hmp, Both_hmp, CS_line, US_line, Both_line
% Row 1: PN
% Row 2: IN

%% Plot each cell type (rows) × heatmaps+lineplots (columns)
for ct = 1:2
    if isempty(results_all{ct})
        continue;
    end

    res = results_all{ct};

    % Get cluster boundaries for responsive clusters only (1, 2, 3)
    Clusters_sorted = res.Clusters(res.leafOrder);

    % Filter to only plot clusters 1, 2, 3
    responsive_mask = ismember(Clusters_sorted, [1, 2, 3]);
    responsive_leafOrder = res.leafOrder(responsive_mask);
    Clusters_sorted_responsive = Clusters_sorted(responsive_mask);
    n_clu = find(diff(Clusters_sorted_responsive) ~= 0);

    % Plot heatmaps (columns 1-3)
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

        % Heatmap - tile position: ct=row, stim=column
        tile_idx = (ct-1)*6 + stim;
        ax = nexttile(t_plots, tile_idx);
        matrix = psth_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        imagesc(g.timeaxis_hmp, 1:size(matrix, 1), matrix);
        clim([clim_min clim_max]);
        colormap(ax, g.colors.Heatmap);
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        % Labels
        if ct == 1
            title(stim_title, 'FontSize', 10);
        end

        % Y-label only on first column with cell type
        if stim == 1
            ylabel(sprintf('%s (neuron #)', cell_types{ct}), 'FontSize', 10);
        end

        % Set yticks to first and last
        n_neurons_plot = size(matrix, 1);
        yticks([1, n_neurons_plot]);

        % Remove yticklabels from all but first column
        if stim ~= 1
            set(gca, 'YTickLabel', []);
        end

        if ct == 2
            xlabel('Time (s)', 'FontSize', 10);
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
        end

        % Set axis font size
        set(gca, 'FontSize', 10);

        % Add cluster lines (black)
        hold on;
        for i = 1:length(n_clu)
            yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
        end

        % Add onset and offset markers for each neuron
        for n = 1:length(responsive_leafOrder)
            idx_n = responsive_leafOrder(n);

            % Get appropriate onset/offset latencies based on stimulus
            if stim == 1  % CS
                onset_lat = res.CS_onset_lat(idx_n);
                offset_lat = res.CS_offset_lat(idx_n);
            elseif stim == 2  % US
                onset_lat = res.US_onset_lat(idx_n);
                offset_lat = res.US_offset_lat(idx_n);
            else  % CS+US
                onset_lat = res.Both_onset_lat(idx_n);
                offset_lat = res.Both_offset_lat(idx_n);
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
    end

    % Plot lineplots (columns 4-6) - stacked by cluster within each stimulus
    plot_idx = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    time_vec = g.timeaxis_hmp;

    % Ensure consistent lengths
    if length(time_vec) ~= length(plot_idx)
        time_vec = linspace(-g.plotwin(1), g.plotwin(2), length(plot_idx));
    end

    % Calculate global y-limits for this cell type (only clusters 1-3)
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

    fprintf('  %s lineplots: y_min = %.2f, y_max = %.2f\n', cell_types{ct}, y_min*1.1, y_max*1.1);

    % Create nested tiledlayouts for each stimulus column (4, 5, 6)
    for stim = 1:3  % CS, US, CS+US
        % Create nested 3×1 layout for this stimulus
        tile_idx = (ct-1)*6 + stim + 3;
        t_nested = tiledlayout(t_plots, 3, 1, 'TileSpacing', 'none', 'Padding', 'tight');
        t_nested.Layout.Tile = tile_idx;

        % Plot each cluster in separate row
        for c = [1 2 3]
            clust_idx = find(res.Clusters == c);

            ax_line = nexttile(t_nested, c);

            % Add title on first cluster of top row
            if ct == 1 && c == 1
                if stim == 1
                    title('CS', 'FontSize', 10, 'FontWeight', 'bold');
                elseif stim == 2
                    title('US', 'FontSize', 10, 'FontWeight', 'bold');
                else
                    title('CS+US', 'FontSize', 10, 'FontWeight', 'bold');
                end
            end

            hold on;
            xline(0, '--k', 'LineWidth', g.xlinewidth);

            if ~isempty(clust_idx)
                if stim == 1
                    psth_mean = mean(res.psth_CS_Hz(clust_idx, plot_idx), 1);
                elseif stim == 2
                    psth_mean = mean(res.psth_US_Hz(clust_idx, plot_idx), 1);
                else
                    psth_mean = mean(res.psth_Both_Hz(clust_idx, plot_idx), 1);
                end
                plot(time_vec, psth_mean, '-', 'Color', cluster_colors(c, :), 'LineWidth', 2.5);
            end

            hold off;
            xlim([time_vec(1) time_vec(end)]);
            ylim([y_min*1.1 y_max*1.1]);

            if ct == 2 && c == 3
                xlabel('Time (s)', 'FontSize', 10);
                set(gca, 'XColor', 'k');
                xticks([-1 0 1]);
            else
                set(gca, 'XTickLabel', []);
                set(gca, 'XColor', 'none');
            end

            set(gca, 'YTickLabel', []);
            set(gca, 'YColor', 'none');
            set(gca, 'FontSize', 10);

            % Add scalebar on bottom cluster (c=3) of CS+US panel
            if c == 3 && stim == 3
                curr_ylim = ylim;
                y_range = curr_ylim(2) - curr_ylim(1);
                scalebar_size = 20;
                x_pos = time_vec(end) - 0.3;
                y_pos = curr_ylim(2) - scalebar_size - 0.05*y_range;

                hold on;
                plot([x_pos x_pos], [y_pos y_pos+scalebar_size], 'k-', 'LineWidth', 3);
                text(x_pos+0.2, y_pos+scalebar_size/2, sprintf('%d Hz', scalebar_size), ...
                    'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
                hold off;
            end
        end
    end
end

% Add colorbar - manually create with full control, positioned at bottom heatmap (IN)
drawnow;  % Ensure all positions are updated
cb_width = 0.008;
cb_left = 0.50;
cb_bottom = 0.65;  % Aligned with IN heatmap (rows 2-3 in 6-row layout)
cb_height = 0.08;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'right');
%ylabel(cb_ax, 'Z-score', 'FontSize', 10);
cb_ax.FontSize = 10;

%% Add metric bar charts (2 rows × 3 columns in bottom section)
% Rows 4-5: US-sel, Multisensory (CS-sel removed)
% Columns span width of figure, divided into 3 metrics

% Create nested tiledlayout for all metrics (3 columns now)
t_metrics = tiledlayout(t, 2, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
t_metrics.Layout.Tile = 7;  % Start at row 4, column 1
t_metrics.Layout.TileSpan = [3 2];  % Span 3 rows × 2 columns

metric_names = {'ΔnSpikes', 'ΔFR (Hz)', 'Response length (ms)'};

for c = [2 3]  % US-selective, Multisensory only (CS-selective removed)
    for metric = 1:3  % ΔnSpikes, ΔFR, Response length
        ax_metric = nexttile(t_metrics, (c-2)*3 + metric);

        % Collect data for PN and IN
        data_PN = [];
        data_IN = [];

        for ct = 1:2
            if isempty(results_all{ct})
                continue;
            end

            res = results_all{ct};
            clust_idx = find(res.Clusters == c);

            fprintf('Cluster %d, Metric %d, CellType %s: %d neurons\n', c, metric, cell_types{ct}, length(clust_idx));

            if ~isempty(clust_idx)
                if metric == 1  % ΔnSpikes between onset and offset (change from baseline)
                    CS_metric = zeros(length(clust_idx), 1);
                    US_metric = zeros(length(clust_idx), 1);
                    Both_metric = zeros(length(clust_idx), 1);

                    baseline_idx = 1:(g.pre_time / g.bin_time);

                    for n = 1:length(clust_idx)
                        idx_n = clust_idx(n);

                        % Calculate baseline spike count per bin
                        baseline_spikes_per_bin = mean(res.psth_CS_Hz(idx_n, baseline_idx)) * g.bin_time;

                        % CS ΔnSpikes
                        if ~isnan(res.CS_onset_lat(idx_n)) && ~isnan(res.CS_offset_lat(idx_n))
                            onset_bin = round(res.CS_onset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            offset_bin = round(res.CS_offset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            num_bins = offset_bin - onset_bin + 1;
                            % Sum spikes in response window
                            response_spikes = sum(res.psth_CS_Hz(idx_n, onset_bin:offset_bin)) * g.bin_time;
                            % Subtract expected baseline spikes
                            CS_metric(n) = response_spikes - (baseline_spikes_per_bin * num_bins);
                        else
                            CS_metric(n) = 0;  % No response detected
                        end

                        % US ΔnSpikes
                        if ~isnan(res.US_onset_lat(idx_n)) && ~isnan(res.US_offset_lat(idx_n))
                            onset_bin = round(res.US_onset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            offset_bin = round(res.US_offset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            num_bins = offset_bin - onset_bin + 1;
                            response_spikes = sum(res.psth_US_Hz(idx_n, onset_bin:offset_bin)) * g.bin_time;
                            US_metric(n) = response_spikes - (baseline_spikes_per_bin * num_bins);
                        else
                            US_metric(n) = 0;
                        end

                        % CS+US ΔnSpikes
                        if ~isnan(res.Both_onset_lat(idx_n)) && ~isnan(res.Both_offset_lat(idx_n))
                            onset_bin = round(res.Both_onset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            offset_bin = round(res.Both_offset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            num_bins = offset_bin - onset_bin + 1;
                            response_spikes = sum(res.psth_Both_Hz(idx_n, onset_bin:offset_bin)) * g.bin_time;
                            Both_metric(n) = response_spikes - (baseline_spikes_per_bin * num_bins);
                        else
                            Both_metric(n) = 0;
                        end
                    end
                elseif metric == 2  % Firing rate change (Hz)
                    CS_metric = zeros(length(clust_idx), 1);
                    US_metric = zeros(length(clust_idx), 1);
                    Both_metric = zeros(length(clust_idx), 1);

                    baseline_idx = 1:(g.pre_time / g.bin_time);

                    for n = 1:length(clust_idx)
                        idx_n = clust_idx(n);

                        % Calculate baseline firing rate for this neuron
                        baseline_fr = mean(res.psth_CS_Hz(idx_n, baseline_idx));

                        % CS peak firing rate change
                        if ~isnan(res.CS_onset_lat(idx_n)) && ~isnan(res.CS_offset_lat(idx_n))
                            onset_bin = round(res.CS_onset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            offset_bin = round(res.CS_offset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            peak_fr = max(res.psth_CS_Hz(idx_n, onset_bin:offset_bin));
                            CS_metric(n) = peak_fr - baseline_fr;
                        else
                            CS_metric(n) = 0;  % No response detected
                        end

                        % US peak firing rate change
                        if ~isnan(res.US_onset_lat(idx_n)) && ~isnan(res.US_offset_lat(idx_n))
                            onset_bin = round(res.US_onset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            offset_bin = round(res.US_offset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            peak_fr = max(res.psth_US_Hz(idx_n, onset_bin:offset_bin));
                            US_metric(n) = peak_fr - baseline_fr;
                        else
                            US_metric(n) = 0;  % No response detected
                        end

                        % CS+US peak firing rate change
                        if ~isnan(res.Both_onset_lat(idx_n)) && ~isnan(res.Both_offset_lat(idx_n))
                            onset_bin = round(res.Both_onset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            offset_bin = round(res.Both_offset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                            peak_fr = max(res.psth_Both_Hz(idx_n, onset_bin:offset_bin));
                            Both_metric(n) = peak_fr - baseline_fr;
                        else
                            Both_metric(n) = 0;  % No response detected
                        end
                    end
                elseif metric == 3  % Response length (time between onset and offset in ms)
                    CS_metric = zeros(length(clust_idx), 1);
                    US_metric = zeros(length(clust_idx), 1);
                    Both_metric = zeros(length(clust_idx), 1);

                    for n = 1:length(clust_idx)
                        idx_n = clust_idx(n);

                        % CS response length
                        if ~isnan(res.CS_onset_lat(idx_n)) && ~isnan(res.CS_offset_lat(idx_n))
                            CS_metric(n) = (res.CS_offset_lat(idx_n) - res.CS_onset_lat(idx_n)) * 1000;  % Convert to ms
                        else
                            CS_metric(n) = 0;  % No response detected
                        end

                        % US response length
                        if ~isnan(res.US_onset_lat(idx_n)) && ~isnan(res.US_offset_lat(idx_n))
                            US_metric(n) = (res.US_offset_lat(idx_n) - res.US_onset_lat(idx_n)) * 1000;
                        else
                            US_metric(n) = 0;
                        end

                        % CS+US response length
                        if ~isnan(res.Both_onset_lat(idx_n)) && ~isnan(res.Both_offset_lat(idx_n))
                            Both_metric(n) = (res.Both_offset_lat(idx_n) - res.Both_onset_lat(idx_n)) * 1000;
                        else
                            Both_metric(n) = 0;
                        end
                    end
                end

                if ct == 1
                    data_PN = [CS_metric, US_metric, Both_metric];
                else
                    data_IN = [CS_metric, US_metric, Both_metric];
                end
            end
        end

        % Plot grouped bars - PN bars [2 3 4], IN bars [7 8 9] with larger gap
        hold on;

        if ~isempty(data_PN) && ~isempty(data_IN)
            means_PN = mean(data_PN, 1, 'omitnan');
            means_IN = mean(data_IN, 1, 'omitnan');
            sems_PN = std(data_PN, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_PN(:,1))));
            sems_IN = std(data_IN, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_IN(:,1))));

            % PN bars at positions 2, 3, 4 (narrower bars) - use cluster color for dark version
            bar_color_PN = cluster_colors(c, :) * 0.7;  % Darker version
            bar([2 3 4], means_PN, 0.6, 'FaceColor', bar_color_PN, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([2 3 4], means_PN, sems_PN, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

            % IN bars at positions 7, 8, 9 (gap at positions 5, 6) - use cluster color for light version
            bar_color_IN = cluster_colors(c, :) + (1 - cluster_colors(c, :)) * 0.5;  % Lighter version
            bar([7 8 9], means_IN, 0.6, 'FaceColor', bar_color_IN, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([7 8 9], means_IN, sems_IN, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

        elseif ~isempty(data_PN)
            means_PN = mean(data_PN, 1, 'omitnan');
            sems_PN = std(data_PN, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_PN(:,1))));
            bar_color_PN = cluster_colors(c, :) * 0.7;  % Darker version
            bar([2 3 4], means_PN, 0.6, 'FaceColor', bar_color_PN, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([2 3 4], means_PN, sems_PN, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
        elseif ~isempty(data_IN)
            means_IN = mean(data_IN, 1, 'omitnan');
            sems_IN = std(data_IN, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_IN(:,1))));
            bar_color_IN = cluster_colors(c, :) + (1 - cluster_colors(c, :)) * 0.5;  % Lighter version
            bar([7 8 9], means_IN, 0.6, 'FaceColor', bar_color_IN, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([7 8 9], means_IN, sems_IN, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
        end

        hold off;

        % Formatting - set ylim FIRST before calculating significance positions
        xlim([0.5 10.5]);
        xticks([2 3 4 7 8 9]);

        % Set fixed y-limits based on metric type
        if metric == 1  % ΔnSpikes
            ylim([0 100]);
            y_ticks = [0 50 100];
        elseif metric == 2  % ΔFR
            ylim([0 300]);
            y_ticks = [0 150 300];
        else  % metric == 3, Response length
            ylim([0 700]);
            y_ticks = [0 350 700];
        end

        % Set y-ticks
        yticks(y_ticks);

        % Statistical comparisons - Wilcoxon signed rank test (paired data)
        % Compare CS vs US, CS vs CS+US, US vs CS+US within each cell type

        % Set fixed spacing for significance levels based on metric type
        if metric == 1  % ΔnSpikes
            sig_spacing = 4;  % Fixed absolute spacing between significance levels
        elseif metric == 2  % ΔFR
            sig_spacing = 12;  % Fixed absolute spacing
        else  % metric == 3, Response length
            sig_spacing = 50;  % Fixed absolute spacing
        end

        % Get current ylim for consistent spacing calculations (AFTER ylim is set)
        curr_ylim = ylim;
        y_range = curr_ylim(2) - curr_ylim(1);

        % PN comparisons
        if ~isempty(data_PN)
            % CS vs US
            [p_PN_CS_US, ~] = signrank(data_PN(:,1), data_PN(:,2));
            % CS vs CS+US
            [p_PN_CS_Both, ~] = signrank(data_PN(:,1), data_PN(:,3));
            % US vs CS+US
            [p_PN_US_Both, ~] = signrank(data_PN(:,2), data_PN(:,3));

            % Add significance markers for PN
            y_max_PN = max(means_PN + sems_PN);
            hold on;

            % Level 1: Short comparisons (CS vs US, US vs CS+US)
            y_pos_level1 = y_max_PN + 0.08 * y_range;

            % CS vs US
            if p_PN_CS_US < 0.001
                sig_text = '***';
            elseif p_PN_CS_US < 0.01
                sig_text = '**';
            elseif p_PN_CS_US < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([2.1 2.9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                text(2.5, y_pos_level1 + 0.01 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10);
            end

            % US vs CS+US
            if p_PN_US_Both < 0.001
                sig_text = '***';
            elseif p_PN_US_Both < 0.01
                sig_text = '**';
            elseif p_PN_US_Both < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([3.1 3.9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                text(3.5, y_pos_level1 + 0.01 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10);
            end

            % Level 2: Long comparison (CS vs CS+US) - fixed absolute spacing above level 1
            y_pos_level2 = y_pos_level1 + sig_spacing;

            % CS vs CS+US
            if p_PN_CS_Both < 0.001
                sig_text = '***';
            elseif p_PN_CS_Both < 0.01
                sig_text = '**';
            elseif p_PN_CS_Both < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([2 4], [y_pos_level2 y_pos_level2], 'k-', 'LineWidth', 1.5);
                text(3, y_pos_level2 + 0.01 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10);
            end

            hold off;
        end

        % IN comparisons
        if ~isempty(data_IN)
            % CS vs US
            [p_IN_CS_US, ~] = signrank(data_IN(:,1), data_IN(:,2));
            % CS vs CS+US
            [p_IN_CS_Both, ~] = signrank(data_IN(:,1), data_IN(:,3));
            % US vs CS+US
            [p_IN_US_Both, ~] = signrank(data_IN(:,2), data_IN(:,3));

            % Add significance markers for IN
            y_max_IN = max(means_IN + sems_IN);
            hold on;

            % Level 1: Short comparisons (CS vs US, US vs CS+US)
            y_pos_level1 = y_max_IN + 0.08 * y_range;

            % CS vs US
            if p_IN_CS_US < 0.001
                sig_text = '***';
            elseif p_IN_CS_US < 0.01
                sig_text = '**';
            elseif p_IN_CS_US < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([7.1 7.9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                text(7.5, y_pos_level1 + 0.01 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10);
            end

            % US vs CS+US
            if p_IN_US_Both < 0.001
                sig_text = '***';
            elseif p_IN_US_Both < 0.01
                sig_text = '**';
            elseif p_IN_US_Both < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([8.1 8.9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                text(8.5, y_pos_level1 + 0.01 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10);
            end

            % Level 2: Long comparison (CS vs CS+US) - fixed absolute spacing above level 1
            y_pos_level2 = y_pos_level1 + sig_spacing;

            % CS vs CS+US
            if p_IN_CS_Both < 0.001
                sig_text = '***';
            elseif p_IN_CS_Both < 0.01
                sig_text = '**';
            elseif p_IN_CS_Both < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([7 9], [y_pos_level2 y_pos_level2], 'k-', 'LineWidth', 1.5);
                text(8, y_pos_level2 + 0.01 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10);
            end

            hold off;
        end

        % Add title on top row
        if c == 2
            title(metric_names{metric}, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end

        % X-tick labels only on bottom row
        if c == 3
            xticklabels({'CS', 'US', 'CS+US', 'CS', 'US', 'CS+US'});
        else
            set(gca, 'XTickLabel', []);
        end

        % Add PN/IN labels on top row
        if c == 2
            curr_ylim = ylim;
            y_pos = curr_ylim(2) * 0.95;  % Position at 95% of y-max
            text(3, y_pos, 'PN', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            text(8, y_pos, 'IN', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end

        % Add cluster name label only on first column
        if metric == 1
            ylabel(cluster_names{c}, 'FontSize', 10, 'FontWeight', 'bold');
        end

        set(gca, 'FontSize', 10);
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
