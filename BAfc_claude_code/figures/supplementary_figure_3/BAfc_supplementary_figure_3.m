% Supplementary figure 3
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

    % Calculate peak FR for CS, US, and CS+US (same as figure 3)
    CS_peak_Hz = max(psth_CS_Hz(:, g.roi), [], 2);
    US_peak_Hz = max(psth_US_Hz(:, g.roi), [], 2);
    Both_peak_Hz = max(psth_Both_Hz(:, g.roi), [], 2);

    % Store results
    results_all{ct}.Clusters = Clusters;
    results_all{ct}.leafOrder = leafOrder;
    results_all{ct}.psth_CS = psth_CS;
    results_all{ct}.psth_US = psth_US;
    results_all{ct}.psth_Both = psth_Both;
    results_all{ct}.psth_CS_Hz = psth_CS_Hz;
    results_all{ct}.psth_US_Hz = psth_US_Hz;
    results_all{ct}.psth_Both_Hz = psth_Both_Hz;
    results_all{ct}.CS_peak_Hz = CS_peak_Hz;
    results_all{ct}.US_peak_Hz = US_peak_Hz;
    results_all{ct}.Both_peak_Hz = Both_peak_Hz;
    results_all{ct}.CS_onset_lat = CS_onset_lat;
    results_all{ct}.CS_offset_lat = CS_offset_lat;
    results_all{ct}.US_onset_lat = US_onset_lat;
    results_all{ct}.US_offset_lat = US_offset_lat;
    results_all{ct}.Both_onset_lat = Both_onset_lat;
    results_all{ct}.Both_offset_lat = Both_offset_lat;
    results_all{ct}.n_neurons = n_neurons;
    results_all{ct}.animals = g.cell_metrics.animal(idx_neurons);
end

%% Storage for Friedman test data
% Structure: kw_data_storage{cell_type, cluster, stimulus}
kw_data_storage = cell(2, 3, 3);  % 2 cell types × 3 clusters × 3 stimuli (CS, US, Both)

% Populate kw_data_storage with peak FR values
for ct = 1:2
    if isempty(results_all{ct})
        continue;
    end

    res = results_all{ct};

    for c = 1:3  % CS-selective, US-selective, Multisensory
        clust_idx = find(res.Clusters == c);

        if ~isempty(clust_idx)
            % Store peak FR values (same as figure 3)
            kw_data_storage{ct, c, 1} = res.CS_peak_Hz(clust_idx);
            kw_data_storage{ct, c, 2} = res.US_peak_Hz(clust_idx);
            kw_data_storage{ct, c, 3} = res.Both_peak_Hz(clust_idx);
        end
    end
end

%% Create figure (2×2 grid: heatmaps and lineplots only)
fig = figure('Units', 'pixels', 'Position', [100, 100, 1000, 600], 'Visible', 'on');
t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

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
annotation(fig, 'textbox', [0.01 0.48 0.05 0.05], 'String', 'C', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% D: IN lineplots (row 2, col 2)
annotation(fig, 'textbox', [0.50 0.48 0.05 0.05], 'String', 'D', ...
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
t_plots.Layout.TileSpan = [2 2];  % Span all rows, all cols

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
cb_bottom = 0.15;  % Aligned with IN heatmap (row 2 in 2-row layout)
cb_height = 0.15;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'right');
cb_ax.FontSize = 10;

%% Perform Friedman tests for Excel export
% Perform Friedman tests for each cell type × cluster combination
kw_results = struct();

cluster_names_kw = {'CS-sel', 'US-sel', 'Multi'};

for ct = 1:2
    if isempty(results_all{ct})
        continue;
    end

    for c = 1:3  % CS-selective, US-selective, Multisensory
        % Check if we have data for this cluster
        if isempty(kw_data_storage{ct, c, 1})
            continue;
        end

        % Prepare data matrix for Friedman test (repeated measures)
        % Rows = subjects (neurons), Columns = conditions (CS, US, Both)
        % Friedman test requires at least 2 neurons
        if length(kw_data_storage{ct, c, 1}) >= 2
            % Ensure data is properly formatted as column vectors
            data_matrix = [kw_data_storage{ct, c, 1}(:), kw_data_storage{ct, c, 2}(:), kw_data_storage{ct, c, 3}(:)];

            % Friedman test (repeated measures)
            [p_friedman, ~, ~] = friedman(data_matrix, 1, 'off');

            % Only perform post-hoc tests if Friedman test is significant
            if p_friedman < 0.05
                % Post-hoc pairwise Wilcoxon signed-rank tests
                p_values = zeros(3, 1);
                p_values(1) = signrank(kw_data_storage{ct, c, 1}, kw_data_storage{ct, c, 2});  % CS vs US
                p_values(2) = signrank(kw_data_storage{ct, c, 1}, kw_data_storage{ct, c, 3});  % CS vs Both
                p_values(3) = signrank(kw_data_storage{ct, c, 2}, kw_data_storage{ct, c, 3});  % US vs Both
            else
                % Set p-values to 1 (non-significant) if Friedman is not significant
                p_values = ones(3, 1);
            end
        else
            % Not enough data for Friedman test (n < 2)
            p_friedman = 1;
            p_values = ones(3, 1);
        end

        % Store results
        kw_results(ct, c).p_friedman = p_friedman;
        kw_results(ct, c).p_values = p_values;
    end
end

%% Export data to Excel
export_figure3_supp_PN_IN_to_excel(results_all, kw_data_storage, kw_results, g, cell_types, cluster_names);
exportgraphics(gcf, 'supplementary_figure_3.png', 'Resolution', 300);
fprintf('\nDone.\n');

%% Helper functions
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

function export_figure3_supp_PN_IN_to_excel(results_all, kw_data_storage, kw_results, g, cell_types, cluster_names)
    output_filename = 'supplementary_figure_3_data.xlsx';

    if exist(output_filename, 'file')
        delete(output_filename);
    end

    stim_names = {'CS', 'US', 'CS+US'};
    all_cluster_names = {'CS-sel', 'US-sel', 'Multi', 'Non-resp', 'Inhibited'};

    %% Sheet 1: Peak FR and Response Length (side by side layout like figure 3)
    sheet_data = {};
    sheet_data{1, 1} = 'LA PN vs IN: Peak FR and Response Length';
    sheet_data{2, 1} = '';
    row = 3;

    % Process each cluster (1=CS-sel, 2=US-sel, 3=Multi)
    for c = 1:3
        sheet_data{row, 1} = sprintf('--- %s ---', cluster_names{c});
        row = row + 1;

        % Process each cell type (PN and IN)
        for ct = 1:2
            if isempty(results_all{ct})
                continue;
            end

            res = results_all{ct};
            clust_idx = find(res.Clusters == c);

            if isempty(clust_idx)
                sheet_data{row, 1} = sprintf('%s: No neurons', cell_types{ct});
                row = row + 1;
                continue;
            end

            % Calculate response length for each stimulus
            CS_resplen = zeros(length(clust_idx), 1);
            US_resplen = zeros(length(clust_idx), 1);
            Both_resplen = zeros(length(clust_idx), 1);

            for n = 1:length(clust_idx)
                idx_n = clust_idx(n);

                % CS response length
                if ~isnan(res.CS_onset_lat(idx_n)) && ~isnan(res.CS_offset_lat(idx_n))
                    CS_resplen(n) = (res.CS_offset_lat(idx_n) - res.CS_onset_lat(idx_n)) * 1000;
                end

                % US response length
                if ~isnan(res.US_onset_lat(idx_n)) && ~isnan(res.US_offset_lat(idx_n))
                    US_resplen(n) = (res.US_offset_lat(idx_n) - res.US_onset_lat(idx_n)) * 1000;
                end

                % CS+US response length
                if ~isnan(res.Both_onset_lat(idx_n)) && ~isnan(res.Both_offset_lat(idx_n))
                    Both_resplen(n) = (res.Both_offset_lat(idx_n) - res.Both_onset_lat(idx_n)) * 1000;
                end
            end

            % Perform Friedman test for response length
            response_length_data = [CS_resplen, US_resplen, Both_resplen];
            if length(clust_idx) >= 2
                [p_friedman_rl, ~, ~] = friedman(response_length_data, 1, 'off');

                % Get post-hoc for response length
                if p_friedman_rl < 0.05
                    p_cs_us_rl = signrank(CS_resplen, US_resplen);
                    p_cs_both_rl = signrank(CS_resplen, Both_resplen);
                    p_us_both_rl = signrank(US_resplen, Both_resplen);
                else
                    p_cs_us_rl = 1;
                    p_cs_both_rl = 1;
                    p_us_both_rl = 1;
                end
            else
                p_friedman_rl = 1;
                p_cs_us_rl = 1;
                p_cs_both_rl = 1;
                p_us_both_rl = 1;
            end

            % Header with both Peak FR and Response Length side by side
            sheet_data{row, 1} = sprintf('%s (n=%d)', cell_types{ct}, length(clust_idx));
            row = row + 1;

            sheet_data{row, 1} = 'Stimulus';
            sheet_data{row, 2} = 'Peak FR Mean (Hz)';
            sheet_data{row, 3} = 'Peak FR SEM (Hz)';
            sheet_data{row, 4} = 'Peak FR Median (Hz)';
            sheet_data{row, 5} = 'Peak FR SD (Hz)';
            sheet_data{row, 6} = '';  % Empty column separator
            sheet_data{row, 7} = 'Resp Length Mean (ms)';
            sheet_data{row, 8} = 'Resp Length SEM (ms)';
            sheet_data{row, 9} = 'Resp Length Median (ms)';
            sheet_data{row, 10} = 'Resp Length SD (ms)';
            row = row + 1;

            % Data for each stimulus (side by side)
            stim_metrics_fr = {kw_data_storage{ct, c, 1}, kw_data_storage{ct, c, 2}, kw_data_storage{ct, c, 3}};
            resp_metrics = {CS_resplen, US_resplen, Both_resplen};

            for stim = 1:3
                data_fr = stim_metrics_fr{stim};
                data_rl = resp_metrics{stim};

                sheet_data{row, 1} = stim_names{stim};
                sheet_data{row, 2} = mean(data_fr, 'omitnan');
                sheet_data{row, 3} = std(data_fr, 'omitnan') / sqrt(sum(~isnan(data_fr)));
                sheet_data{row, 4} = median(data_fr, 'omitnan');
                sheet_data{row, 5} = std(data_fr, 'omitnan');
                sheet_data{row, 6} = '';  % Empty column separator
                sheet_data{row, 7} = mean(data_rl, 'omitnan');
                sheet_data{row, 8} = std(data_rl, 'omitnan') / sqrt(sum(~isnan(data_rl)));
                sheet_data{row, 9} = median(data_rl, 'omitnan');
                sheet_data{row, 10} = std(data_rl, 'omitnan');
                row = row + 1;
            end

            % Add statistical test results side by side
            sheet_data{row, 1} = '';
            row = row + 1;

            % Statistical test headers
            sheet_data{row, 1} = 'Statistical Tests:';
            sheet_data{row, 2} = 'Peak FR';
            sheet_data{row, 6} = '';  % Empty column separator
            sheet_data{row, 7} = 'Response Length';
            row = row + 1;

            % Friedman test results
            if ~isempty(kw_results) && numel(kw_results) >= (ct + (c-1)*2)
                p_friedman = kw_results(ct, c).p_friedman;
                p_values = kw_results(ct, c).p_values;

                sheet_data{row, 1} = 'Friedman test p-value:';
                sheet_data{row, 2} = p_friedman;
                sheet_data{row, 3} = format_significance(p_friedman);
                sheet_data{row, 6} = '';  % Empty column separator
                sheet_data{row, 7} = p_friedman_rl;
                sheet_data{row, 8} = format_significance(p_friedman_rl);
                row = row + 1;

                % Post-hoc header
                if p_friedman < 0.05 || p_friedman_rl < 0.05
                    sheet_data{row, 1} = 'Post-hoc (Wilcoxon signed-rank):';
                    row = row + 1;

                    % CS vs US
                    sheet_data{row, 1} = '  CS vs US:';
                    if p_friedman < 0.05
                        sheet_data{row, 2} = p_values(1);
                        sheet_data{row, 3} = format_significance(p_values(1));
                    else
                        sheet_data{row, 2} = 'n/a';
                        sheet_data{row, 3} = '';
                    end
                    sheet_data{row, 6} = '';  % Empty column separator
                    if p_friedman_rl < 0.05
                        sheet_data{row, 7} = p_cs_us_rl;
                        sheet_data{row, 8} = format_significance(p_cs_us_rl);
                    else
                        sheet_data{row, 7} = 'n/a';
                        sheet_data{row, 8} = '';
                    end
                    row = row + 1;

                    % CS vs CS+US
                    sheet_data{row, 1} = '  CS vs CS+US:';
                    if p_friedman < 0.05
                        sheet_data{row, 2} = p_values(2);
                        sheet_data{row, 3} = format_significance(p_values(2));
                    else
                        sheet_data{row, 2} = 'n/a';
                        sheet_data{row, 3} = '';
                    end
                    sheet_data{row, 6} = '';  % Empty column separator
                    if p_friedman_rl < 0.05
                        sheet_data{row, 7} = p_cs_both_rl;
                        sheet_data{row, 8} = format_significance(p_cs_both_rl);
                    else
                        sheet_data{row, 7} = 'n/a';
                        sheet_data{row, 8} = '';
                    end
                    row = row + 1;

                    % US vs CS+US
                    sheet_data{row, 1} = '  US vs CS+US:';
                    if p_friedman < 0.05
                        sheet_data{row, 2} = p_values(3);
                        sheet_data{row, 3} = format_significance(p_values(3));
                    else
                        sheet_data{row, 2} = 'n/a';
                        sheet_data{row, 3} = '';
                    end
                    sheet_data{row, 6} = '';  % Empty column separator
                    if p_friedman_rl < 0.05
                        sheet_data{row, 7} = p_us_both_rl;
                        sheet_data{row, 8} = format_significance(p_us_both_rl);
                    else
                        sheet_data{row, 7} = 'n/a';
                        sheet_data{row, 8} = '';
                    end
                    row = row + 1;
                else
                    sheet_data{row, 1} = '  (Post-hoc not performed - Both Friedman tests n.s.)';
                    row = row + 1;
                end
            end

            sheet_data{row, 1} = '';
            row = row + 1;
        end

        sheet_data{row, 1} = '';
        row = row + 1;
    end

    writecell(sheet_data, output_filename, 'Sheet', 'Summary_Stats');

    %% Sheet 2: Raw Data - Individual Neuron Values
    raw_data = {};
    raw_data{1, 1} = 'LA PN vs IN: Individual Neuron Peak FR and Response Length Values';
    raw_data{2, 1} = 'These are the individual data points used to calculate means and SEMs';
    raw_data{3, 1} = 'Local #';
    raw_data{3, 2} = 'Global Index';
    raw_data{3, 3} = 'Animal ID';
    raw_data{3, 4} = 'Cell Type';
    raw_data{3, 5} = 'Cluster';
    raw_data{3, 6} = '';
    raw_data{3, 7} = 'Peak FR CS (Hz)';
    raw_data{3, 8} = 'Peak FR US (Hz)';
    raw_data{3, 9} = 'Peak FR CS+US (Hz)';
    raw_data{3, 10} = '';
    raw_data{3, 11} = 'Response Length CS (ms)';
    raw_data{3, 12} = 'Response Length US (ms)';
    raw_data{3, 13} = 'Response Length CS+US (ms)';

    row = 4;

    % Get cell_metrics for animal IDs
    cell_metrics = g.cell_metrics;

    % Process each cell type
    for ct = 1:2
        if isempty(results_all{ct})
            continue;
        end

        res = results_all{ct};

        % Get indices for this cell type in the full cell_metrics
        idx_neurons = strcmp(cell_metrics.brainRegion, 'LA') & strcmp(cell_metrics.putativeCellType, cell_types{ct});
        global_indices = find(idx_neurons);

        baseline_idx = 1:(g.pre_time / g.bin_time);

        % Process each neuron in sorted order (using leafOrder)
        for i = 1:length(res.leafOrder)
            n = res.leafOrder(i);  % Get neuron index in sorted order
            global_idx = global_indices(n);
            animal_id = cell_metrics.animal{global_idx};

            % Calculate metrics
            baseline_fr = mean(res.psth_CS_Hz(n, baseline_idx));

            % Delta FR
            CS_fr = 0;
            US_fr = 0;
            Both_fr = 0;
            CS_resplen = 0;
            US_resplen = 0;
            Both_resplen = 0;

            if ~isnan(res.CS_onset_lat(n)) && ~isnan(res.CS_offset_lat(n))
                onset_bin = round(res.CS_onset_lat(n) / g.bin_time) + 1 + g.roi(1) - 1;
                offset_bin = round(res.CS_offset_lat(n) / g.bin_time) + 1 + g.roi(1) - 1;
                peak_fr = max(res.psth_CS_Hz(n, onset_bin:offset_bin));
                CS_fr = peak_fr - baseline_fr;
                CS_resplen = (res.CS_offset_lat(n) - res.CS_onset_lat(n)) * 1000;
            end

            if ~isnan(res.US_onset_lat(n)) && ~isnan(res.US_offset_lat(n))
                onset_bin = round(res.US_onset_lat(n) / g.bin_time) + 1 + g.roi(1) - 1;
                offset_bin = round(res.US_offset_lat(n) / g.bin_time) + 1 + g.roi(1) - 1;
                peak_fr = max(res.psth_US_Hz(n, onset_bin:offset_bin));
                US_fr = peak_fr - baseline_fr;
                US_resplen = (res.US_offset_lat(n) - res.US_onset_lat(n)) * 1000;
            end

            if ~isnan(res.Both_onset_lat(n)) && ~isnan(res.Both_offset_lat(n))
                onset_bin = round(res.Both_onset_lat(n) / g.bin_time) + 1 + g.roi(1) - 1;
                offset_bin = round(res.Both_offset_lat(n) / g.bin_time) + 1 + g.roi(1) - 1;
                peak_fr = max(res.psth_Both_Hz(n, onset_bin:offset_bin));
                Both_fr = peak_fr - baseline_fr;
                Both_resplen = (res.Both_offset_lat(n) - res.Both_onset_lat(n)) * 1000;
            end

            raw_data{row, 1} = i;  % Sorted order (1, 2, 3, ...)
            raw_data{row, 2} = global_idx;
            raw_data{row, 3} = animal_id;
            raw_data{row, 4} = cell_types{ct};
            raw_data{row, 5} = all_cluster_names{res.Clusters(n)};
            raw_data{row, 6} = '';
            % Use peak FR from kw_data_storage (if available for this cluster)
            cluster_id = res.Clusters(n);
            if cluster_id <= 3 && ~isempty(kw_data_storage{ct, cluster_id, 1})
                % Find this neuron's position in the cluster
                clust_neurons = find(res.Clusters == cluster_id);
                pos_in_cluster = find(clust_neurons == n);
                if ~isempty(pos_in_cluster)
                    raw_data{row, 7} = kw_data_storage{ct, cluster_id, 1}(pos_in_cluster);
                    raw_data{row, 8} = kw_data_storage{ct, cluster_id, 2}(pos_in_cluster);
                    raw_data{row, 9} = kw_data_storage{ct, cluster_id, 3}(pos_in_cluster);
                else
                    raw_data{row, 7} = CS_fr;
                    raw_data{row, 8} = US_fr;
                    raw_data{row, 9} = Both_fr;
                end
            else
                raw_data{row, 7} = CS_fr;
                raw_data{row, 8} = US_fr;
                raw_data{row, 9} = Both_fr;
            end
            raw_data{row, 10} = '';
            raw_data{row, 11} = CS_resplen;
            raw_data{row, 12} = US_resplen;
            raw_data{row, 13} = Both_resplen;

            row = row + 1;
        end
    end

    writecell(raw_data, output_filename, 'Sheet', 'RawData_PeakFR_RespLen');
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
