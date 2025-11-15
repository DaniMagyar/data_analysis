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

% Brain regions and cell types to analyze
brain_regions = {'LA', 'Astria'};
cell_type_filter = {'all', 'all'};

cluster_colors = [
    0.8 0.2 0.2;    % CS-selective
    0.2 0.4 0.8;    % US-selective
    0.6 0.2 0.6;    % Multisensory
];

cluster_names = {'CS-sel', 'US-sel', 'CS&US'};

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
    results_all{br}.Clusters = Clusters;
    results_all{br}.leafOrder = leafOrder;
    results_all{br}.psth_CS = psth_CS;
    results_all{br}.psth_US = psth_US;
    results_all{br}.psth_Both = psth_Both;
    results_all{br}.psth_CS_Hz = psth_CS_Hz;
    results_all{br}.psth_US_Hz = psth_US_Hz;
    results_all{br}.psth_Both_Hz = psth_Both_Hz;
    results_all{br}.CS_onset_lat = CS_onset_lat;
    results_all{br}.CS_offset_lat = CS_offset_lat;
    results_all{br}.US_onset_lat = US_onset_lat;
    results_all{br}.US_offset_lat = US_offset_lat;
    results_all{br}.Both_onset_lat = Both_onset_lat;
    results_all{br}.Both_offset_lat = Both_offset_lat;
    results_all{br}.n_neurons = n_neurons;
end

%% Create main figure (heatmaps, lineplots, and delta Peak FR bars)
fig = figure('Units', 'pixels', 'Position', [100, 100, 1000, 1000], 'Visible', 'on');
t = tiledlayout(fig, 5, 7, 'TileSpacing', 'compact', 'Padding', 'compact');

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

%% Create single 2×6 tiledlayout for all heatmaps and lineplots
t_plots = tiledlayout(t, 2, 6, 'TileSpacing', 'compact', 'Padding', 'tight');
t_plots.Layout.Tile = 1;  % Start at row 1, col 1
t_plots.Layout.TileSpan = [4 6];  % Span rows 1-4, cols 1-6

% Column arrangement: CS_hmp, US_hmp, Both_hmp, CS_line, US_line, Both_line
% Row 1: LA
% Row 2: Astria
% Column 7: Delta Peak FR bars (nested 3×1 for each region)

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

        % Heatmap - tile position: br=row, stim=column
        tile_idx = (br-1)*6 + stim;
        ax = nexttile(t_plots, tile_idx);
        matrix = psth_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        imagesc(g.timeaxis_hmp, 1:size(matrix, 1), matrix);
        clim([clim_min clim_max]);
        colormap(ax, g.colors.Heatmap);
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        % Labels
        if br == 1
            title(stim_title, 'FontSize', g.fontSize2, 'FontWeight', 'bold');
        end

        % Set yticks to first and last
        n_neurons_plot = size(matrix, 1);
        yticks([1, n_neurons_plot]);

        % Add brain region as ylabel (only on first column)
        if stim == 1
            % Display AStria instead of Astria
            if strcmp(brain_regions{br}, 'Astria')
                ylabel('AStria neurons', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
            else
                ylabel([brain_regions{br} ' neurons'], 'FontSize', g.fontSize2, 'FontWeight', 'bold');
            end
        else
            % Remove yticklabels for columns 2 and 3
            set(gca, 'YTickLabel', []);
        end

        if br == 2
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
        end

        % Set axis font size
        set(gca, 'FontSize', g.fontSize2);

        % Add cluster lines (black)
        hold on;
        for i = 1:length(n_clu)
            yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
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

    % Create nested tiledlayouts for each stimulus column (4, 5, 6)
    for stim = 1:3  % CS, US, CS+US
        % Create nested 3×1 layout for this stimulus
        tile_idx = (br-1)*6 + stim + 3;
        t_nested = tiledlayout(t_plots, 3, 1, 'TileSpacing', 'none', 'Padding', 'tight');
        t_nested.Layout.Tile = tile_idx;

        % Plot each cluster in separate row
        for c = [1 2 3]
            clust_idx = find(res.Clusters == c);

            ax_line = nexttile(t_nested, c);

            % Add title on first cluster of top row
            if br == 1 && c == 1
                if stim == 1
                    title('CS', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
                elseif stim == 2
                    title('US', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
                else
                    title('CS+US', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
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

            % Add scalebar on bottom cluster (c=3) of CS panel (first lineplot)
            if c == 3 && stim == 1
                curr_ylim = ylim;
                y_range = curr_ylim(2) - curr_ylim(1);
                scalebar_size = 20;
                x_pos = time_vec(end) - 1.0;  % Shifted more to the left
                y_pos = curr_ylim(1) + scalebar_size + 0.05*y_range;  % Positioned near bottom

                hold on;
                plot([x_pos x_pos], [y_pos y_pos+scalebar_size], 'k-', 'LineWidth', 3);
                text(x_pos+0.2, y_pos+scalebar_size/2, sprintf('%d Hz', scalebar_size), ...
                    'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
                hold off;
            end
        end
    end
end

%% Add Delta Peak FR bar charts in column 7
% Create 2×1 container for LA and Astria bar sections
t_bars_container = tiledlayout(t, 2, 1, 'TileSpacing', 'compact', 'Padding', 'tight');
t_bars_container.Layout.Tile = 7;  % Column 7
t_bars_container.Layout.TileSpan = [4 1];  % Span rows 1-4

for br = 1:2
    if isempty(results_all{br})
        continue;
    end

    res = results_all{br};

    % Create nested 3×1 layout for this region's bars (matching lineplot structure)
    t_bars_nested = tiledlayout(t_bars_container, 3, 1, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_bars_nested.Layout.Tile = br;  % LA=1, Astria=2

    % Plot each cluster in separate row
    for c = [1 2 3]  % CS-selective, US-selective, Multisensory
        ax_bar = nexttile(t_bars_nested, c);

        clust_idx = find(res.Clusters == c);

        if ~isempty(clust_idx)
            % Calculate Delta Peak FR (metric == 2 code from bar chart section)
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
                    CS_metric(n) = 0;
                end

                % US peak firing rate change
                if ~isnan(res.US_onset_lat(idx_n)) && ~isnan(res.US_offset_lat(idx_n))
                    onset_bin = round(res.US_onset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                    offset_bin = round(res.US_offset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                    peak_fr = max(res.psth_US_Hz(idx_n, onset_bin:offset_bin));
                    US_metric(n) = peak_fr - baseline_fr;
                else
                    US_metric(n) = 0;
                end

                % CS+US peak firing rate change
                if ~isnan(res.Both_onset_lat(idx_n)) && ~isnan(res.Both_offset_lat(idx_n))
                    onset_bin = round(res.Both_onset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                    offset_bin = round(res.Both_offset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                    peak_fr = max(res.psth_Both_Hz(idx_n, onset_bin:offset_bin));
                    Both_metric(n) = peak_fr - baseline_fr;
                else
                    Both_metric(n) = 0;
                end
            end

            % Calculate means and SEMs
            means_data = mean([CS_metric, US_metric, Both_metric], 1, 'omitnan');
            sems_data = std([CS_metric, US_metric, Both_metric], 0, 1, 'omitnan') ./ sqrt(sum(~isnan(CS_metric)));

            % Plot bars
            hold on;
            % Use same color as pie chart
            bar_color = cluster_colors(c, :);  % Match pie chart colors

            bar([1 2 3], means_data, 0.4, 'FaceColor', bar_color, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([1 2 3], means_data, sems_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

            % Statistical comparisons - Wilcoxon signed rank test
            [p_CS_US, ~] = signrank(CS_metric, US_metric);
            [p_CS_Both, ~] = signrank(CS_metric, Both_metric);
            [p_US_Both, ~] = signrank(US_metric, Both_metric);

            % Add significance markers
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
                text(1.5, y_pos_level1 + 0.05 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
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
                text(2.5, y_pos_level1 + 0.05 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            % Level 2: Long comparison (CS vs CS+US)
            y_pos_level2 = y_pos_level1 + 12;  % Fixed spacing for Delta Peak FR

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
                text(2, y_pos_level2 + 0.05 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            hold off;
        end

        % Formatting
        xlim([0.5 3.5]);
        xticks([1 2 3]);
        ylim([0 100]);
        yticks([0 50 100]);

        % Add title on top cluster of LA
        if br == 1 && c == 1
            title('ΔFR (Hz)', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end

        % X-tick labels only on bottom cluster of Astria (br==2)
        if c == 3 && br == 2
            xticklabels({'CS', 'US', 'CS+US'});
        else
            set(gca, 'XTickLabel', []);
        end

        set(gca, 'FontSize', g.fontSize2);
        box off;
    end
end

%% Add row 5 with two nested tiledlayouts (pie charts and bar graphs)
% Create main container for row 5
t_row5 = tiledlayout(t, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
t_row5.Layout.Tile = 29;  % Start at row 5, col 1 (in 5x7 grid: (5-1)*7+1=29)
t_row5.Layout.TileSpan = [1 7];  % Span row 5, cols 1-7

% Nested tiledlayout 1: Pie charts (1×2)
t_pie = tiledlayout(t_row5, 1, 2, 'TileSpacing', 'tight', 'Padding', 'none');
t_pie.Layout.Tile = 1;

% Calculate cluster proportions for each region
for br = 1:2
    if isempty(results_all{br})
        continue;
    end

    res = results_all{br};

    % Count neurons in each cluster (only clusters 1, 2, 3)
    n_CS_sel = sum(res.Clusters == 1);
    n_US_sel = sum(res.Clusters == 2);
    n_Multi = sum(res.Clusters == 3);

    % Create pie chart
    ax_pie = nexttile(t_pie, br);

    pie_data = [n_CS_sel, n_US_sel, n_Multi];
    total_n = n_CS_sel + n_US_sel + n_Multi;

    % Calculate percentages
    pct_CS = (n_CS_sel / total_n) * 100;
    pct_US = (n_US_sel / total_n) * 100;
    pct_Multi = (n_Multi / total_n) * 100;

    pie_labels = {sprintf('%.0f%%', pct_CS), ...
                  sprintf('%.0f%%', pct_US), ...
                  sprintf('%.0f%%', pct_Multi)};

    p = pie(ax_pie, pie_data, pie_labels);

    % Set colors
    for i = 1:2:length(p)  % Every other element is a patch
        patch_idx = (i+1)/2;
        set(p(i), 'FaceColor', cluster_colors(patch_idx, :));
        set(p(i), 'EdgeColor', 'k');
        set(p(i), 'LineWidth', 1);
    end

    % Set text properties and move inside pie - FontSize 12
    for i = 2:2:length(p)  % Text elements
        set(p(i), 'FontSize', 12);
        set(p(i), 'FontWeight', 'bold');
        set(p(i), 'Color', 'w');  % White font color
        % Move text closer to center (inside the pie)
        pos = get(p(i), 'Position');
        set(p(i), 'Position', pos * 0.15);  % Move toward center
    end

    % Add title
    if strcmp(brain_regions{br}, 'Astria')
        title('AStria', 'FontSize', 10, 'FontWeight', 'bold');
    else
        title(brain_regions{br}, 'FontSize', 10, 'FontWeight', 'bold');
    end

    % Add legend to the right of the first pie chart (LA)
    if br == 1
        % Create invisible plot objects with square markers for legend
        hold(ax_pie, 'on');
        h1 = plot(ax_pie, NaN, NaN, 's', 'MarkerSize', 12, 'MarkerFaceColor', cluster_colors(1, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        h2 = plot(ax_pie, NaN, NaN, 's', 'MarkerSize', 12, 'MarkerFaceColor', cluster_colors(2, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        h3 = plot(ax_pie, NaN, NaN, 's', 'MarkerSize', 12, 'MarkerFaceColor', cluster_colors(3, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        lgd = legend(ax_pie, [h1, h2, h3], {'CS-sel', 'US-sel', 'CS&US'}, ...
                     'Location', 'eastoutside', 'FontSize', 10, 'Box', 'off');
        lgd.ItemTokenSize = [30, 30];  % Increase marker size in legend
    end

end

% Add title to pie chart section
title(t_pie, 'Response categories', 'FontSize', 12, 'FontWeight', 'bold');

% Nested tiledlayout 2: Bar graphs (1×3)
t_bars = tiledlayout(t_row5, 1, 3, 'TileSpacing', 'none', 'Padding', 'tight');
t_bars.Layout.Tile = 2;

%% Add comparison bar plots (CS, US, CS+US LA vs Astria Delta Peak FR)
stim_names = {'CS', 'US', 'CS+US'};

for stim = 1:3  % CS, US, CS+US
    ax_comp = nexttile(t_bars, stim);

    % Collect data from LA and Astria
    data_regions = [];
    region_labels = {};

    for br = 1:2
        if isempty(results_all{br})
            continue;
        end

        res = results_all{br};

        % Select appropriate PSTH and latencies
        if stim == 1  % CS
            psth_Hz = res.psth_CS_Hz;
            onset_lat = res.CS_onset_lat;
            offset_lat = res.CS_offset_lat;
        elseif stim == 2  % US
            psth_Hz = res.psth_US_Hz;
            onset_lat = res.US_onset_lat;
            offset_lat = res.US_offset_lat;
        else  % CS+US
            psth_Hz = res.psth_Both_Hz;
            onset_lat = res.Both_onset_lat;
            offset_lat = res.Both_offset_lat;
        end

        % Find responsive neurons (have onset and offset detected)
        responsive_idx = find(~isnan(onset_lat) & ~isnan(offset_lat));

        if ~isempty(responsive_idx)
            % Calculate Delta Peak FR for each responsive neuron
            delta_peak_fr = zeros(length(responsive_idx), 1);
            baseline_idx = 1:(g.pre_time / g.bin_time);

            for n = 1:length(responsive_idx)
                idx_n = responsive_idx(n);

                % Calculate baseline firing rate
                baseline_fr = mean(psth_Hz(idx_n, baseline_idx));

                % Calculate peak firing rate in response window
                onset_bin = round(onset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                offset_bin = round(offset_lat(idx_n) / g.bin_time) + 1 + g.roi(1) - 1;
                peak_fr = max(psth_Hz(idx_n, onset_bin:offset_bin));

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
        errorbar([1 2], means_data, sems_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

        % Statistical comparison - Wilcoxon rank-sum test (unpaired)
        [p_val, ~] = ranksum(LA_data, Astria_data);

        % Add significance marker
        y_max = max(means_data + sems_data);
        curr_ylim = ylim;
        y_range = curr_ylim(2) - curr_ylim(1);
        y_pos = y_max + 0.1 * y_range;

        if p_val < 0.001
            sig_text = '***';
        elseif p_val < 0.01
            sig_text = '**';
        elseif p_val < 0.05
            sig_text = '*';
        else
            sig_text = 'n.s.';
        end

        plot([1 2], [y_pos y_pos], 'k-', 'LineWidth', 1.5);
        text(1.5, y_pos + 0.15 * y_range, sig_text, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontSize', 10);

        hold off;

        % Formatting
        xlim([0.5 2.5]);
        ylim([0 100]);
        xticks([1 2]);
        xticklabels({'LA', 'AStria'});

        % Y-label and Y-tick labels only on first panel
        if stim == 1
            ylabel('ΔFR (Hz)', 'FontSize', 10, 'Interpreter', 'tex');
        else
            set(gca, 'YTickLabel', []);
        end

        title(stim_names{stim}, 'FontSize', 10, 'FontWeight', 'bold');
        set(gca, 'FontSize', 10);
        box off;
    end
end

% Add title to bar graph section
title(t_bars, 'Across region comparison', 'FontSize', 12, 'FontWeight', 'bold');

% Add colorbar - manually create with full control, positioned at bottom heatmap (Astria)
drawnow;  % Ensure all positions are updated
cb_width = 0.008;
cb_left = 0.44;
cb_bottom = 0.27;  % Aligned with Astria heatmap (row 2 in 2-row layout)
cb_height = 0.08;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'right');
%ylabel(cb_ax, 'Z-score', 'FontSize', g.fontSize2);
cb_ax.FontSize = g.fontSize2;

%% Chi-square test for cluster proportions between LA and Astria
fprintf('\n=== Chi-square Test: LA vs Astria ===\n');

% Build contingency table: 2 regions × 3 clusters (CS-sel, US-sel, CS+US)
contingency_table = zeros(2, 3);
for br = 1:2
    if ~isempty(results_all{br})
        for c = 1:3
            contingency_table(br, c) = sum(results_all{br}.Clusters == c);
        end
    end
end

fprintf('Contingency Table (regions × clusters):\n');
fprintf('%12s', 'Region', 'CS-sel', 'US-sel', 'CS+US', 'Total');
fprintf('\n');
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
        all_clusters = [all_clusters; results_all{br}.Clusters];
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

%% Create separate figure for bar charts
fig_bars = figure('Units', 'pixels', 'Position', [100, -100, 1000, 1000], 'Visible', 'on');
t_bars = tiledlayout(fig_bars, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Create nested tiledlayout for all metrics (3 columns now)
t_metrics = t_bars;

metric_names = {'ΔnSpikes', 'ΔFR (Hz)', 'Response length (ms)'};

for c = [1 2 3]  % CS-selective, US-selective, Multisensory
    for metric = 1:3  % ΔnSpikes, ΔPeak FR, Response length
        % Offset by 3 to move bars to rows 2-4
        ax_metric = nexttile(t_metrics, c*3 + metric);

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
                elseif metric == 2  % Peak firing rate change (Hz)
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

            % LA bars at positions 2, 3, 4 (narrower bars) - use cluster color (darker)
            bar_color_LA = cluster_colors(c, :) * 0.7;  % Darker version
            bar([2 3 4], means_LA, 0.6, 'FaceColor', bar_color_LA, 'EdgeColor', 'k', 'LineWidth', 1);
            errorbar([2 3 4], means_LA, sems_LA, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);

            % Astria bars at positions 7, 8, 9 (gap at positions 5, 6) - use cluster color (lighter)
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

        % Statistical comparisons - Wilcoxon signed rank test (paired data)
        % Compare CS vs US, CS vs CS+US, US vs CS+US within each region

        % Set fixed spacing for significance levels based on metric type
        if metric == 1  % ΔnSpikes
            sig_spacing = 4;  % Fixed absolute spacing between significance levels
        elseif metric == 2  % ΔPeak FR
            sig_spacing = 12;  % Fixed absolute spacing (10% of 120)
        else  % metric == 3, Response length
            sig_spacing = 50;  % Fixed absolute spacing (10% of 500)
        end

        % LA comparisons
        if ~isempty(data_LA)
            % CS vs US
            [p_LA_CS_US, ~] = signrank(data_LA(:,1), data_LA(:,2));
            % CS vs CS+US
            [p_LA_CS_Both, ~] = signrank(data_LA(:,1), data_LA(:,3));
            % US vs CS+US
            [p_LA_US_Both, ~] = signrank(data_LA(:,2), data_LA(:,3));

            % Add significance markers for LA
            y_max_LA = max(means_LA + sems_LA);
            curr_ylim = ylim;
            y_range = curr_ylim(2) - curr_ylim(1);
            hold on;

            % Level 1: Short comparisons (CS vs US, US vs CS+US)
            % Position above data with fixed absolute distance
            y_pos_level1 = y_max_LA + 0.08 * y_range;

            % CS vs US
            if p_LA_CS_US < 0.001
                sig_text = '***';
            elseif p_LA_CS_US < 0.01
                sig_text = '**';
            elseif p_LA_CS_US < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([2.1 2.9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                text(2.5, y_pos_level1 + 0.05 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            % US vs CS+US
            if p_LA_US_Both < 0.001
                sig_text = '***';
            elseif p_LA_US_Both < 0.01
                sig_text = '**';
            elseif p_LA_US_Both < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([3.1 3.9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                text(3.5, y_pos_level1 + 0.05 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            % Level 2: Long comparison (CS vs CS+US) - fixed absolute spacing above level 1
            y_pos_level2 = y_pos_level1 + sig_spacing;

            % CS vs CS+US
            if p_LA_CS_Both < 0.001
                sig_text = '***';
            elseif p_LA_CS_Both < 0.01
                sig_text = '**';
            elseif p_LA_CS_Both < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([2 4], [y_pos_level2 y_pos_level2], 'k-', 'LineWidth', 1.5);
                text(3, y_pos_level2 + 0.05 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            hold off;
        end

        % Astria comparisons
        if ~isempty(data_Astria)
            % CS vs US
            [p_Astria_CS_US, ~] = signrank(data_Astria(:,1), data_Astria(:,2));
            % CS vs CS+US
            [p_Astria_CS_Both, ~] = signrank(data_Astria(:,1), data_Astria(:,3));
            % US vs CS+US
            [p_Astria_US_Both, ~] = signrank(data_Astria(:,2), data_Astria(:,3));

            % Add significance markers for Astria
            y_max_Astria = max(means_Astria + sems_Astria);
            curr_ylim = ylim;
            y_range = curr_ylim(2) - curr_ylim(1);
            hold on;

            % Level 1: Short comparisons (CS vs US, US vs CS+US)
            % Position above data with fixed absolute distance
            y_pos_level1 = y_max_Astria + 0.08 * y_range;

            % CS vs US
            if p_Astria_CS_US < 0.001
                sig_text = '***';
            elseif p_Astria_CS_US < 0.01
                sig_text = '**';
            elseif p_Astria_CS_US < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([7.1 7.9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                text(7.5, y_pos_level1 + 0.05 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            % US vs CS+US
            if p_Astria_US_Both < 0.001
                sig_text = '***';
            elseif p_Astria_US_Both < 0.01
                sig_text = '**';
            elseif p_Astria_US_Both < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([8.1 8.9], [y_pos_level1 y_pos_level1], 'k-', 'LineWidth', 1.5);
                text(8.5, y_pos_level1 + 0.05 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            % Level 2: Long comparison (CS vs CS+US) - fixed absolute spacing above level 1
            y_pos_level2 = y_pos_level1 + sig_spacing;

            % CS vs CS+US
            if p_Astria_CS_Both < 0.001
                sig_text = '***';
            elseif p_Astria_CS_Both < 0.01
                sig_text = '**';
            elseif p_Astria_CS_Both < 0.05
                sig_text = '*';
            else
                sig_text = '';
            end
            if ~isempty(sig_text)
                plot([7 9], [y_pos_level2 y_pos_level2], 'k-', 'LineWidth', 1.5);
                text(8, y_pos_level2 + 0.05 * y_range, sig_text, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', g.fontSize2);
            end

            hold off;
        end

        % Formatting
        xlim([0.5 10.5]);
        xticks([2 3 4 7 8 9]);

        % Set fixed y-limits based on metric type
        if metric == 1  % ΔnSpikes
            ylim([0 32]);
            y_ticks = [0 16 32];
        elseif metric == 2  % ΔPeak FR
            ylim([0 100]);
            y_ticks = [0 50 100];
        else  % metric == 3, Response length
            ylim([0 500]);
            y_ticks = [0 250 500];
        end

        % Add title on top row (now row 2, so c == 1)
        if c == 1
            title(metric_names{metric}, 'FontSize', g.fontSize1, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end

        % X-tick labels only on bottom row (now c == 3)
        if c == 3
            xticklabels({'CS', 'US', 'CS+US', 'CS', 'US', 'CS+US'});
        else
            set(gca, 'XTickLabel', []);
        end

        % Set y-ticks
        yticks(y_ticks);

        % Add LA/AStria labels on top row (now row 2, so c == 1)
        if c == 1
            curr_ylim = ylim;
            y_pos = curr_ylim(2) * 0.95;  % Position at 95% of y-max
            text(3, y_pos, 'LA', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            text(8, y_pos, 'AStria', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end

        % Add cluster name label only on first column
        if metric == 1
            ylabel(cluster_names{c}, 'FontSize', g.fontSize2, 'FontWeight', 'bold');
        end

        set(gca, 'FontSize', g.fontSize2);
        box off;
    end
end

%% Add chi-square results to row 1 (3 panels)
% Panel 1: Permutation distribution
ax_perm = nexttile(t_bars, 1);
hold on;
histogram(chi2_perm, 50, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'Normalization', 'probability');
xline(chi2_obs, 'r-', 'LineWidth', 2);
hold off;
xlabel('\chi^2 statistic', 'FontSize', g.fontSize2, 'Interpreter', 'tex');
ylabel('Probability', 'FontSize', g.fontSize2);
title('Permutation Test', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
legend({'Null', 'Observed'}, 'Location', 'northwest', 'FontSize', 8);
set(gca, 'FontSize', g.fontSize2);
box off;

% Panel 2: Contingency table with values
ax_table = nexttile(t_bars, 2);
axis off;

% Display contingency table as text
y_start = 0.7;
y_step = 0.2;

% Header row
text(0.05, y_start, 'Region', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
text(0.35, y_start, 'CS-sel', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
text(0.58, y_start, 'US-sel', 'FontSize', g.fontSize2, 'FontWeight', 'bold');
text(0.80, y_start, 'CS+US', 'FontSize', g.fontSize2, 'FontWeight', 'bold');

% Data rows
y_pos = y_start - y_step;
text(0.05, y_pos, 'LA', 'FontSize', g.fontSize2);
text(0.35, y_pos, sprintf('%d', contingency_table(1, 1)), 'FontSize', g.fontSize2, 'HorizontalAlignment', 'center');
text(0.58, y_pos, sprintf('%d', contingency_table(1, 2)), 'FontSize', g.fontSize2, 'HorizontalAlignment', 'center');
text(0.80, y_pos, sprintf('%d', contingency_table(1, 3)), 'FontSize', g.fontSize2, 'HorizontalAlignment', 'center');

y_pos = y_pos - y_step;
text(0.05, y_pos, 'AStria', 'FontSize', g.fontSize2);
text(0.35, y_pos, sprintf('%d', contingency_table(2, 1)), 'FontSize', g.fontSize2, 'HorizontalAlignment', 'center');
text(0.58, y_pos, sprintf('%d', contingency_table(2, 2)), 'FontSize', g.fontSize2, 'HorizontalAlignment', 'center');
text(0.80, y_pos, sprintf('%d', contingency_table(2, 3)), 'FontSize', g.fontSize2, 'HorizontalAlignment', 'center');

% Add title with chi-square value
title(sprintf('Contingency Table (\\chi^2=%.2f)', chi2_obs), 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Interpreter', 'tex');

% Panel 3: Effect size and statistics
ax_stats = nexttile(t_bars, 3);
axis off;

% Display statistics with better spacing
y_start = 0.7;
y_step = 0.18;
y_pos = y_start;

text(0.05, y_pos, sprintf('\\chi^2: %.2f', chi2_obs), 'FontSize', g.fontSize2, 'Interpreter', 'tex');
y_pos = y_pos - y_step;
text(0.05, y_pos, sprintf('p: %.4f', p_perm), 'FontSize', g.fontSize2);
y_pos = y_pos - y_step;
text(0.05, y_pos, sprintf('V: %.3f', cramers_v), 'FontSize', g.fontSize2);
y_pos = y_pos - y_step;

if p_perm < 0.001
    sig_text = '***';
    text_color = [0.8 0.2 0.2];
elseif p_perm < 0.01
    sig_text = '**';
    text_color = [0.8 0.2 0.2];
elseif p_perm < 0.05
    sig_text = '*';
    text_color = [0.8 0.2 0.2];
else
    sig_text = 'n.s.';
    text_color = [0.4 0.4 0.4];
end
text(0.05, y_pos, sig_text, 'FontSize', 14, 'FontWeight', 'bold', 'Color', text_color);

title('Statistical Summary', 'FontSize', g.fontSize2, 'FontWeight', 'bold');

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
