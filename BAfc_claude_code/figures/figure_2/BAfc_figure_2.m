% BAfc_figure_2 - OPTIMIZED VERSION
% Compare 4 brain regions (LA, BA, Astria, CeA) with heatmaps and cluster firing rates
% Optimizations:
% - Eliminated duplicate sorting code (was repeated at lines 195-243 and 382-441)
% - Consolidated baseline calculations
% - Created helper function for cluster sorting
% - Vectorized contingency table building
% - Reduced redundant calculations in US metrics section
% - Consolidated cluster line plotting

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

ttl = {'triptest_sound_only','triptest_shocks_only'};
hmptitles = {'CS', 'US'};

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
g.excitation_threshold = 1.5;
g.inhibition_fr_drop = 0.50;
g.inhibition_window_short = 0.2;
g.inhibition_window_long = 1.0;
g.use_percentile = true;
g.clim_percentile = 99;
g.onset_threshold = g.excitation_threshold;
g.bin_time = 0.001;
g.min_consec_bins = max(1, round(0.020 / g.bin_time));
g.smoothvalue = 201;
g.plotwin = [2 2];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
g.min_cluster_percent = 2.1;

% Pre-calculate commonly used indices
g.baseline_idx = 1:(g.pre_time / g.bin_time);
g.inhibition_idx_short = g.pre_time/g.bin_time+1:(g.pre_time+g.inhibition_window_short)/g.bin_time;
g.inhibition_idx_long = g.pre_time/g.bin_time+1:(g.pre_time+g.inhibition_window_long)/g.bin_time;
g.plot_idx = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);

brain_regions = {'LA', 'BA', 'Astria', 'CeA'};
cell_type_filter = {'all', 'all', 'all', 'all'};

cluster_colors = [
    0.8 0.2 0.2;    % CS-selective
    0.2 0.4 0.8;    % US-selective
    0.6 0.2 0.6;    % Multisensory
    0.6 0.6 0.6;    % Non-responsive
    0.2 0.6 0.3     % Inhibited
];

cluster_names = {'CS-sel', 'US-sel', 'Multi', 'Non-resp', 'Inhibited'};

%% Calculate PSTHs once
psthZ_full = cell(1, numel(ttl));
psthHz_full = cell(1, numel(ttl));

% Calculate Savitzky-Golay filter delay correction
filter_delay = floor(g.smoothvalue / 2);  % Half the filter width (symmetric kernel)

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
    baseline_mean = mean(psth_spx_og(:, g.baseline_idx), 2);
    baseline_std = std(psth_spx_og(:, g.baseline_idx), 0, 2);
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;
    psth_spx_smooth = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);

    % Correct for filter delay by shifting forward (move data earlier in time)
    psth_spx_corrected = zeros(size(psth_spx_smooth));
    psth_spx_corrected(:, filter_delay+1:end) = psth_spx_smooth(:, 1:end-filter_delay);
    psth_spx_corrected(:, 1:filter_delay) = repmat(psth_spx_smooth(:, 1), 1, filter_delay);
    psthZ_full{hmp} = psth_spx_corrected;
end

%% Process each brain region
results_all = cell(1, 4);

for br = 1:4

    % Get neuron indices
    if strcmp(cell_type_filter{br}, 'PN')
        idx_neurons = strcmp(g.cell_metrics.brainRegion, brain_regions{br}) & strcmp(g.cell_metrics.putativeCellType, 'PN');
    else
        idx_neurons = strcmp(g.cell_metrics.brainRegion, brain_regions{br});
    end

    n_neurons = sum(idx_neurons);

    if n_neurons == 0
        continue;
    end

    % Store indices and animal IDs for later use
    results_all{br}.idx_neurons = find(idx_neurons);
    results_all{br}.animals = g.cell_metrics.animal(idx_neurons);

    % Extract PSTHs
    psth_CS = psthZ_full{1}(idx_neurons, :);
    psth_US = psthZ_full{2}(idx_neurons, :);
    psth_CS_Hz = psthHz_full{1}(idx_neurons, :);
    psth_US_Hz = psthHz_full{2}(idx_neurons, :);

    % Calculate responses
    CS_peak = max(psth_CS(:, g.roi), [], 2);
    US_peak = max(psth_US(:, g.roi), [], 2);
    CS_peak_Hz = max(psth_CS_Hz(:, g.roi), [], 2);
    US_peak_Hz = max(psth_US_Hz(:, g.roi), [], 2);

    CS_baseline_fr = mean(psth_CS_Hz(:, g.baseline_idx), 2);
    US_baseline_fr = mean(psth_US_Hz(:, g.baseline_idx), 2);

    % Rule 1: Short window inhibition
    CS_test_fr_short = mean(psth_CS_Hz(:, g.inhibition_idx_short), 2);
    US_test_fr_short = mean(psth_US_Hz(:, g.inhibition_idx_short), 2);
    CS_fr_drop_short = (CS_baseline_fr - CS_test_fr_short) ./ (CS_baseline_fr + eps);
    US_fr_drop_short = (US_baseline_fr - US_test_fr_short) ./ (US_baseline_fr + eps);

    % Rule 2: Long window inhibition
    CS_test_fr_long = mean(psth_CS_Hz(:, g.inhibition_idx_long), 2);
    US_test_fr_long = mean(psth_US_Hz(:, g.inhibition_idx_long), 2);
    CS_fr_drop_long = (CS_baseline_fr - CS_test_fr_long) ./ (CS_baseline_fr + eps);
    US_fr_drop_long = (US_baseline_fr - US_test_fr_long) ./ (US_baseline_fr + eps);

    % Classify neurons
    CS_excited = CS_peak >= g.excitation_threshold;
    US_excited = US_peak >= g.excitation_threshold;

    % Two-rule inhibition
    CS_inhibited_rule1 = CS_fr_drop_short >= g.inhibition_fr_drop;
    CS_inhibited_rule2 = CS_fr_drop_long >= g.inhibition_fr_drop;
    CS_inhibited = CS_inhibited_rule1 | CS_inhibited_rule2;

    US_inhibited_rule1 = US_fr_drop_short >= g.inhibition_fr_drop;
    US_inhibited_rule2 = US_fr_drop_long >= g.inhibition_fr_drop;
    US_inhibited = US_inhibited_rule1 | US_inhibited_rule2;

    Clusters = zeros(n_neurons, 1);
    Clusters(CS_excited & ~US_excited) = 1;
    Clusters(US_excited & ~CS_excited) = 2;
    Clusters(CS_excited & US_excited) = 3;
    Clusters(~CS_excited & ~US_excited & (CS_inhibited | US_inhibited)) = 5;
    Clusters(~CS_excited & ~US_excited & ~CS_inhibited & ~US_inhibited) = 4;

    % Compute latencies
    CS_onset_lat = nan(n_neurons, 1);
    CS_offset_lat = nan(n_neurons, 1);
    US_onset_lat = nan(n_neurons, 1);
    US_offset_lat = nan(n_neurons, 1);

    for n = 1:n_neurons
        [CS_onset_lat(n), CS_offset_lat(n)] = compute_onset_offset_latency(psth_CS(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [US_onset_lat(n), US_offset_lat(n)] = compute_onset_offset_latency(psth_US(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
    end

    % Store results (will sort after merging)
    results_all{br}.Clusters_original = Clusters;
    results_all{br}.Clusters = Clusters;
    results_all{br}.CS_onset_lat = CS_onset_lat;
    results_all{br}.CS_offset_lat = CS_offset_lat;
    results_all{br}.US_onset_lat = US_onset_lat;
    results_all{br}.US_offset_lat = US_offset_lat;
    results_all{br}.CS_peak = CS_peak;
    results_all{br}.US_peak = US_peak;
    results_all{br}.CS_peak_Hz = CS_peak_Hz;
    results_all{br}.US_peak_Hz = US_peak_Hz;
    results_all{br}.psth_CS = psth_CS;
    results_all{br}.psth_US = psth_US;
    results_all{br}.psth_CS_Hz = psth_CS_Hz;
    results_all{br}.psth_US_Hz = psth_US_Hz;
    results_all{br}.n_neurons = n_neurons;
end

%% Merge small clusters within each brain region
for br = 1:4
    if isempty(results_all{br})
        continue;
    end

    Clusters = results_all{br}.Clusters_original;
    psth_CS = results_all{br}.psth_CS;
    psth_US = results_all{br}.psth_US;
    n_neurons = results_all{br}.n_neurons;

    min_cluster_size = ceil(n_neurons * g.min_cluster_percent / 100);

    % Check cluster sizes
    unique_clusters = unique(Clusters);
    small_clusters_in_region = [];

    for c = unique_clusters'
        n_in_cluster = sum(Clusters == c);
        if n_in_cluster < min_cluster_size
            small_clusters_in_region = [small_clusters_in_region; c];
        end
    end

    % Merge small clusters
    if ~isempty(small_clusters_in_region)
        Clusters_merged = Clusters;

        % Calculate centroids for all clusters
        centroids = zeros(5, size(psth_CS, 2) + size(psth_US, 2));
        for c = 1:5
            if sum(Clusters == c) > 0
                psth_concat = [psth_CS(Clusters == c, :), psth_US(Clusters == c, :)];
                centroids(c, :) = mean(psth_concat, 1);
            end
        end

        for sc = small_clusters_in_region'
            large_clusters = setdiff(unique_clusters, small_clusters_in_region);

            if isempty(large_clusters)
                Clusters_merged(Clusters == sc) = 4;
                continue;
            end

            % Calculate distances to large clusters
            distances = zeros(length(large_clusters), 1);
            for i = 1:length(large_clusters)
                lc = large_clusters(i);
                distances(i) = norm(centroids(sc, :) - centroids(lc, :));
            end

            [~, min_idx] = min(distances);
            closest_cluster = large_clusters(min_idx);

            Clusters_merged(Clusters == sc) = closest_cluster;
        end

        results_all{br}.Clusters = Clusters_merged;
    end
end

%% Build contingency table and sort neurons (CONSOLIDATED)

% Build contingency table using vectorized approach
contingency_table = zeros(4, 5);
region_names_stat = {'LA', 'BA', 'AStria', 'CeA'};

for br = 1:4
    if ~isempty(results_all{br})
        % Vectorized cluster counting
        for c = 1:5
            contingency_table(br, c) = sum(results_all{br}.Clusters == c);
        end

        % Sort neurons using helper function
        leafOrder = sort_neurons_by_cluster(results_all{br}, g);
        results_all{br}.leafOrder = leafOrder;
    end
end

contingency_table_for_stats = contingency_table;

%% Create figure
fig = figure('Position', [100, 100, 1000, 700], 'Units', 'pixels');
t = tiledlayout(fig, 4, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add panel label B
annotation(fig, 'textbox', [0.01 0.95 0.05 0.05], 'String', 'B', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Determine global color limits
all_values = [];
for br = 1:4
    if ~isempty(results_all{br})
        psth_sorted = results_all{br}.psth_CS(results_all{br}.leafOrder, :);
        matrix = psth_sorted(:, g.plot_idx);
        all_values = [all_values; matrix(:)];

        psth_sorted = results_all{br}.psth_US(results_all{br}.leafOrder, :);
        matrix = psth_sorted(:, g.plot_idx);
        all_values = [all_values; matrix(:)];
    end
end

clim_min = prctile(all_values, (100 - g.clim_percentile) / 2);
clim_max = prctile(all_values, 100 - (100 - g.clim_percentile) / 2);

%% Plot each brain region
for br = 1:4
    if isempty(results_all{br})
        continue;
    end

    res = results_all{br};

    % Nested tiledlayout for CS and US heatmaps
    t_heatmaps = tiledlayout(t, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_heatmaps.Layout.Tile = (br-1)*5 + 1;
    t_heatmaps.Layout.TileSpan = [1 2];

    % Pre-calculate cluster boundaries (used by both heatmaps)
    Clusters_sorted = res.Clusters(res.leafOrder);
    n_clu = find(diff(Clusters_sorted) ~= 0);

    % Heatmap CS
    ax1 = nexttile(t_heatmaps, 1);
    psth_sorted = res.psth_CS(res.leafOrder, :);
    matrix = psth_sorted(:, g.plot_idx);
    imagesc(g.timeaxis_hmp, 1:size(matrix, 1), matrix);
    clim([clim_min clim_max]);
    colormap(ax1, g.colors.Heatmap);
    xline(0, '--k', 'LineWidth', g.xlinewidth);

    n_neurons_plot = size(matrix, 1);
    yticks([1, n_neurons_plot]);

    ylabel_names = {'LA', 'BA', 'AStria', 'CeA'};
    ylabel(sprintf('%s neurons', ylabel_names{br}), 'FontSize', g.fontSize2, 'FontWeight', 'bold');

    if br == 1
        title('CS', 'FontSize', g.fontSize1);
    end
    if br == 4
        xlabel('Time (s)', 'FontSize', g.fontSize2);
        xticks([-1 0 1]);
    else
        set(gca, 'XTickLabel', []);
    end

    set(gca, 'FontSize', g.fontSize2);

    % Add cluster lines (consolidated function)
    add_cluster_lines(n_clu);

    % Heatmap US
    ax2 = nexttile(t_heatmaps, 2);
    psth_sorted = res.psth_US(res.leafOrder, :);
    matrix = psth_sorted(:, g.plot_idx);
    imagesc(g.timeaxis_hmp, 1:size(matrix, 1), matrix);
    clim([clim_min clim_max]);
    colormap(ax2, g.colors.Heatmap);
    xline(0, '--k', 'LineWidth', g.xlinewidth);

    n_neurons_plot = size(matrix, 1);
    yticks([1, n_neurons_plot]);
    set(gca, 'YTickLabel', []);

    if br == 1
        title('US', 'FontSize', g.fontSize1);
    end
    if br == 4
        xlabel('Time (s)', 'FontSize', g.fontSize2);
        xticks([-1 0 1]);
    else
        set(gca, 'XTickLabel', []);
    end

    set(gca, 'FontSize', g.fontSize2);

    % Add cluster lines (reuse same boundaries)
    add_cluster_lines(n_clu);

    % Stacked bar chart - column 3
    ax = nexttile(t, (br-1)*5 + 3);

    cluster_order = [1, 2, 3, 4, 5];
    cluster_counts = zeros(5, 1);
    for i = 1:5
        c = cluster_order(i);
        cluster_counts(i) = sum(res.Clusters == c);
    end

    cluster_props = cluster_counts / res.n_neurons * 100;

    b = bar(0.5, flipud(cluster_props)', 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.5);

    for i = 1:5
        c = cluster_order(6-i);
        b(i).CData = cluster_colors(c, :);
    end

    xlim([0 2]);
    ylim([0 100]);
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    box off;

    % Add percentage labels
    cumulative = 0;
    for i = 1:5
        c = cluster_order(6-i);
        prop = cluster_props(c);
        if round(prop) > 0
            text(0.5, cumulative + prop/2, sprintf('%.0f%%', prop), ...
                'HorizontalAlignment', 'center', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Color', 'w');
        end
        cumulative = cumulative + prop;
    end

    % Add colored text labels
    y_positions = [85, 70, 50, 30, 10];
    for i = 1:5
        c = cluster_order(i);
        if sum(res.Clusters == c) > 0
            text(1.3, y_positions(i), cluster_names{c}, ...
                'FontSize', g.fontSize2, 'FontWeight', 'bold', ...
                'Color', cluster_colors(c, :), 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle');
        end
    end

    % Cluster PSTHs - nested tiledlayout in columns 4-5
    t_outer = tiledlayout(t, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_outer.Layout.Tile = (br-1)*5 + 4;
    t_outer.Layout.TileSpan = [1 2];

    % CS lineplots
    t_nested_CS = tiledlayout(t_outer, 5, 1, 'TileSpacing', 'none', 'Padding', 'tight');
    t_nested_CS.Layout.Tile = 1;
    if br == 1
        t_nested_CS.Title.String = 'CS';
        t_nested_CS.Title.FontSize = g.fontSize1;
        t_nested_CS.Title.FontWeight = 'bold';
    end

    % US lineplots
    t_nested_US = tiledlayout(t_outer, 5, 1, 'TileSpacing', 'none', 'Padding', 'tight');
    t_nested_US.Layout.Tile = 2;
    if br == 1
        t_nested_US.Title.String = 'US';
        t_nested_US.Title.FontSize = g.fontSize1;
        t_nested_US.Title.FontWeight = 'bold';
    end

    time_vec = g.timeaxis_hmp;
    if length(time_vec) ~= length(g.plot_idx)
        time_vec = linspace(-g.plotwin(1), g.plotwin(2), length(g.plot_idx));
    end

    % Calculate global y-limits for this region
    y_max_resp = 0;
    y_min_resp = 0;
    y_max_inhib = 0;
    y_min_inhib = 0;

    for c = [1 2 3 4]
        clust_idx = find(res.Clusters == c);
        if ~isempty(clust_idx)
            psth_CS_mean = mean(res.psth_CS_Hz(clust_idx, g.plot_idx), 1);
            psth_US_mean = mean(res.psth_US_Hz(clust_idx, g.plot_idx), 1);
            y_max_resp = max([y_max_resp, max(psth_CS_mean), max(psth_US_mean)]);
            y_min_resp = min([y_min_resp, min(psth_CS_mean), min(psth_US_mean)]);
        end
    end

    clust_idx = find(res.Clusters == 5);
    if ~isempty(clust_idx)
        psth_CS_mean = mean(res.psth_CS_Hz(clust_idx, g.plot_idx), 1);
        psth_US_mean = mean(res.psth_US_Hz(clust_idx, g.plot_idx), 1);
        y_max_inhib = max([y_max_inhib, max(psth_CS_mean), max(psth_US_mean)]);
        y_min_inhib = min([y_min_inhib, min(psth_CS_mean), min(psth_US_mean)]);
    end

    % Plot all 5 clusters
    for c = [1 2 3 4 5]
        clust_idx = find(res.Clusters == c);

        % CS panel
        ax_CS = nexttile(t_nested_CS, c);
        hold on;
        xline(0, '--k', 'LineWidth', g.xlinewidth);
        if ~isempty(clust_idx)
            psth_CS_mean = mean(res.psth_CS_Hz(clust_idx, g.plot_idx), 1);
            plot(time_vec, psth_CS_mean, '-', 'Color', cluster_colors(c, :), 'LineWidth', 2.5);
        end
        hold off;
        xlim([time_vec(1) time_vec(end)]);
        ylim(ternary(c == 5, [y_min_inhib*2.2 y_max_inhib*2.2], [y_min_resp*1.1 y_max_resp*1.1]));

        if br == 4 && c == 5
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            set(gca, 'XColor', 'k');
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
            set(gca, 'XColor', 'none');
        end
        set(gca, 'YTickLabel', [], 'YColor', 'none', 'FontSize', g.fontSize2);

        % US panel
        ax_US = nexttile(t_nested_US, c);
        hold on;
        xline(0, '--k', 'LineWidth', g.xlinewidth);
        if ~isempty(clust_idx)
            psth_US_mean = mean(res.psth_US_Hz(clust_idx, g.plot_idx), 1);
            plot(time_vec, psth_US_mean, '-', 'Color', cluster_colors(c, :), 'LineWidth', 2.5);
        end
        hold off;
        xlim([time_vec(1) time_vec(end)]);
        ylim(ternary(c == 5, [y_min_inhib*2.2 y_max_inhib*2.2], [y_min_resp*1.1 y_max_resp*1.1]));

        if br == 4 && c == 5
            xlabel('Time (s)', 'FontSize', g.fontSize2);
            set(gca, 'XColor', 'k');
            xticks([-1 0 1]);
        else
            set(gca, 'XTickLabel', []);
            set(gca, 'XColor', 'none');
        end
        set(gca, 'YTickLabel', [], 'YColor', 'none', 'FontSize', g.fontSize2);

        % Add scalebars
        if c == 4
            curr_ylim = ylim;
            y_range = curr_ylim(2) - curr_ylim(1);
            scalebar_size = 20;
            x_pos = time_vec(end) - 0.85;  % Moved slightly more left
            y_pos = curr_ylim(2) - scalebar_size - 0.05*y_range;

            hold on;
            plot([x_pos x_pos], [y_pos y_pos+scalebar_size], 'k-', 'LineWidth', 3);
            text(x_pos+0.1, y_pos+scalebar_size/2, sprintf('%dHz', scalebar_size), ...
                'FontSize', g.fontSize2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            hold off;
        end

        if c == 5
            curr_ylim = ylim;
            y_range = curr_ylim(2) - curr_ylim(1);
            scalebar_size = 5;
            x_pos = time_vec(end) - 0.85;  % Moved slightly more left
            y_pos = curr_ylim(2) - scalebar_size - 0.05*y_range;

            hold on;
            plot([x_pos x_pos], [y_pos y_pos+scalebar_size], 'k-', 'LineWidth', 3);
            text(x_pos+0.1, y_pos+scalebar_size/2, sprintf('%dHz', scalebar_size), ...
                'FontSize', g.fontSize2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            hold off;
        end
    end
end

% Add colorbar
drawnow;
cb_width = 0.008;
cb_left = 0.41;
cb_bottom = 0.07;
cb_height = 0.08;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'right');
%ylabel(cb_ax, 'Z-score', 'FontSize', g.fontSize2);
cb_ax.FontSize = g.fontSize2;

%% Statistical Analysis: Chi-square test
[chi2_obs, p_parametric] = calculate_chi_square(contingency_table_for_stats);

% Cramér's V for effect size
n_total = sum(contingency_table_for_stats(:));
df_cramer = min(size(contingency_table_for_stats, 1) - 1, size(contingency_table_for_stats, 2) - 1);
cramers_v = sqrt(chi2_obs / (n_total * df_cramer));

% Permutation test
n_permutations = 10000;
chi2_perm = zeros(n_permutations, 1);

% Pool cluster assignments
all_clusters = [];
region_indices = [];
for br = 1:4
    if ~isempty(results_all{br})
        n_neurons = results_all{br}.n_neurons;
        all_clusters = [all_clusters; results_all{br}.Clusters];
        region_indices = [region_indices; br * ones(n_neurons, 1)];
    end
end

for perm = 1:n_permutations
    shuffled_clusters = all_clusters(randperm(length(all_clusters)));

    perm_table = zeros(4, 5);
    for br = 1:4
        br_mask = region_indices == br;
        br_clusters = shuffled_clusters(br_mask);
        for c = 1:5
            perm_table(br, c) = sum(br_clusters == c);
        end
    end

    chi2_perm(perm) = calculate_chi_square(perm_table);
end

p_perm = sum(chi2_perm >= chi2_obs) / n_permutations;

%% Pairwise region similarity matrix

n_regions = 4;
chi2_distance_matrix = zeros(n_regions, n_regions);
p_value_matrix = ones(n_regions, n_regions);
cramers_v_matrix = zeros(n_regions, n_regions);

pairwise_p_values = [];
pairwise_pairs = [];

for r1 = 1:n_regions
    for r2 = r1+1:n_regions
        pair_table = contingency_table_for_stats([r1, r2], :);

        col_sums = sum(pair_table, 1);
        pair_table = pair_table(:, col_sums > 0);

        if size(pair_table, 2) > 0 && all(sum(pair_table, 2) > 0)
            [chi2_pair, p_pair] = calculate_chi_square(pair_table);
            chi2_distance_matrix(r1, r2) = chi2_pair;
            chi2_distance_matrix(r2, r1) = chi2_pair;
            p_value_matrix(r1, r2) = p_pair;
            p_value_matrix(r2, r1) = p_pair;

            n_pair = sum(pair_table(:));
            df_pair = min(size(pair_table, 1) - 1, size(pair_table, 2) - 1);
            v_pair = sqrt(chi2_pair / (n_pair * df_pair));
            cramers_v_matrix(r1, r2) = v_pair;
            cramers_v_matrix(r2, r1) = v_pair;

            pairwise_p_values = [pairwise_p_values; p_pair];
            pairwise_pairs = [pairwise_pairs; r1, r2];
        end
    end
end

% FDR correction
[~, sort_idx] = sort(pairwise_p_values);
n_tests = length(pairwise_p_values);
p_adjusted = ones(n_tests, 1);

for i = 1:n_tests
    idx = sort_idx(i);
    p_adjusted(idx) = min(pairwise_p_values(idx) * n_tests / i, 1);
end

for i = n_tests-1:-1:1
    idx = sort_idx(i);
    idx_next = sort_idx(i+1);
    p_adjusted(idx) = min(p_adjusted(idx), p_adjusted(idx_next));
end

% Reconstruct matrix
p_adjusted_matrix = ones(n_regions, n_regions);
for i = 1:n_tests
    r1 = pairwise_pairs(i, 1);
    r2 = pairwise_pairs(i, 2);
    p_adjusted_matrix(r1, r2) = p_adjusted(i);
    p_adjusted_matrix(r2, r1) = p_adjusted(i);
end

%% US metrics figure
fig_US = figure('Position', [100, 100, 1000, 300], 'Units', 'pixels');
t_US = tiledlayout(fig_US, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

bar_color = [0.2 0.4 0.7];  % Blue color

% Column 1: Chi-square p-values matrix
ax_pval = nexttile(t_US, 1);

% Add panel label C
annotation(fig_US, 'textbox', [0.01 0.97 0.05 0.05], 'String', 'C', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

color_matrix = ones(4, 4, 3);
for r1 = 1:4
    for r2 = 1:4
        if r1 ~= r2
            p_val = p_value_matrix(r1, r2);
            if p_val < 0.001
                color_matrix(r1, r2, :) = [0.65 0.65 0.65];  % Lighter grey for p<0.001
            elseif p_val < 0.05
                color_matrix(r1, r2, :) = [0.80 0.80 0.80];  % Light grey for 0.001<=p<0.05
            end
        end
    end
end

imagesc(ax_pval, color_matrix);
set(ax_pval, 'XTick', 1:4, 'XTickLabel', region_names_stat);
set(ax_pval, 'YTick', 1:4, 'YTickLabel', region_names_stat);
set(ax_pval, 'FontSize', g.fontSize2);
axis(ax_pval, 'square');
title('Chi-square P-values', 'FontSize', 12, 'FontWeight', 'bold');

% Add p-values and significance
hold(ax_pval, 'on');
for r1 = 1:4
    for r2 = 1:4
        if r1 ~= r2
            p_val = p_value_matrix(r1, r2);

            if p_val < 0.001
                val_str = '<0.001';
                sig_str = '***';
            elseif p_val < 0.01
                val_str = sprintf('%.3f', p_val);
                sig_str = '**';
            elseif p_val < 0.05
                val_str = sprintf('%.3f', p_val);
                sig_str = '*';
            else
                val_str = sprintf('%.3f', p_val);
                sig_str = 'n.s.';
            end

            text(ax_pval, r2, r1, {val_str, sig_str}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Color', 'k');
        else
            text(ax_pval, r2, r1, '-', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
end
hold(ax_pval, 'off');

% Columns 2-3: Nested tiledlayout for the 2 bar graphs
t_US_metrics = tiledlayout(t_US, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
t_US_metrics.Layout.Tile = 2;
t_US_metrics.Layout.TileSpan = [1 2];
title(t_US_metrics, 'Comparison of US-evoked excitatory responses', 'FontSize', 12, 'FontWeight', 'bold');

% Add panel label D
annotation(fig_US, 'textbox', [0.35 0.97 0.05 0.05], 'String', 'D', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

clusters_to_use = [2, 3];
metric_names = {'FR (Hz)', 'Response length (ms)'};

% Pre-calculate all metrics for all regions (use already calculated peak values)
all_metric_data = cell(2, 3);

for metric = 1:2
    for br = 1:3
        if isempty(results_all{br})
            all_metric_data{metric, br} = [];
            continue;
        end

        res = results_all{br};

        % Get US-responsive indices
        clust_idx = [];
        for c = clusters_to_use
            clust_idx = [clust_idx; find(res.Clusters == c)];
        end

        if isempty(clust_idx)
            all_metric_data{metric, br} = [];
            continue;
        end

        if metric == 1  % Peak FR - use already calculated US_peak_Hz
            metric_values = res.US_peak_Hz(clust_idx);
        else  % Response length
            metric_values = zeros(length(clust_idx), 1);
            for n = 1:length(clust_idx)
                idx_n = clust_idx(n);
                onset_lat = res.US_onset_lat(idx_n);
                offset_lat = res.US_offset_lat(idx_n);

                if ~isnan(onset_lat) && ~isnan(offset_lat)
                    metric_values(n) = (offset_lat - onset_lat) * 1000;
                else
                    metric_values(n) = 0;
                end
            end
        end

        all_metric_data{metric, br} = metric_values;
    end
end

% Plot each metric in the nested tiledlayout
for metric = 1:2
    ax = nexttile(t_US_metrics, metric);
    hold on;

    positions = [1 2 3];
    means = zeros(3, 1);
    sems = zeros(3, 1);

    for br = 1:3
        if ~isempty(all_metric_data{metric, br})
            means(br) = mean(all_metric_data{metric, br});
            sems(br) = std(all_metric_data{metric, br}) / sqrt(length(all_metric_data{metric, br}));
        end
    end

    % Plot bars
    for br = 1:3
        if means(br) > 0 || sems(br) > 0
            bar(positions(br), means(br), 0.4, 'FaceColor', bar_color, 'EdgeColor', 'k', 'LineWidth', 1);
        end
    end

    % Plot individual data points as empty grey circles
    for br = 1:3
        if ~isempty(all_metric_data{metric, br})
            n_points = length(all_metric_data{metric, br});
            % Add small jitter to x-position for visibility
            x_jitter = positions(br) + (rand(n_points, 1) - 0.5) * 0.15;
            scatter(x_jitter, all_metric_data{metric, br}, 16, [0.5 0.5 0.5], 'LineWidth', 0.5);
        end
    end

    % Error bars
    for br = 1:3
        if means(br) > 0 || sems(br) > 0
            errorbar(positions(br), means(br), sems(br), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
        end
    end

    % Statistical comparisons
    % Use 95th percentile to avoid compression from extreme outliers
    all_data_values = [];
    for br = 1:3
        if ~isempty(all_metric_data{metric, br})
            all_data_values = [all_data_values; all_metric_data{metric, br}];
        end
    end
    if ~isempty(all_data_values)
        data_95th = prctile(all_data_values, 95);
        y_max = max(max(means + sems), data_95th);
    else
        y_max = max(means + sems);
    end
    y_range = y_max * 1.5;  % Add 50% headroom for significance markers

    comparison_pairs = [1 2; 1 3; 2 3];
    sig_results = [];

    for i = 1:size(comparison_pairs, 1)
        br1 = comparison_pairs(i, 1);
        br2 = comparison_pairs(i, 2);

        if ~isempty(all_metric_data{metric, br1}) && ~isempty(all_metric_data{metric, br2})
            data1 = all_metric_data{metric, br1};
            data2 = all_metric_data{metric, br2};

            if length(data1) > 1 && length(data2) > 1
                [p_val, ~] = ranksum(data1, data2);

                if p_val < 0.05
                    sig_results = [sig_results; br1, br2, p_val];
                end
            end
        end
    end

    % Add significance markers
    if ~isempty(sig_results)
        level_height = y_range * 0.12;
        y_start = y_max + y_range * 0.1;

        for i = 1:size(sig_results, 1)
            br1 = sig_results(i, 1);
            br2 = sig_results(i, 2);
            p_val = sig_results(i, 3);

            y_pos = y_start + (i-1) * level_height;

            plot([positions(br1) positions(br2)], [y_pos y_pos], 'k-', 'LineWidth', 1.5);

            if p_val < 0.001
                sig_text = '***';
            elseif p_val < 0.01
                sig_text = '**';
            else
                sig_text = '*';
            end

            text((positions(br1) + positions(br2)) / 2, y_pos, sig_text, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', g.fontSize2);
        end
    end

    hold off;

    xlim([0.3 3.7]);
    ylim([0 y_range]);
    xticks(positions);
    xticklabels({'LA', 'BA', 'AStria'});
    title(metric_names{metric}, 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Interpreter', 'tex');
    set(gca, 'FontSize', g.fontSize2);
    box off;
end

%% Kruskal-Wallis test for US metrics
kw_results = struct();
metric_names_kw = {'PeakFR', 'ResponseLength'};

for metric = 1:2
    % Prepare data for Kruskal-Wallis test
    group_labels = {};
    all_values = [];

    for br = 1:3
        if ~isempty(all_metric_data{metric, br})
            all_values = [all_values; all_metric_data{metric, br}];
            group_labels = [group_labels; repmat({brain_regions{br}}, length(all_metric_data{metric, br}), 1)];
        end
    end

    if ~isempty(all_values)
        [p_kw, tbl, stats] = kruskalwallis(all_values, group_labels, 'off');

        % Post-hoc Dunn test (only if K-W is significant)
        region_pairs = {[1 2], [1 3], [2 3]};
        pair_names = {'LA vs BA', 'LA vs AStria', 'BA vs AStria'};
        p_values = zeros(3, 1);

        % Only do post-hoc if K-W is significant
        if p_kw >= 0.05
            p_values = NaN(3, 1);
        else
            % Perform Dunn test using multcompare with 'dunn-sidak' method
            comparison_results = multcompare(stats, 'CType', 'dunn-sidak', 'Display', 'off');

            % Extract p-values for our specific comparisons
            % comparison_results format: [group1, group2, lower_ci, estimate, upper_ci, p_value]
            for pair = 1:3
                br1 = region_pairs{pair}(1);
                br2 = region_pairs{pair}(2);

                % Find the row in comparison_results that matches this pair
                match_idx = find((comparison_results(:,1) == br1 & comparison_results(:,2) == br2) | ...
                                 (comparison_results(:,1) == br2 & comparison_results(:,2) == br1));

                if ~isempty(match_idx)
                    p_values(pair) = comparison_results(match_idx, 6);
                else
                    p_values(pair) = NaN;
                end
            end
        end

        % Store results
        kw_results.(metric_names_kw{metric}).p_kw = p_kw;
        kw_results.(metric_names_kw{metric}).p_posthoc = p_values;
        kw_results.(metric_names_kw{metric}).pair_names = pair_names;
    end
end

%% Export statistics
export_figure2_to_excel(results_all, contingency_table_for_stats, chi2_obs, p_perm, cramers_v, ...
    p_value_matrix, cramers_v_matrix, p_adjusted_matrix, kw_results, all_metric_data, ...
    brain_regions, cluster_names, metric_names, g);
exportgraphics(gcf, 'figure_2_stats.png', 'Resolution', 300);
close(gcf)
exportgraphics(gcf, 'figure_2.png', 'Resolution', 300);
close(gcf)

%% Helper functions

function leafOrder = sort_neurons_by_cluster(res, g)
    % Consolidated sorting function to eliminate duplicate code
    Clusters = res.Clusters;
    CS_onset_lat = res.CS_onset_lat;
    CS_offset_lat = res.CS_offset_lat;
    US_onset_lat = res.US_onset_lat;
    US_offset_lat = res.US_offset_lat;
    psth_CS = res.psth_CS;
    psth_US = res.psth_US;

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
        else  % Non-responsive (4) or Inhibited (5)
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
end

function add_cluster_lines(n_clu)
    % Consolidated cluster line plotting
    hold on;
    for i = 1:length(n_clu)
        yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
    end
    hold off;
end

function export_figure2_to_excel(results_all, contingency_table, chi2_obs, p_perm, cramers_v, ...
    p_value_matrix, cramers_v_matrix, p_adjusted_matrix, kw_results, all_metric_data, ...
    brain_regions, cluster_names, metric_names, g)

    output_filename = 'figure_2_data.xlsx';

    if exist(output_filename, 'file')
        delete(output_filename);
    end

    % Extract unique animals
    all_animals = {};
    for br = 1:4
        if ~isempty(results_all{br})
            animals_br = results_all{br}.animals;
            if ~iscell(animals_br)
                animals_br = {animals_br};
            end
            all_animals = [all_animals; animals_br(:)];
        end
    end
    unique_animals = unique(all_animals);

    %% Sheet 1: Summary - Sample sizes
    sheet1 = {};
    sheet1{1,1} = 'FIGURE 2: SAMPLE SIZES';
    sheet1{2,1} = ['Generated: ' datestr(now)];
    sheet1{3,1} = '';
    sheet1{4,1} = 'Number of animals (N)';
    sheet1{4,2} = length(unique_animals);
    sheet1{5,1} = 'Animal IDs';
    sheet1{5,2} = strjoin(unique_animals, ', ');
    sheet1{6,1} = '';

    row_idx = 7;
    sheet1{row_idx,1} = 'Region';
    sheet1{row_idx,2} = 'n neurons';
    row_idx = row_idx + 1;

    total_neurons = 0;
    for br = 1:4
        sheet1{row_idx,1} = brain_regions{br};
        if ~isempty(results_all{br})
            n_neurons = results_all{br}.n_neurons;
            sheet1{row_idx,2} = n_neurons;
            total_neurons = total_neurons + n_neurons;
        else
            sheet1{row_idx,2} = 0;
        end
        row_idx = row_idx + 1;
    end

    sheet1{row_idx,1} = 'Total';
    sheet1{row_idx,2} = total_neurons;

    writecell(sheet1, output_filename, 'Sheet', 'Summary_SampleSizes');

    %% Sheet 2: Panel B - Cluster distributions
    sheet2 = {};
    sheet2{1,1} = 'PANEL B: CLUSTER DISTRIBUTIONS';
    sheet2{2,1} = '';
    sheet2{3,1} = 'Contingency Table';
    sheet2{4,1} = 'Region';
    for c = 1:5
        sheet2{4, c+1} = cluster_names{c};
    end
    sheet2{4,7} = 'Total';

    row_idx = 5;
    for br = 1:4
        sheet2{row_idx,1} = brain_regions{br};
        for c = 1:5
            sheet2{row_idx, c+1} = contingency_table(br, c);
        end
        sheet2{row_idx,7} = sum(contingency_table(br, :));
        row_idx = row_idx + 1;
    end

    sheet2{row_idx,1} = 'Total';
    for c = 1:5
        sheet2{row_idx, c+1} = sum(contingency_table(:, c));
    end
    sheet2{row_idx,7} = sum(contingency_table(:));

    row_idx = row_idx + 2;
    sheet2{row_idx,1} = 'Cluster Proportions by Region (%)';
    row_idx = row_idx + 1;
    sheet2{row_idx,1} = 'Region';
    for c = 1:5
        sheet2{row_idx, c+1} = cluster_names{c};
    end

    row_idx = row_idx + 1;
    for br = 1:4
        sheet2{row_idx,1} = brain_regions{br};
        total_region = sum(contingency_table(br, :));
        if total_region > 0
            for c = 1:5
                pct = 100 * contingency_table(br, c) / total_region;
                sheet2{row_idx, c+1} = pct;
            end
        end
        row_idx = row_idx + 1;
    end

    writecell(sheet2, output_filename, 'Sheet', 'PanelB_Clusters');

    %% Sheet 3: Panel C - Chi-square test
    sheet3 = {};
    sheet3{1,1} = 'PANEL C: CHI-SQUARE TEST';
    sheet3{2,1} = '';
    sheet3{3,1} = 'Overall Chi-square Test';
    sheet3{4,1} = 'Chi-square statistic';
    sheet3{4,2} = chi2_obs;
    sheet3{5,1} = 'Permutation p-value (10,000 permutations)';
    sheet3{5,2} = p_perm;
    sheet3{5,3} = format_significance(p_perm);
    sheet3{6,1} = 'Cramér''s V (effect size)';
    sheet3{6,2} = cramers_v;
    sheet3{7,1} = '';

    sheet3{8,1} = 'Pairwise Comparisons (FDR-corrected)';
    sheet3{9,1} = 'Comparison';
    sheet3{9,2} = 'Raw p-value';
    sheet3{9,3} = 'Adjusted p-value';
    sheet3{9,4} = 'Significance';
    sheet3{9,5} = 'Cramér''s V';

    region_pairs = {[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]};
    pair_names = {'LA vs BA', 'LA vs AStria', 'LA vs CeA', 'BA vs AStria', 'BA vs CeA', 'AStria vs CeA'};

    row_idx = 10;
    for i = 1:length(region_pairs)
        r1 = region_pairs{i}(1);
        r2 = region_pairs{i}(2);
        p_raw = p_value_matrix(r1, r2);
        p_adj = p_adjusted_matrix(r1, r2);
        v_pair = cramers_v_matrix(r1, r2);

        sheet3{row_idx,1} = pair_names{i};
        sheet3{row_idx,2} = p_raw;
        sheet3{row_idx,3} = p_adj;
        sheet3{row_idx,4} = format_significance(p_adj);
        sheet3{row_idx,5} = v_pair;
        row_idx = row_idx + 1;
    end

    writecell(sheet3, output_filename, 'Sheet', 'PanelC_ChiSquare');

    %% Sheet 4: Panel D - US Metrics
    metric_names_full = {'Peak FR (Hz)', 'Response Length (ms)'};
    metric_names_kw = {'PeakFR', 'ResponseLength'};

    sheet4 = {};
    sheet4{1,1} = 'PANEL D: US-EVOKED RESPONSE METRICS';
    sheet4{2,1} = '';

    row_idx = 3;
    for metric = 1:2
        sheet4{row_idx,1} = ['--- ' metric_names_full{metric} ' ---'];
        row_idx = row_idx + 1;

        sheet4{row_idx,1} = 'Sample sizes (US-responsive: US-sel + Multi)';
        row_idx = row_idx + 1;
        sheet4{row_idx,1} = 'Region';
        sheet4{row_idx,2} = 'n';
        row_idx = row_idx + 1;

        for br = 1:3
            sheet4{row_idx,1} = brain_regions{br};
            if ~isempty(all_metric_data{metric, br})
                sheet4{row_idx,2} = length(all_metric_data{metric, br});
            else
                sheet4{row_idx,2} = 0;
            end
            row_idx = row_idx + 1;
        end

        sheet4{row_idx,1} = '';
        row_idx = row_idx + 1;
        sheet4{row_idx,1} = 'Descriptive Statistics';
        row_idx = row_idx + 1;
        sheet4{row_idx,1} = 'Region';
        sheet4{row_idx,2} = 'Mean';
        sheet4{row_idx,3} = 'SEM';
        sheet4{row_idx,4} = 'Median';
        sheet4{row_idx,5} = 'SD';
        row_idx = row_idx + 1;

        for br = 1:3
            if ~isempty(all_metric_data{metric, br})
                data = all_metric_data{metric, br};
                sheet4{row_idx,1} = brain_regions{br};
                sheet4{row_idx,2} = mean(data);
                sheet4{row_idx,3} = std(data)/sqrt(length(data));
                sheet4{row_idx,4} = median(data);
                sheet4{row_idx,5} = std(data);
                row_idx = row_idx + 1;
            end
        end

        if isfield(kw_results, metric_names_kw{metric})
            kw_data = kw_results.(metric_names_kw{metric});

            sheet4{row_idx,1} = '';
            row_idx = row_idx + 1;
            sheet4{row_idx,1} = 'Kruskal-Wallis Test';
            row_idx = row_idx + 1;
            sheet4{row_idx,1} = 'p-value';
            sheet4{row_idx,2} = kw_data.p_kw;
            sheet4{row_idx,3} = format_significance(kw_data.p_kw);
            row_idx = row_idx + 1;

            sheet4{row_idx,1} = '';
            row_idx = row_idx + 1;
            sheet4{row_idx,1} = 'Post-hoc Pairwise Comparisons (Dunn-Šidák)';
            row_idx = row_idx + 1;
            sheet4{row_idx,1} = 'Comparison';
            sheet4{row_idx,2} = 'p-value';
            sheet4{row_idx,3} = 'Significance';
            row_idx = row_idx + 1;

            for pair = 1:3
                if ~isnan(kw_data.p_posthoc(pair))
                    sheet4{row_idx,1} = kw_data.pair_names{pair};
                    sheet4{row_idx,2} = kw_data.p_posthoc(pair);
                    sheet4{row_idx,3} = format_significance(kw_data.p_posthoc(pair));
                    row_idx = row_idx + 1;
                end
            end
        end

        sheet4{row_idx,1} = '';
        row_idx = row_idx + 1;
    end

    writecell(sheet4, output_filename, 'Sheet', 'PanelD_USMetrics');

    %% Sheet 5: Cluster latencies
    sheet5 = {};
    sheet5{1,1} = 'CLUSTER-SPECIFIC LATENCIES';
    sheet5{2,1} = '';

    row_idx = 3;
    for br = 1:4
        if isempty(results_all{br})
            continue;
        end

        sheet5{row_idx,1} = ['--- ' brain_regions{br} ' ---'];
        row_idx = row_idx + 1;

        res = results_all{br};

        for c = 1:5
            clust_idx = res.Clusters == c;
            n_cluster = sum(clust_idx);

            if n_cluster == 0
                continue;
            end

            sheet5{row_idx,1} = [cluster_names{c} ' (n=' num2str(n_cluster) ')'];
            row_idx = row_idx + 1;

            % CS latencies (for CS-sel and Multi)
            if c == 1 || c == 3
                cs_onsets = res.CS_onset_lat(clust_idx);
                cs_onsets = cs_onsets(~isnan(cs_onsets));
                if ~isempty(cs_onsets)
                    sheet5{row_idx,1} = '  CS onset latency (ms)';
                    sheet5{row_idx,2} = mean(cs_onsets)*1000;
                    sheet5{row_idx,3} = std(cs_onsets)*1000;
                    sheet5{row_idx,4} = median(cs_onsets)*1000;
                    sheet5{row_idx,5} = length(cs_onsets);
                    row_idx = row_idx + 1;
                end

                cs_offsets = res.CS_offset_lat(clust_idx);
                cs_offsets = cs_offsets(~isnan(cs_offsets));
                if ~isempty(cs_offsets)
                    sheet5{row_idx,1} = '  CS offset latency (ms)';
                    sheet5{row_idx,2} = mean(cs_offsets)*1000;
                    sheet5{row_idx,3} = std(cs_offsets)*1000;
                    sheet5{row_idx,4} = median(cs_offsets)*1000;
                    sheet5{row_idx,5} = length(cs_offsets);
                    row_idx = row_idx + 1;
                end
            end

            % US latencies (for US-sel and Multi)
            if c == 2 || c == 3
                us_onsets = res.US_onset_lat(clust_idx);
                us_onsets = us_onsets(~isnan(us_onsets));
                if ~isempty(us_onsets)
                    sheet5{row_idx,1} = '  US onset latency (ms)';
                    sheet5{row_idx,2} = mean(us_onsets)*1000;
                    sheet5{row_idx,3} = std(us_onsets)*1000;
                    sheet5{row_idx,4} = median(us_onsets)*1000;
                    sheet5{row_idx,5} = length(us_onsets);
                    row_idx = row_idx + 1;
                end

                us_offsets = res.US_offset_lat(clust_idx);
                us_offsets = us_offsets(~isnan(us_offsets));
                if ~isempty(us_offsets)
                    sheet5{row_idx,1} = '  US offset latency (ms)';
                    sheet5{row_idx,2} = mean(us_offsets)*1000;
                    sheet5{row_idx,3} = std(us_offsets)*1000;
                    sheet5{row_idx,4} = median(us_offsets)*1000;
                    sheet5{row_idx,5} = length(us_offsets);
                    row_idx = row_idx + 1;
                end
            end
        end

        sheet5{row_idx,1} = '';
        row_idx = row_idx + 1;
    end

    writecell(sheet5, output_filename, 'Sheet', 'ClusterLatencies');

    %% Sheet 6: Raw Data - All neurons with US metrics
    sheet6 = {};
    sheet6{1,1} = 'RAW DATA: ALL NEURONS';
    sheet6{2,1} = 'Individual neuron values (all clusters). CS and US peak values in z-score and Hz. Response Length calculated for US-responsive neurons (US-sel + Multi).';
    sheet6{3,1} = '';

    row_idx = 4;

    for br = 1:3
        if isempty(results_all{br})
            continue;
        end

        sheet6{row_idx,1} = ['--- ' brain_regions{br} ' ---'];
        row_idx = row_idx + 1;

        % Header
        sheet6{row_idx,1} = 'Local #';
        sheet6{row_idx,2} = 'Global Index';
        sheet6{row_idx,3} = 'Animal ID';
        sheet6{row_idx,4} = 'Cluster';
        sheet6{row_idx,5} = 'CS Peak (z-score)';
        sheet6{row_idx,6} = 'US Peak (z-score)';
        sheet6{row_idx,7} = 'CS Peak (Hz)';
        sheet6{row_idx,8} = 'US Peak (Hz)';
        sheet6{row_idx,9} = '';  % Empty separator
        sheet6{row_idx,10} = 'Response Length (ms)';
        row_idx = row_idx + 1;

        res = results_all{br};
        n_neurons = res.n_neurons;

        % Calculate Response Length for ALL neurons (not just US-responsive)
        resp_len_all = nan(n_neurons, 1);

        for n = 1:n_neurons
            onset_lat = res.US_onset_lat(n);
            offset_lat = res.US_offset_lat(n);

            if ~isnan(onset_lat) && ~isnan(offset_lat)
                % Response length
                resp_len_all(n) = (offset_lat - onset_lat) * 1000;
            end
        end

        % Loop through all neurons sorted by cluster
        cluster_order = [1, 2, 3, 4, 5];
        local_num = 1;
        for c = cluster_order
            clust_idx = find(res.Clusters == c);
            if isempty(clust_idx)
                continue;
            end

            for i = 1:length(clust_idx)
                local_idx = clust_idx(i);
                global_idx = res.idx_neurons(local_idx);
                animal_id = g.cell_metrics.animal{global_idx};
                cluster_name = cluster_names{c};

                sheet6{row_idx,1} = local_num;
                sheet6{row_idx,2} = global_idx;
                sheet6{row_idx,3} = animal_id;
                sheet6{row_idx,4} = cluster_name;
                sheet6{row_idx,5} = res.CS_peak(local_idx);
                sheet6{row_idx,6} = res.US_peak(local_idx);
                sheet6{row_idx,7} = res.CS_peak_Hz(local_idx);
                sheet6{row_idx,8} = res.US_peak_Hz(local_idx);
                sheet6{row_idx,9} = '';  % Empty separator
                sheet6{row_idx,10} = resp_len_all(local_idx);
                row_idx = row_idx + 1;
                local_num = local_num + 1;
            end
        end

        % Add summary stats for US-responsive neurons only (US-sel + Multi)
        us_resp_idx = res.Clusters == 2 | res.Clusters == 3;
        resp_len_resp = resp_len_all(us_resp_idx);
        resp_len_resp = resp_len_resp(~isnan(resp_len_resp));

        sheet6{row_idx,1} = '';
        row_idx = row_idx + 1;
        sheet6{row_idx,1} = 'Summary (US-responsive only):';
        row_idx = row_idx + 1;
        sheet6{row_idx,1} = 'Mean:';
        if ~isempty(resp_len_resp)
            sheet6{row_idx,10} = mean(resp_len_resp);
        end
        row_idx = row_idx + 1;
        sheet6{row_idx,1} = 'SEM:';
        if ~isempty(resp_len_resp)
            sheet6{row_idx,10} = std(resp_len_resp)/sqrt(length(resp_len_resp));
        end
        row_idx = row_idx + 1;
        sheet6{row_idx,1} = 'Median:';
        if ~isempty(resp_len_resp)
            sheet6{row_idx,10} = median(resp_len_resp);
        end
        row_idx = row_idx + 1;
        sheet6{row_idx,1} = 'SD:';
        if ~isempty(resp_len_resp)
            sheet6{row_idx,10} = std(resp_len_resp);
        end
        row_idx = row_idx + 1;

        row_idx = row_idx + 2;
    end

    writecell(sheet6, output_filename, 'Sheet', 'RawData_AllNeurons');
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

function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

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
            % No offset detected - use end of window
            offset_lat = (length(seg) - 1) * bin_time;
        else
            offset_idx = onset_idx + offset_idx_relative - 1;
            offset_lat = (offset_idx - 1) * bin_time;
        end
    end
end
