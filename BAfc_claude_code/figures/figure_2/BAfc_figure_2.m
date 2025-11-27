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
g.excitation_threshold = 2;
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

fprintf('\nInhibition detection (two-rule system):\n');
fprintf('  Rule 1: %.0f%% FR drop in 0-%.3f sec window\n', g.inhibition_fr_drop*100, g.inhibition_window_short);
fprintf('  Rule 2: %.0f%% FR drop in 0-%.3f sec window\n', g.inhibition_fr_drop*100, g.inhibition_window_long);
fprintf('  (Classified as inhibited if EITHER rule is satisfied)\n');

for br = 1:4
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

    fprintf('  CS inhibited: %d total (Rule1: %d, Rule2: %d, Both: %d)\n', ...
        sum(CS_inhibited), sum(CS_inhibited_rule1), sum(CS_inhibited_rule2), ...
        sum(CS_inhibited_rule1 & CS_inhibited_rule2));
    fprintf('  US inhibited: %d total (Rule1: %d, Rule2: %d, Both: %d)\n', ...
        sum(US_inhibited), sum(US_inhibited_rule1), sum(US_inhibited_rule2), ...
        sum(US_inhibited_rule1 & US_inhibited_rule2));

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
    results_all{br}.psth_CS = psth_CS;
    results_all{br}.psth_US = psth_US;
    results_all{br}.psth_CS_Hz = psth_CS_Hz;
    results_all{br}.psth_US_Hz = psth_US_Hz;
    results_all{br}.n_neurons = n_neurons;
end

%% Merge small clusters within each brain region
fprintf('\n--- Merging Small Clusters Within Each Region ---\n');
fprintf('Minimum cluster size: %.1f%% of neurons in each region\n', g.min_cluster_percent);

for br = 1:4
    if isempty(results_all{br})
        continue;
    end

    Clusters = results_all{br}.Clusters_original;
    psth_CS = results_all{br}.psth_CS;
    psth_US = results_all{br}.psth_US;
    n_neurons = results_all{br}.n_neurons;

    min_cluster_size = ceil(n_neurons * g.min_cluster_percent / 100);
    fprintf('\n%s (n=%d): min cluster size = %d neurons (%.1f%%)\n', ...
        brain_regions{br}, n_neurons, min_cluster_size, g.min_cluster_percent);

    % Check cluster sizes
    unique_clusters = unique(Clusters);
    small_clusters_in_region = [];

    for c = unique_clusters'
        n_in_cluster = sum(Clusters == c);
        fprintf('  %s: %d neurons (%.1f%%)', cluster_names{c}, n_in_cluster, n_in_cluster/n_neurons*100);
        if n_in_cluster < min_cluster_size
            fprintf(' (SMALL - will merge)\n');
            small_clusters_in_region = [small_clusters_in_region; c];
        else
            fprintf('\n');
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
                fprintf('  %s: All clusters small, merging %s into Non-resp\n', ...
                    brain_regions{br}, cluster_names{sc});
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

            fprintf('  %s: Merging %s (n=%d) into %s (distance=%.2f)\n', ...
                brain_regions{br}, cluster_names{sc}, sum(Clusters == sc), ...
                cluster_names{closest_cluster}, distances(min_idx));

            Clusters_merged(Clusters == sc) = closest_cluster;
        end

        results_all{br}.Clusters = Clusters_merged;
    else
        fprintf('  No merging needed\n');
    end
end

%% Build contingency table and sort neurons (CONSOLIDATED)
fprintf('\nBuilding contingency table and sorting neurons...\n');

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

fprintf('\nContingency Table After Merging (regions × clusters):\n');
fprintf('%-10s', 'Region');
for c = 1:5
    fprintf('%12s', cluster_names{c});
end
fprintf('%12s\n', 'Total');

for br = 1:4
    fprintf('%-10s', region_names_stat{br});
    for c = 1:5
        fprintf('%12d', contingency_table(br, c));
    end
    fprintf('%12d\n', sum(contingency_table(br, :)));
end

fprintf('%-10s', 'Total');
for c = 1:5
    fprintf('%12d', sum(contingency_table(:, c)));
end
fprintf('%12d\n', sum(contingency_table(:)));

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
    fprintf('    %s: Unique clusters in heatmap: %s\n', brain_regions{br}, mat2str(unique(Clusters_sorted)));

    % Print all neurons in heatmap order
    fprintf('\n=== %s Heatmap Neurons (in display order) ===\n', brain_regions{br});
    for n = 1:length(res.leafOrder)
        global_idx = res.idx_neurons(res.leafOrder(n));
        animal_id = g.cell_metrics.animal{global_idx};
        cell_id = g.cell_metrics.cellID(global_idx);
        brain_region = g.cell_metrics.brainRegion{global_idx};
        cluster_id = res.Clusters(res.leafOrder(n));
        fprintf('  Heatmap row %3d: Animal=%s, CellID=%d, Region=%s, Cluster=%d (%s)\n', ...
            n, animal_id, cell_id, brain_region, cluster_id, cluster_names{cluster_id});
    end
    fprintf('\n');

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
    fprintf('  %s: y_resp=[%.2f, %.2f], y_inhib=[%.2f, %.2f]\n', ...
        brain_regions{br}, y_min_resp*1.1, y_max_resp*1.1, y_min_inhib*1.1, y_max_inhib*1.1);

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
fprintf('\n=== Chi-square Test ===\n');

[chi2_obs, p_parametric] = calculate_chi_square(contingency_table_for_stats);
fprintf('Observed chi-square statistic: %.4f\n', chi2_obs);
fprintf('Parametric p-value: %.4f\n', p_parametric);

% Cramér's V for effect size
n_total = sum(contingency_table_for_stats(:));
df_cramer = min(size(contingency_table_for_stats, 1) - 1, size(contingency_table_for_stats, 2) - 1);
cramers_v = sqrt(chi2_obs / (n_total * df_cramer));
fprintf('Cramér''s V (effect size): %.4f\n', cramers_v);

% Check expected counts
row_totals = sum(contingency_table_for_stats, 2);
col_totals = sum(contingency_table_for_stats, 1);
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

fprintf('Running %d permutations...\n', n_permutations);
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

    if mod(perm, 1000) == 0
        fprintf('  Completed %d/%d permutations\n', perm, n_permutations);
    end
end

p_perm = sum(chi2_perm >= chi2_obs) / n_permutations;
fprintf('\nPermutation p-value: %.4f\n', p_perm);
fprintf('Significance (p < 0.05): %s\n', ternary(p_perm < 0.05, 'YES', 'NO'));

%% Pairwise region similarity matrix
fprintf('\n=== Pairwise Region Similarity (Chi-square Distance) ===\n');

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

            fprintf('%s vs %s: chi2 = %.4f, p = %.4f, Cramér''s V = %.4f\n', ...
                region_names_stat{r1}, region_names_stat{r2}, chi2_pair, p_pair, v_pair);
        end
    end
end

% FDR correction
fprintf('\n--- FDR Correction (Benjamini-Hochberg) ---\n');
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

fprintf('Pairwise comparisons with FDR-corrected p-values:\n');
for i = 1:n_tests
    r1 = pairwise_pairs(i, 1);
    r2 = pairwise_pairs(i, 2);
    fprintf('%s vs %s: p = %.4f, p_adj = %.4f, %s\n', ...
        region_names_stat{r1}, region_names_stat{r2}, ...
        pairwise_p_values(i), p_adjusted(i), ...
        ternary(p_adjusted(i) < 0.05, 'significant', 'n.s.'));
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
fprintf('\n=== Creating US Metrics Figure ===\n');

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
metric_names = {'\DeltaFR (Hz)', 'Response length (ms)'};
metrics_to_plot = [2, 3];  % Only plot DeltaFR and Response length

% Pre-calculate all metrics for all regions (using already computed latencies)
all_metric_data = cell(3, 3);

for metric = 1:3
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

        % Calculate metric using pre-computed latencies
        metric_values = zeros(length(clust_idx), 1);

        for n = 1:length(clust_idx)
            idx_n = clust_idx(n);
            onset_lat = res.US_onset_lat(idx_n);
            offset_lat = res.US_offset_lat(idx_n);

            if ~isnan(onset_lat) && ~isnan(offset_lat)
                onset_bin = round(onset_lat / g.bin_time) + 1 + g.roi(1) - 1;
                offset_bin = round(offset_lat / g.bin_time) + 1 + g.roi(1) - 1;

                if metric == 1  % Delta nSpikes
                    baseline_fr = mean(res.psth_US_Hz(idx_n, g.baseline_idx));
                    response_fr = mean(res.psth_US_Hz(idx_n, onset_bin:offset_bin));
                    response_duration = (offset_bin - onset_bin + 1) * g.bin_time;
                    metric_values(n) = (response_fr - baseline_fr) * response_duration;
                elseif metric == 2  % Delta Peak FR
                    baseline_fr = mean(res.psth_US_Hz(idx_n, g.baseline_idx));
                    peak_fr = max(res.psth_US_Hz(idx_n, onset_bin:offset_bin));
                    metric_values(n) = peak_fr - baseline_fr;
                else  % Response length
                    metric_values(n) = (offset_lat - onset_lat) * 1000;
                end
            else
                metric_values(n) = 0;
            end
        end

        all_metric_data{metric, br} = metric_values;
    end
end

% Plot each metric in the nested tiledlayout
for plot_idx = 1:length(metrics_to_plot)
    metric = metrics_to_plot(plot_idx);
    ax = nexttile(t_US_metrics, plot_idx);
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
    title(metric_names{plot_idx}, 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Interpreter', 'tex');
    set(gca, 'FontSize', g.fontSize2);
    box off;
end

%% Kruskal-Wallis test for US metrics
fprintf('\n=== Kruskal-Wallis Tests for US Metrics ===\n');

kw_results = struct();
metric_names_kw = {'DeltanSpikes', 'DeltaFR', 'ResponseLength'};

for metric = 1:3
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

        fprintf('\n%s:\n', metric_names_kw{metric});
        fprintf('  Kruskal-Wallis p = %.4f\n', p_kw);

        % Post-hoc pairwise Mann-Whitney U tests (only if K-W is significant for DeltaFR)
        region_pairs = {[1 2], [1 3], [2 3]};
        pair_names = {'LA vs BA', 'LA vs AStria', 'BA vs AStria'};
        p_values = zeros(3, 1);

        % For DeltaFR (metric == 2), only do post-hoc if K-W is significant
        if metric == 2 && p_kw >= 0.05
            fprintf('  K-W non-significant, skipping post-hoc comparisons\n');
            p_values = NaN(3, 1);
        else
            for pair = 1:3
                br1 = region_pairs{pair}(1);
                br2 = region_pairs{pair}(2);

                if ~isempty(all_metric_data{metric, br1}) && ~isempty(all_metric_data{metric, br2})
                    p_values(pair) = ranksum(all_metric_data{metric, br1}, all_metric_data{metric, br2});
                    fprintf('    %s: p = %.4f\n', pair_names{pair}, p_values(pair));
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

%% Combined statistical figures
fprintf('\n=== Creating Combined Statistical Figures ===\n');

fig_stats = figure('Position', [200, 100, 1000, 500]);
t_stats = tiledlayout(fig_stats, 2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add panel label A for supplementary figure
annotation(fig_stats, 'textbox', [0.01 0.97 0.05 0.05], 'String', 'A', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Panel 1: Permutation distribution
ax1 = nexttile(t_stats, 1);
histogram(ax1, chi2_perm, 50, 'Normalization', 'probability', 'FaceColor', [0.7 0.7 0.7]);
hold on;
xline(chi2_obs, 'r-', 'LineWidth', 2);
xlabel('Chi-square statistic', 'FontSize', g.fontSize2);
ylabel('Probability', 'FontSize', g.fontSize2);
title(sprintf('Permutation Test (p=%.4f)', p_perm), 'FontSize', g.fontSize1);
legend({'Permuted', 'Observed'}, 'Location', 'northeast');
hold off;

% Panels 2-3: Contingency table (spanning 2 columns)
ax_table = nexttile(t_stats, 2, [1 2]);
axis off;

% Add panel label B for contingency table
annotation(fig_stats, 'textbox', [0.25 0.97 0.05 0.05], 'String', 'B', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

table_text = {'Region', 'CS-sel', 'US-sel', 'Multi', 'Non-resp', 'Inhib'};
y_pos = 0.85;  % Start lower to avoid overlap with title
y_step = 0.15;

% Define x positions for each column to ensure proper spacing
x_positions = [0.05, 0.22, 0.38, 0.54, 0.70, 0.86];

% Header row
for col = 1:6
    text(ax_table, x_positions(col), y_pos, table_text{col}, ...
        'FontSize', g.fontSize2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Data rows
y_pos = y_pos - y_step;
for br = 1:4
    text(ax_table, x_positions(1), y_pos, region_names_stat{br}, ...
        'FontSize', g.fontSize2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    for c = 1:5
        text(ax_table, x_positions(c+1), y_pos, sprintf('%d', contingency_table(br, c)), ...
            'FontSize', g.fontSize2, 'HorizontalAlignment', 'center');
    end
    y_pos = y_pos - y_step;
end

xlim(ax_table, [0 1.0]);
ylim(ax_table, [0 1.0]);
title(ax_table, sprintf('Contingency Table (\\chi^2=%.2f)', chi2_obs), 'FontSize', g.fontSize1, 'FontWeight', 'bold');

% Panels 4-5: Cramér's V matrix (spanning 2 columns)
ax3 = nexttile(t_stats, 4, [1 2]);

% Add panel label C for Cramer's V
annotation(fig_stats, 'textbox', [0.61 0.97 0.05 0.05], 'String', 'C', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

imagesc(ax3, cramers_v_matrix);
colormap(ax3, parula);
cb3 = colorbar(ax3);
ylabel(cb3, 'Cramér''s V', 'FontSize', g.fontSize2);
set(ax3, 'CLim', [0 0.5]);
set(ax3, 'XTick', 1:4, 'XTickLabel', region_names_stat);
set(ax3, 'YTick', 1:4, 'YTickLabel', region_names_stat);
set(ax3, 'FontSize', g.fontSize2);
axis(ax3, 'square');
title('Cramér''s V (Effect Size)', 'FontSize', g.fontSize1);

hold(ax3, 'on');
for r1 = 1:4
    for r2 = 1:4
        if r1 ~= r2
            val_str = sprintf('%.3f', cramers_v_matrix(r1, r2));
            p_val = p_adjusted_matrix(r1, r2);

            if p_val < 0.001
                sig_str = '***';
            elseif p_val < 0.01
                sig_str = '**';
            elseif p_val < 0.05
                sig_str = '*';
            else
                sig_str = 'n.s.';
            end

            text(ax3, r2, r1, {val_str, sig_str}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Color', 'k');
        else
            text(ax3, r2, r1, '-', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', g.fontSize1, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
end
hold(ax3, 'off');

%% Row 2: Kruskal-Wallis test results
% Add panel label D for Kruskal-Wallis results
annotation(fig_stats, 'textbox', [0.01 0.48 0.05 0.05], 'String', 'D', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Show DeltaFR (metric 2) and Response Length (metric 3)
metrics_to_show = [2, 3];
metric_names_full = {'\DeltanSpikes', '\DeltaFR (Hz)', 'Response length (ms)'};

for idx = 1:length(metrics_to_show)
    metric = metrics_to_show(idx);
    ax_kw = nexttile(t_stats, 5 + idx);

    if isfield(kw_results, metric_names_kw{metric})
        kw_data = kw_results.(metric_names_kw{metric});

        % Create a simple table display
        axis(ax_kw, 'off');

        % Title with overall Kruskal-Wallis p-value
        text(ax_kw, 0.5, 0.95, sprintf('K-W: p=%.4f', kw_data.p_kw), ...
            'FontSize', g.fontSize2, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

        % Post-hoc comparisons (only show if K-W is significant, i.e., not all NaN)
        if ~all(isnan(kw_data.p_posthoc))
            y_pos = 0.75;
            y_step = 0.20;

            for pair = 1:3
                if ~isnan(kw_data.p_posthoc(pair))
                    p_val = kw_data.p_posthoc(pair);

                    if p_val < 0.001
                        sig_str = '***';
                    elseif p_val < 0.01
                        sig_str = '**';
                    elseif p_val < 0.05
                        sig_str = '*';
                    else
                        sig_str = 'ns';
                    end

                    text(ax_kw, 0.5, y_pos, sprintf('%s: %.4f %s', kw_data.pair_names{pair}, p_val, sig_str), ...
                        'FontSize', g.fontSize2 - 1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
                    y_pos = y_pos - y_step;
                end
            end
        end

        xlim(ax_kw, [0 1]);
        ylim(ax_kw, [0 1]);
        title(ax_kw, metric_names_full{metric}, 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Interpreter', 'tex');
    end
end

%% Export statistics
export_figure2_stats(results_all, contingency_table_for_stats, chi2_obs, p_perm, cramers_v, ...
    p_value_matrix, cramers_v_matrix, p_adjusted_matrix, kw_results, all_metric_data, ...
    brain_regions, cluster_names, metric_names, g);

fprintf('\nDone.\n');

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

function export_figure2_stats(results_all, contingency_table, chi2_obs, p_perm, cramers_v, ...
    p_value_matrix, cramers_v_matrix, p_adjusted_matrix, kw_results, all_metric_data, ...
    brain_regions, cluster_names, metric_names, g)

    fid = fopen('figure_2_stats.txt', 'w');

    fprintf(fid, '========================================\n');
    fprintf(fid, 'FIGURE 2 STATISTICS\n');
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, '========================================\n\n');

    %% Sample sizes per region
    fprintf(fid, '### SAMPLE SIZES BY REGION ###\n\n');

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

    fprintf(fid, 'Number of animals (N): %d\n', length(unique_animals));
    fprintf(fid, 'Animal IDs: %s\n\n', strjoin(unique_animals, ', '));

    total_neurons = 0;
    for br = 1:4
        if ~isempty(results_all{br})
            n_neurons = results_all{br}.n_neurons;
            total_neurons = total_neurons + n_neurons;
            fprintf(fid, '%s: n = %d neurons\n', brain_regions{br}, n_neurons);
        else
            fprintf(fid, '%s: n = 0 neurons\n', brain_regions{br});
        end
    end
    fprintf(fid, 'Total neurons: %d\n\n', total_neurons);

    %% Cluster distributions
    fprintf(fid, '### PANEL B: CLUSTER DISTRIBUTIONS BY REGION ###\n\n');

    fprintf(fid, 'Contingency Table (neurons per cluster):\n');
    fprintf(fid, '%-10s', 'Region');
    for c = 1:5
        fprintf(fid, '%-12s', cluster_names{c});
    end
    fprintf(fid, '%-10s\n', 'Total');

    for br = 1:4
        fprintf(fid, '%-10s', brain_regions{br});
        for c = 1:5
            fprintf(fid, '%-12d', contingency_table(br, c));
        end
        fprintf(fid, '%-10d\n', sum(contingency_table(br, :)));
    end
    fprintf(fid, '\n');

    % Proportions
    fprintf(fid, 'Cluster Proportions by Region:\n');
    for br = 1:4
        fprintf(fid, '%s:\n', brain_regions{br});
        total_region = sum(contingency_table(br, :));
        if total_region > 0
            for c = 1:5
                n_cluster = contingency_table(br, c);
                pct = 100 * n_cluster / total_region;
                fprintf(fid, '  %s: %d (%.1f%%)\n', cluster_names{c}, n_cluster, pct);
            end
        end
        fprintf(fid, '\n');
    end

    %% Chi-square test
    fprintf(fid, '### PANEL C (SUPPLEMENTARY): CHI-SQUARE TEST ###\n\n');

    fprintf(fid, 'Overall Chi-square Test:\n');
    fprintf(fid, '  Chi-square statistic: %.4f\n', chi2_obs);
    fprintf(fid, '  Permutation p-value (10,000 permutations): %.4f\n', p_perm);
    fprintf(fid, '  Cramér''s V (effect size): %.4f\n', cramers_v);
    fprintf(fid, '  Significant (p < 0.05): %s\n\n', char(string(p_perm < 0.05)));

    fprintf(fid, 'Pairwise Chi-square Tests (Bonferroni-Holm corrected):\n');
    region_pairs = {[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]};
    pair_names = {'LA vs BA', 'LA vs AStria', 'LA vs CeA', 'BA vs AStria', 'BA vs CeA', 'AStria vs CeA'};

    for i = 1:length(region_pairs)
        r1 = region_pairs{i}(1);
        r2 = region_pairs{i}(2);
        chi2_pair = cramers_v_matrix(r1, r2);
        p_raw = p_value_matrix(r1, r2);
        p_adj = p_adjusted_matrix(r1, r2);
        v_pair = cramers_v_matrix(r1, r2);

        fprintf(fid, '  %s:\n', pair_names{i});
        fprintf(fid, '    Raw p-value: %.4f\n', p_raw);
        fprintf(fid, '    Adjusted p-value (Bonferroni-Holm): %.4f\n', p_adj);
        fprintf(fid, '    Cramér''s V: %.4f\n', v_pair);
        fprintf(fid, '    Significant (adjusted p < 0.05): %s\n', char(string(p_adj < 0.05)));
    end
    fprintf(fid, '\n');

    %% US metrics
    fprintf(fid, '### PANEL D: US-EVOKED RESPONSE METRICS ###\n\n');

    metric_names_full = {'Delta nSpikes', 'Delta Peak FR (Hz)', 'Response Length (ms)'};
    metric_names_kw = {'DeltanSpikes', 'DeltaFR', 'ResponseLength'};

    for metric = 1:3
        fprintf(fid, '--- %s ---\n', metric_names_full{metric});

        fprintf(fid, 'Sample sizes (US-responsive neurons: US-sel + Multi):\n');
        for br = 1:3
            if ~isempty(all_metric_data{metric, br})
                fprintf(fid, '  %s: n = %d\n', brain_regions{br}, length(all_metric_data{metric, br}));
            else
                fprintf(fid, '  %s: n = 0\n', brain_regions{br});
            end
        end

        fprintf(fid, '\nDescriptive Statistics:\n');
        for br = 1:3
            if ~isempty(all_metric_data{metric, br})
                data = all_metric_data{metric, br};
                fprintf(fid, '  %s: %.3f ± %.3f (mean ± SEM), median = %.3f, IQR = [%.3f, %.3f]\n', ...
                    brain_regions{br}, mean(data), std(data)/sqrt(length(data)), ...
                    median(data), prctile(data, 25), prctile(data, 75));
            end
        end

        if isfield(kw_results, metric_names_kw{metric})
            kw_data = kw_results.(metric_names_kw{metric});

            fprintf(fid, '\nStatistical Test: Kruskal-Wallis test\n');
            fprintf(fid, '  p-value: %.4f\n', kw_data.p_kw);
            fprintf(fid, '  Significant (p < 0.05): %s\n', char(string(kw_data.p_kw < 0.05)));

            fprintf(fid, '\nPost-hoc Pairwise Comparisons (Mann-Whitney U test):\n');
            for pair = 1:3
                if ~isnan(kw_data.p_posthoc(pair))
                    fprintf(fid, '  %s: p = %.4f', kw_data.pair_names{pair}, kw_data.p_posthoc(pair));
                    if kw_data.p_posthoc(pair) < 0.001
                        fprintf(fid, ' ***\n');
                    elseif kw_data.p_posthoc(pair) < 0.01
                        fprintf(fid, ' **\n');
                    elseif kw_data.p_posthoc(pair) < 0.05
                        fprintf(fid, ' *\n');
                    else
                        fprintf(fid, ' (n.s.)\n');
                    end
                end
            end
        end
        fprintf(fid, '\n');
    end

    %% Cluster-specific statistics
    fprintf(fid, '### CLUSTER-SPECIFIC STATISTICS ###\n\n');

    for br = 1:4
        if isempty(results_all{br})
            continue;
        end

        fprintf(fid, '--- %s ---\n', brain_regions{br});
        res = results_all{br};

        for c = 1:5
            clust_idx = res.Clusters == c;
            n_cluster = sum(clust_idx);

            if n_cluster == 0
                continue;
            end

            fprintf(fid, '\n%s (n = %d):\n', cluster_names{c}, n_cluster);

            % CS latencies (for CS-sel and Multi)
            if c == 1 || c == 3
                cs_onsets = res.CS_onset_lat(clust_idx);
                cs_onsets = cs_onsets(~isnan(cs_onsets));
                if ~isempty(cs_onsets)
                    fprintf(fid, '  CS onset latency: %.3f ± %.3f ms (mean ± SD), median = %.3f ms, n = %d\n', ...
                        mean(cs_onsets)*1000, std(cs_onsets)*1000, median(cs_onsets)*1000, length(cs_onsets));
                end

                cs_offsets = res.CS_offset_lat(clust_idx);
                cs_offsets = cs_offsets(~isnan(cs_offsets));
                if ~isempty(cs_offsets)
                    fprintf(fid, '  CS offset latency: %.3f ± %.3f ms (mean ± SD), median = %.3f ms, n = %d\n', ...
                        mean(cs_offsets)*1000, std(cs_offsets)*1000, median(cs_offsets)*1000, length(cs_offsets));
                end
            end

            % US latencies (for US-sel and Multi)
            if c == 2 || c == 3
                us_onsets = res.US_onset_lat(clust_idx);
                us_onsets = us_onsets(~isnan(us_onsets));
                if ~isempty(us_onsets)
                    fprintf(fid, '  US onset latency: %.3f ± %.3f ms (mean ± SD), median = %.3f ms, n = %d\n', ...
                        mean(us_onsets)*1000, std(us_onsets)*1000, median(us_onsets)*1000, length(us_onsets));
                end

                us_offsets = res.US_offset_lat(clust_idx);
                us_offsets = us_offsets(~isnan(us_offsets));
                if ~isempty(us_offsets)
                    fprintf(fid, '  US offset latency: %.3f ± %.3f ms (mean ± SD), median = %.3f ms, n = %d\n', ...
                        mean(us_offsets)*1000, std(us_offsets)*1000, median(us_offsets)*1000, length(us_offsets));
                end
            end
        end
        fprintf(fid, '\n');
    end

    fprintf(fid, '========================================\n');
    fprintf(fid, 'END OF FIGURE 2 STATISTICS\n');
    fprintf(fid, '========================================\n');

    fclose(fid);
    fprintf('Statistics exported to: figure_2_stats.txt\n');
end
