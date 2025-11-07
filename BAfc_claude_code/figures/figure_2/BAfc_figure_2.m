% BAfc_figure_2
% Compare 4 brain regions (LA, BA, Astria, CeA) with heatmaps and cluster firing rates
% Left: 4 heatmaps stacked vertically
% Right: 4 rows of firing rate plots for each cluster

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
g.inhibition_fr_drop = 0.50;  % Minimum FR drop threshold
g.inhibition_window_short = 0.2;  % Short time window (sec) for Rule 1
g.inhibition_window_long = 1.0;   % Long time window (sec) for Rule 2 (entire test period)
g.use_percentile = true;
g.clim_percentile = 99;
g.onset_threshold = g.excitation_threshold;
g.bin_time = 0.001;
g.min_consec_bins = max(1, round(0.020 / g.bin_time));
g.smoothvalue = 201;
g.plotwin = [2 2];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
g.min_cluster_percent = 2;  % Minimum cluster size as percentage of neurons in the region (for merging)

% Brain regions and cell types to analyze
brain_regions = {'LA', 'BA', 'Astria', 'CeA'};
cell_type_filter = {'all', 'all', 'all', 'all'};  % All neurons

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
results_all = cell(1, 4);
merge_clusters = false;  % Will be set to true if merging is needed

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

    % Extract PSTHs
    psth_CS = psthZ_full{1}(idx_neurons, :);
    psth_US = psthZ_full{2}(idx_neurons, :);
    psth_CS_Hz = psthHz_full{1}(idx_neurons, :);
    psth_US_Hz = psthHz_full{2}(idx_neurons, :);

    % Calculate responses
    CS_peak = max(psth_CS(:, g.roi), [], 2);
    US_peak = max(psth_US(:, g.roi), [], 2);

    baseline_idx = 1:(g.pre_time / g.bin_time);
    CS_baseline_fr = mean(psth_CS_Hz(:, baseline_idx), 2);
    US_baseline_fr = mean(psth_US_Hz(:, baseline_idx), 2);

    % Rule 1: FR drop in short window (0 to g.inhibition_window_short sec)
    inhibition_idx_short = g.pre_time/g.bin_time+1:(g.pre_time+g.inhibition_window_short)/g.bin_time;
    CS_test_fr_short = mean(psth_CS_Hz(:, inhibition_idx_short), 2);
    US_test_fr_short = mean(psth_US_Hz(:, inhibition_idx_short), 2);
    CS_fr_drop_short = (CS_baseline_fr - CS_test_fr_short) ./ (CS_baseline_fr + eps);
    US_fr_drop_short = (US_baseline_fr - US_test_fr_short) ./ (US_baseline_fr + eps);

    % Rule 2: FR drop in long window (0 to g.inhibition_window_long sec)
    inhibition_idx_long = g.pre_time/g.bin_time+1:(g.pre_time+g.inhibition_window_long)/g.bin_time;
    CS_test_fr_long = mean(psth_CS_Hz(:, inhibition_idx_long), 2);
    US_test_fr_long = mean(psth_US_Hz(:, inhibition_idx_long), 2);
    CS_fr_drop_long = (CS_baseline_fr - CS_test_fr_long) ./ (CS_baseline_fr + eps);
    US_fr_drop_long = (US_baseline_fr - US_test_fr_long) ./ (US_baseline_fr + eps);

    % Classify neurons
    CS_excited = CS_peak >= g.excitation_threshold;
    US_excited = US_peak >= g.excitation_threshold;

    % Two-rule inhibition detection (either rule satisfied = inhibited)
    CS_inhibited_rule1 = CS_fr_drop_short >= g.inhibition_fr_drop;
    CS_inhibited_rule2 = CS_fr_drop_long >= g.inhibition_fr_drop;
    CS_inhibited = CS_inhibited_rule1 | CS_inhibited_rule2;

    US_inhibited_rule1 = US_fr_drop_short >= g.inhibition_fr_drop;
    US_inhibited_rule2 = US_fr_drop_long >= g.inhibition_fr_drop;
    US_inhibited = US_inhibited_rule1 | US_inhibited_rule2;

    % Report inhibition detection breakdown
    fprintf('  CS inhibited: %d total (Rule1: %d, Rule2: %d, Both: %d)\n', ...
        sum(CS_inhibited), sum(CS_inhibited_rule1), sum(CS_inhibited_rule2), ...
        sum(CS_inhibited_rule1 & CS_inhibited_rule2));
    fprintf('  US inhibited: %d total (Rule1: %d, Rule2: %d, Both: %d)\n', ...
        sum(US_inhibited), sum(US_inhibited_rule1), sum(US_inhibited_rule2), ...
        sum(US_inhibited_rule1 & US_inhibited_rule2));

    Clusters = zeros(n_neurons, 1);
    Clusters(CS_excited & ~US_excited) = 1;  % CS-selective
    Clusters(US_excited & ~CS_excited) = 2;  % US-selective
    Clusters(CS_excited & US_excited) = 3;   % Multisensory
    Clusters(~CS_excited & ~US_excited & (CS_inhibited | US_inhibited)) = 5;  % Inhibited
    Clusters(~CS_excited & ~US_excited & ~CS_inhibited & ~US_inhibited) = 4;  % Non-responsive

    % Store original clusters
    results_all{br}.Clusters_original = Clusters;

    % Compute latencies
    CS_onset_lat = nan(n_neurons, 1);
    CS_offset_lat = nan(n_neurons, 1);
    US_onset_lat = nan(n_neurons, 1);
    US_offset_lat = nan(n_neurons, 1);

    for n = 1:n_neurons
        [CS_onset_lat(n), CS_offset_lat(n)] = compute_onset_offset_latency(psth_CS(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [US_onset_lat(n), US_offset_lat(n)] = compute_onset_offset_latency(psth_US(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
    end

    % Sort neurons
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

    % Store results (leafOrder will be recalculated after potential merging)
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

    % Calculate minimum cluster size for this region
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

    % Merge small clusters to nearest cluster based on centroid distance
    if ~isempty(small_clusters_in_region)
        Clusters_merged = Clusters;

        % Calculate centroids for all clusters (using CS and US concatenated)
        centroids = zeros(5, size(psth_CS, 2) + size(psth_US, 2));
        for c = 1:5
            if sum(Clusters == c) > 0
                psth_concat = [psth_CS(Clusters == c, :), psth_US(Clusters == c, :)];
                centroids(c, :) = mean(psth_concat, 1);
            end
        end

        for sc = small_clusters_in_region'
            % Find non-small clusters (potential merge targets)
            large_clusters = setdiff(unique_clusters, small_clusters_in_region);

            if isempty(large_clusters)
                % All clusters are small - merge into cluster 4 (Non-responsive)
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

            % Find closest cluster
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
        results_all{br}.Clusters = Clusters;
    end
end

%% Build contingency table (after merging)
% Build contingency table: regions × clusters
% Rows: LA, BA, Astria, CeA
% Columns: CS-sel, US-sel, Multi, Non-resp, Inhibited
contingency_table = zeros(4, 5);
region_names_stat = {'LA', 'BA', 'AStria', 'CeA'};

for br = 1:4
    if ~isempty(results_all{br})
        for c = 1:5
            contingency_table(br, c) = sum(results_all{br}.Clusters == c);
        end
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

% Use this table for statistical tests
contingency_table_for_stats = contingency_table;
merge_clusters = false;  % No longer need separate merge flag

%% Recalculate leafOrder based on (potentially merged) clusters
fprintf('\nRecalculating neuron sorting order...\n');
for br = 1:4
    if ~isempty(results_all{br})
        Clusters = results_all{br}.Clusters;
        CS_onset_lat = results_all{br}.CS_onset_lat;
        CS_offset_lat = results_all{br}.CS_offset_lat;
        US_onset_lat = results_all{br}.US_onset_lat;
        US_offset_lat = results_all{br}.US_offset_lat;
        psth_CS = results_all{br}.psth_CS;
        psth_US = results_all{br}.psth_US;

        % Sort neurons
        leafOrder = [];
        cluster_order = [1, 2, 3, 4, 5];  % All 5 clusters

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
            elseif c == 4  % Non-responsive (potentially includes merged inhibited)
                mean_zscore = mean([mean(psth_CS(clust_idx, g.roi), 2), mean(psth_US(clust_idx, g.roi), 2)], 2);
                [~, sort_idx] = sort(mean_zscore, 'descend');
                leafOrder = [leafOrder; clust_idx(sort_idx)];
                continue;
            else  % Inhibited (cluster 5, only if not merged)
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

        results_all{br}.leafOrder = leafOrder;
    end
end

%% Create figure
fig = figure('Position', [100, 100, 1000, 700], 'Units', 'pixels');
t = tiledlayout(fig, 4, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

% Determine global color limits
all_values = [];
for br = 1:4
    if ~isempty(results_all{br})
        psth_sorted = results_all{br}.psth_CS(results_all{br}.leafOrder, :);
        matrix = psth_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        all_values = [all_values; matrix(:)];

        psth_sorted = results_all{br}.psth_US(results_all{br}.leafOrder, :);
        matrix = psth_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
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

    % Nested tiledlayout for CS and US heatmaps in columns 1-2 (left, spanning 2 columns)
    t_heatmaps = tiledlayout(t, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_heatmaps.Layout.Tile = (br-1)*5 + 1;
    t_heatmaps.Layout.TileSpan = [1 2];

    % Heatmap CS (left)
    ax1 = nexttile(t_heatmaps, 1);
    psth_sorted = res.psth_CS(res.leafOrder, :);
    matrix = psth_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    imagesc(g.timeaxis_hmp, 1:size(matrix, 1), matrix);
    clim([clim_min clim_max]);
    colormap(ax1, g.colors.Heatmap);
    xline(0, '--k', 'LineWidth', g.xlinewidth);

    % Set yticks to first and last
    n_neurons_plot = size(matrix, 1);
    yticks([1, n_neurons_plot]);

    % Add brain region as ylabel
    ylabel_names = {'LA', 'BA', 'AStria', 'CeA'};
    ylabel(ylabel_names{br}, 'FontSize', g.fontSize2, 'FontWeight', 'bold');

    if br == 1
        title('CS', 'FontSize', g.fontSize1);
    end
    if br == 4
        xlabel('Time (s)', 'FontSize', g.fontSize2);
        xticks([-1 0 1]);
    else
        set(gca, 'XTickLabel', []);
    end

    % Set axis font size
    set(gca, 'FontSize', g.fontSize2);

    % Add cluster lines (black)
    hold on;
    Clusters_sorted = res.Clusters(res.leafOrder);
    n_clu = find(diff(Clusters_sorted) ~= 0);
    fprintf('    %s: Unique clusters in heatmap: %s\n', brain_regions{br}, mat2str(unique(Clusters_sorted)));
    for i = 1:length(n_clu)
        yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
    end
    hold off;

    % Heatmap US (right)
    ax2 = nexttile(t_heatmaps, 2);
    psth_sorted = res.psth_US(res.leafOrder, :);
    matrix = psth_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    imagesc(g.timeaxis_hmp, 1:size(matrix, 1), matrix);
    clim([clim_min clim_max]);
    colormap(ax2, g.colors.Heatmap);
    xline(0, '--k', 'LineWidth', g.xlinewidth);

    % Set yticks to first and last (same as CS)
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

    % Set axis font size
    set(gca, 'FontSize', g.fontSize2);

    % Add cluster lines (black)
    hold on;
    for i = 1:length(n_clu)
        yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
    end
    hold off;

    % Stacked bar chart - column 3 (middle, half width)
    ax = nexttile(t, (br-1)*5 + 3);

    % Determine which clusters exist in this brain region
    unique_clusters = unique(res.Clusters);
    cluster_order = [1, 2, 3, 4, 5];  % All 5 clusters
    n_clusters_display = 5;

    % Count neurons in each cluster (in order they appear in heatmap)
    cluster_counts = zeros(n_clusters_display, 1);
    for i = 1:n_clusters_display
        c = cluster_order(i);
        cluster_counts(i) = sum(res.Clusters == c);
    end

    % Calculate proportions
    cluster_props = cluster_counts / res.n_neurons * 100;

    % Create vertical stacked bar (following Gergo style)
    % Reverse order: bar stacks bottom-up, we want last cluster at bottom, cluster 1 at top
    b = bar(0.5, flipud(cluster_props)', 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.5);

    % Color each segment (reversed order to match)
    for i = 1:n_clusters_display
        c = cluster_order(n_clusters_display+1-i);
        b(i).CData = cluster_colors(c, :);
    end

    % Format axis
    xlim([0 2]);
    ylim([0 100]);

    % Hide all axes
    set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
    box off;

    % Add percentage labels on each segment
    cumulative = 0;
    for i = 1:n_clusters_display
        c = cluster_order(n_clusters_display+1-i);  % Reversed cluster index
        prop = cluster_props(c);
        if round(prop) > 0
            text(0.5, cumulative + prop/2, sprintf('%.0f%%', prop), ...
                'HorizontalAlignment', 'center', 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Color', 'w');
        end
        cumulative = cumulative + prop;
    end

    % Add colored text labels as legend (all rows, matching lineplot positions)
    % Position text labels to the right of bar
    y_positions = [85, 70, 50, 30, 10];  % 5 clusters

    for i = 1:n_clusters_display
        c = cluster_order(i);
        % Only show label if cluster exists in this brain region
        if sum(res.Clusters == c) > 0
            text(1.3, y_positions(i), cluster_names{c}, ...
                'FontSize', g.fontSize2, 'FontWeight', 'bold', ...
                'Color', cluster_colors(c, :), 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle');
        end
    end


    % Cluster PSTHs - nested tiledlayout in columns 4-5 (right, spanning 2 columns)
    % Use 1x2 grid to create CS and US titles
    t_outer = tiledlayout(t, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_outer.Layout.Tile = (br-1)*5 + 4;
    t_outer.Layout.TileSpan = [1 2];

    % CS lineplots (left half)
    t_nested_CS = tiledlayout(t_outer, 5, 1, 'TileSpacing', 'none', 'Padding', 'tight');
    t_nested_CS.Layout.Tile = 1;
    if br == 1
        t_nested_CS.Title.String = 'CS';
        t_nested_CS.Title.FontSize = g.fontSize1;
        t_nested_CS.Title.FontWeight = 'bold';
    end

    % US lineplots (right half)
    t_nested_US = tiledlayout(t_outer, 5, 1, 'TileSpacing', 'none', 'Padding', 'tight');
    t_nested_US.Layout.Tile = 2;
    if br == 1
        t_nested_US.Title.String = 'US';
        t_nested_US.Title.FontSize = g.fontSize1;
        t_nested_US.Title.FontWeight = 'bold';
    end

    plot_idx = round((g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    time_vec = g.timeaxis_hmp;

    % Ensure consistent lengths
    if length(time_vec) ~= length(plot_idx)
        time_vec = linspace(-g.plotwin(1), g.plotwin(2), length(plot_idx));
    end

    % First pass: calculate global y-limits for this brain region (using firing rate Hz)
    % Separate limits for responsive clusters (1-4) and inhibited (5, if not merged)
    y_max_resp = 0;
    y_min_resp = 0;
    y_max_inhib = 0;
    y_min_inhib = 0;

    for c = [1 2 3 4]  % Responsive clusters
        clust_idx = find(res.Clusters == c);
        if ~isempty(clust_idx)
            psth_CS_mean = mean(res.psth_CS_Hz(clust_idx, plot_idx), 1);
            psth_US_mean = mean(res.psth_US_Hz(clust_idx, plot_idx), 1);
            y_max_resp = max([y_max_resp, max(psth_CS_mean), max(psth_US_mean)]);
            y_min_resp = min([y_min_resp, min(psth_CS_mean), min(psth_US_mean)]);
        end
    end

    % Inhibited cluster
    clust_idx = find(res.Clusters == 5);
    if ~isempty(clust_idx)
        psth_CS_mean = mean(res.psth_CS_Hz(clust_idx, plot_idx), 1);
        psth_US_mean = mean(res.psth_US_Hz(clust_idx, plot_idx), 1);
        y_max_inhib = max([y_max_inhib, max(psth_CS_mean), max(psth_US_mean)]);
        y_min_inhib = min([y_min_inhib, min(psth_CS_mean), min(psth_US_mean)]);
    end
    fprintf('  %s: y_resp=[%.2f, %.2f], y_inhib=[%.2f, %.2f]\n', ...
        brain_regions{br}, y_min_resp*1.1, y_max_resp*1.1, y_min_inhib*1.1, y_max_inhib*1.1);

    % Plot all 5 clusters
    clusters_to_plot = [1 2 3 4 5];

    for c = clusters_to_plot
        clust_idx = find(res.Clusters == c);

        % CS panel (left)
        ax_CS = nexttile(t_nested_CS, c);

        % Always draw the xline and set up axis
        hold on;
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        if ~isempty(clust_idx)
            psth_CS_mean = mean(res.psth_CS_Hz(clust_idx, plot_idx), 1);
            plot(time_vec, psth_CS_mean, '-', 'Color', cluster_colors(c, :), 'LineWidth', 2.5);
        end

        hold off;
        xlim([time_vec(1) time_vec(end)]);
        % Use different y-limits for inhibited cluster (doubled)
        if c == 5
            ylim([y_min_inhib*2.2 y_max_inhib*2.2]);
        else
            ylim([y_min_resp*1.1 y_max_resp*1.1]);
        end

        % Only show x-axis on last cluster of last row
        if br == 4 && c == 5
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

        % US panel (right)
        ax_US = nexttile(t_nested_US, c);

        % Always draw the xline and set up axis
        hold on;
        xline(0, '--k', 'LineWidth', g.xlinewidth);

        if ~isempty(clust_idx)
            psth_US_mean = mean(res.psth_US_Hz(clust_idx, plot_idx), 1);
            plot(time_vec, psth_US_mean, '-', 'Color', cluster_colors(c, :), 'LineWidth', 2.5);
        end

        hold off;
        xlim([time_vec(1) time_vec(end)]);
        % Use different y-limits for inhibited cluster (doubled)
        if c == 5
            ylim([y_min_inhib*2.2 y_max_inhib*2.2]);
        else
            ylim([y_min_resp*1.1 y_max_resp*1.1]);
        end

        % Only show x-axis on last cluster of last row
        if br == 4 && c == 5
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

        % Add scalebar on clusters
        if c == 4  % Add scalebar to cluster 4 (non-responsive) for responsive scale
            curr_ylim = ylim;
            y_range = curr_ylim(2) - curr_ylim(1);
            scalebar_size = 20;
            x_pos = time_vec(end) - 0.3;
            y_pos = curr_ylim(2) - scalebar_size - 0.05*y_range;  % Position near top

            hold on;
            plot([x_pos x_pos], [y_pos y_pos+scalebar_size], 'k-', 'LineWidth', 3);
            text(x_pos+0.2, y_pos+scalebar_size/2, sprintf('%d Hz', scalebar_size), ...
                'FontSize', g.fontSize2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            hold off;
        end

        if c == 5  % Add separate scalebar to cluster 5 (inhibited)
            curr_ylim = ylim;
            y_range = curr_ylim(2) - curr_ylim(1);
            scalebar_size = 5;  % Smaller scalebar for inhibited neurons
            x_pos = time_vec(end) - 0.3;
            y_pos = curr_ylim(2) - scalebar_size - 0.05*y_range;  % Position near top

            hold on;
            plot([x_pos x_pos], [y_pos y_pos+scalebar_size], 'k-', 'LineWidth', 3);
            text(x_pos+0.2, y_pos+scalebar_size/2, sprintf('%d Hz', scalebar_size), ...
                'FontSize', g.fontSize2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            hold off;
        end
    end
end

% Add colorbar - manually create with full control
drawnow;  % Ensure all positions are updated
cb_width = 0.010;
cb_left = 0.033;
cb_bottom = 0.01;
cb_height = 0.10;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'left');
ylabel(cb_ax, 'Z-score', 'FontSize', g.fontSize2);
cb_ax.FontSize = g.fontSize2;

%% Statistical Analysis: Chi-square test
fprintf('\n=== Chi-square Test ===\n');

% Calculate observed chi-square statistic
[chi2_obs, p_parametric] = calculate_chi_square(contingency_table_for_stats);
fprintf('Observed chi-square statistic: %.4f\n', chi2_obs);
fprintf('Parametric p-value: %.4f\n', p_parametric);

% Cramér's V for effect size (global)
n_total = sum(contingency_table_for_stats(:));
df_cramer = min(size(contingency_table_for_stats, 1) - 1, size(contingency_table_for_stats, 2) - 1);
cramers_v = sqrt(chi2_obs / (n_total * df_cramer));
fprintf('Cramér''s V (effect size): %.4f\n', cramers_v);

% Check for small expected counts
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

% Pool all cluster assignments
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
    % Shuffle cluster assignments
    shuffled_clusters = all_clusters(randperm(length(all_clusters)));

    % Build permuted contingency table
    perm_table = zeros(4, 5);
    for br = 1:4
        br_mask = region_indices == br;
        br_clusters = shuffled_clusters(br_mask);
        for c = 1:5
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
fprintf('Significance (p < 0.05): %s\n', ternary(p_perm < 0.05, 'YES', 'NO'));


%% Pairwise region similarity matrix (chi-square distance)
fprintf('\n=== Pairwise Region Similarity (Chi-square Distance) ===\n');

% Calculate pairwise chi-square distances
n_regions = 4;
chi2_distance_matrix = zeros(n_regions, n_regions);
p_value_matrix = ones(n_regions, n_regions);
cramers_v_matrix = zeros(n_regions, n_regions);

% Collect p-values for FDR correction
pairwise_p_values = [];
pairwise_pairs = [];

for r1 = 1:n_regions
    for r2 = r1+1:n_regions
        % Extract 2-region contingency table (from merged table if applicable)
        pair_table = contingency_table_for_stats([r1, r2], :);

        % Remove empty columns (clusters with 0 neurons in both regions)
        col_sums = sum(pair_table, 1);
        pair_table = pair_table(:, col_sums > 0);

        if size(pair_table, 2) > 0 && all(sum(pair_table, 2) > 0)
            % Calculate chi-square distance
            [chi2_pair, p_pair] = calculate_chi_square(pair_table);
            chi2_distance_matrix(r1, r2) = chi2_pair;
            chi2_distance_matrix(r2, r1) = chi2_pair;
            p_value_matrix(r1, r2) = p_pair;
            p_value_matrix(r2, r1) = p_pair;

            % Calculate Cramér's V for pairwise effect size
            n_pair = sum(pair_table(:));
            df_pair = min(size(pair_table, 1) - 1, size(pair_table, 2) - 1);
            v_pair = sqrt(chi2_pair / (n_pair * df_pair));
            cramers_v_matrix(r1, r2) = v_pair;
            cramers_v_matrix(r2, r1) = v_pair;

            % Store for FDR correction
            pairwise_p_values = [pairwise_p_values; p_pair];
            pairwise_pairs = [pairwise_pairs; r1, r2];

            fprintf('%s vs %s: chi2 = %.4f, p = %.4f, Cramér''s V = %.4f\n', ...
                region_names_stat{r1}, region_names_stat{r2}, chi2_pair, p_pair, v_pair);
        end
    end
end

% Benjamini-Hochberg FDR correction for pairwise tests
fprintf('\n--- FDR Correction (Benjamini-Hochberg) ---\n');
[~, sort_idx] = sort(pairwise_p_values);
n_tests = length(pairwise_p_values);
p_adjusted = ones(n_tests, 1);

for i = 1:n_tests
    idx = sort_idx(i);
    % BH critical value: p(i) <= (i/m) * alpha
    p_adjusted(idx) = min(pairwise_p_values(idx) * n_tests / i, 1);
end

% Ensure monotonicity (p_adj should be non-decreasing)
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

% Reconstruct p_adjusted_matrix from pairwise results (needed for row 5)
p_adjusted_matrix = ones(n_regions, n_regions);
for i = 1:n_tests
    r1 = pairwise_pairs(i, 1);
    r2 = pairwise_pairs(i, 2);
    p_adjusted_matrix(r1, r2) = p_adjusted(i);
    p_adjusted_matrix(r2, r1) = p_adjusted(i);
end



%% Create separate figure: Chi-square p-values + US metrics bar charts
fprintf('\n=== Creating US Metrics Figure ===\n');

fig_US = figure('Position', [100, 100, 1000, 250], 'Units', 'pixels');
t_US = tiledlayout(fig_US, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

bar_color = [0.6 0.6 0.6];  % Simple grey

%% Column 1: Chi-square p-values matrix
ax_pval = nexttile(t_US, 1);

% Create color matrix with three levels:
% - White for n.s. (p>=0.05)
% - Light grey for 0.001<=p<0.05
% - Dark grey for p<0.001
color_matrix = ones(4, 4, 3);  % Start with white
for r1 = 1:4
    for r2 = 1:4
        if r1 ~= r2
            p_val = p_value_matrix(r1, r2);
            if p_val < 0.001
                color_matrix(r1, r2, :) = [0.50 0.50 0.50];  % Dark grey for p<0.001
            elseif p_val < 0.05
                color_matrix(r1, r2, :) = [0.80 0.80 0.80];  % Light grey for 0.001<=p<0.05
            end
            % Otherwise stays white for n.s. (p>=0.05)
        end
    end
end

imagesc(ax_pval, color_matrix);

% Set axis properties
set(ax_pval, 'XTick', 1:4, 'XTickLabel', region_names_stat);
set(ax_pval, 'YTick', 1:4, 'YTickLabel', region_names_stat);
set(ax_pval, 'FontSize', g.fontSize2);
axis(ax_pval, 'square');
title('Chi-square P-values', 'FontSize', g.fontSize2, 'FontWeight', 'bold');

% Add p-values and significance stars
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

%% Columns 2-4: US response metrics (US-selective and Multisensory neurons only)
clusters_to_use = [2, 3];  % US-selective and Multisensory
metric_names = {'\DeltanSpikes', '\DeltaPeak FR (Hz)', 'Response length (ms)'};
baseline_idx = 1:(g.pre_time / g.bin_time);

% Calculate metrics for all regions
all_metric_data = cell(3, 3);  % 3 metrics × 3 regions

for metric = 1:3
    for br = 1:3
        if isempty(results_all{br})
            all_metric_data{metric, br} = [];
            continue;
        end

        res = results_all{br};

        % Get indices for US-responsive clusters only
        clust_idx = [];
        for c = clusters_to_use
            clust_idx = [clust_idx; find(res.Clusters == c)];
        end

        if isempty(clust_idx)
            all_metric_data{metric, br} = [];
            continue;
        end

        % Calculate metric for each neuron
        metric_values = zeros(length(clust_idx), 1);

        for n = 1:length(clust_idx)
            idx_n = clust_idx(n);

            % Get US PSTH and latencies
            onset_lat = res.US_onset_lat(idx_n);
            offset_lat = res.US_offset_lat(idx_n);

            if ~isnan(onset_lat) && ~isnan(offset_lat)
                onset_bin = round(onset_lat / g.bin_time) + 1 + g.roi(1) - 1;
                offset_bin = round(offset_lat / g.bin_time) + 1 + g.roi(1) - 1;

                if metric == 1  % Delta nSpikes
                    baseline_fr = mean(res.psth_US_Hz(idx_n, baseline_idx));
                    response_fr = mean(res.psth_US_Hz(idx_n, onset_bin:offset_bin));
                    response_duration = (offset_bin - onset_bin + 1) * g.bin_time;
                    metric_values(n) = (response_fr - baseline_fr) * response_duration;

                elseif metric == 2  % Delta Peak FR
                    baseline_fr = mean(res.psth_US_Hz(idx_n, baseline_idx));
                    peak_fr = max(res.psth_US_Hz(idx_n, onset_bin:offset_bin));
                    metric_values(n) = peak_fr - baseline_fr;

                else  % Response length (ms)
                    metric_values(n) = (offset_lat - onset_lat) * 1000;  % Convert to ms
                end
            else
                metric_values(n) = 0;
            end
        end

        all_metric_data{metric, br} = metric_values;
    end
end

% Plot each metric directly in main layout
for metric = 1:3
    ax = nexttile(t_US, metric + 1);
    hold on;

    % Bar positions: 1=LA, 2=BA, 3=Astria
    positions = [1 2 3];

    % Calculate means and SEMs
    means = zeros(3, 1);
    sems = zeros(3, 1);

    for br = 1:3
        if ~isempty(all_metric_data{metric, br})
            means(br) = mean(all_metric_data{metric, br});
            sems(br) = std(all_metric_data{metric, br}) / sqrt(length(all_metric_data{metric, br}));
        end
    end

    % Plot bars with simple grey color
    for br = 1:3
        if means(br) > 0 || sems(br) > 0
            bar(positions(br), means(br), 0.4, 'FaceColor', bar_color, 'EdgeColor', 'k', 'LineWidth', 1);
        end
    end

    % Add error bars
    for br = 1:3
        if means(br) > 0 || sems(br) > 0
            errorbar(positions(br), means(br), sems(br), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
        end
    end

    % Statistical comparisons using Wilcoxon rank-sum test
    y_max = max(means + sems);
    y_range = y_max * 2.0;  % Leave room for significance markers (increased from 1.5)

    % Pairwise comparisons (only for 3 regions: LA, BA, Astria)
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

    % Add significance markers (only significant ones)
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

    % Formatting
    xlim([0.3 3.7]);  % More space on sides
    ylim([0 y_range]);  % Set explicit ylim
    xticks(positions);
    xticklabels({'LA', 'BA', 'AStria'});
    title(metric_names{metric}, 'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Interpreter', 'tex');
    set(gca, 'FontSize', g.fontSize2);
    box off;
end

%% Combined statistical figures (Supplementary)
fprintf('\n=== Creating Combined Statistical Figures ===\n');

fig_stats = figure('Position', [200, 100, 1000, 300]);
t_stats = tiledlayout(fig_stats, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

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

% Panel 2: Chi-square contingency table
ax_table = nexttile(t_stats, 2);
axis off;

% Create table text
table_text = {'Region', 'CS-sel', 'US-sel', 'Multi', 'Non-resp', 'Inhib'};
y_pos = 1.0;
y_step = 0.15;

% Header
for col = 1:6
    text(ax_table, (col-1)*0.16, y_pos, table_text{col}, ...
        'FontSize', g.fontSize2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Data rows
y_pos = y_pos - y_step;
for br = 1:4
    text(ax_table, 0, y_pos, region_names_stat{br}, ...
        'FontSize', g.fontSize2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    for c = 1:5
        text(ax_table, (c)*0.16, y_pos, sprintf('%d', contingency_table(br, c)), ...
            'FontSize', g.fontSize2, 'HorizontalAlignment', 'center');
    end
    y_pos = y_pos - y_step;
end

xlim(ax_table, [-0.1 1.0]);
ylim(ax_table, [0 1.2]);
title(ax_table, sprintf('Contingency Table (\\chi^2=%.2f)', chi2_obs), 'FontSize', g.fontSize1, 'FontWeight', 'bold');

% Panel 3: Cramér's V matrix
ax3 = nexttile(t_stats, 3);
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

% Add values and significance
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

            if cramers_v_matrix(r1, r2) > 0.3
                text_color = 'w';
            else
                text_color = 'k';
            end

            text(ax3, r2, r1, {val_str, sig_str}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', g.fontSize2, 'FontWeight', 'bold', 'Color', text_color);
        else
            text(ax3, r2, r1, '-', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', g.fontSize1, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
end
hold(ax3, 'off');

% Panel 4: Benjamini-Hochberg FDR correction
ax4 = nexttile(t_stats, 4);
hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1);  % Unity line
scatter(pairwise_p_values, p_adjusted, 100, 'filled', 'MarkerFaceColor', [0.3 0.3 0.7]);

% Mark significant/non-significant
sig_mask = p_adjusted < 0.05;
scatter(pairwise_p_values(sig_mask), p_adjusted(sig_mask), 100, 'filled', ...
    'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

plot([0 1], [0.05 0.05], 'r--', 'LineWidth', 1.5);  % FDR threshold
hold off;
xlabel('Uncorrected p-value', 'FontSize', g.fontSize2);
ylabel('FDR-corrected p-value', 'FontSize', g.fontSize2);
title('FDR Correction', 'FontSize', g.fontSize1, 'FontWeight', 'bold');
axis square;
xlim([0 max(pairwise_p_values)*1.1]);
ylim([0 max(p_adjusted)*1.1]);
legend({'Unity', 'Non-sig', 'Significant', 'FDR=0.05'}, 'Location', 'northwest', 'FontSize', g.fontSize2);
set(gca, 'FontSize', g.fontSize2);

fprintf('\nDone.\n');
