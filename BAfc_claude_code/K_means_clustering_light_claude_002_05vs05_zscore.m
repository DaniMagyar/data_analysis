% function BAfc_figure_4_improved
% K-means clustering on z-scored PSTH data
clear all, close all

recordings = {...
    'MD298_001_kilosort',...
    'MD299_001_kilosort',...
    'MD300_001_kilosort',...
    'MD304_001_kilosort',...
    'MD307_001_kilosort',...
    'MD309_001_kilosort',...
    'MD315_001_kilosort',...
    'MD316_002_kilosort',...
    'MD317_001_kilosort',...
    'MD318_001_kilosort',...
    'MD318_002_kilosort',...
    'MD319_003_kilosort'};


 ttl = {'triptest_sound_only', 'triptest_sound_only_light'};
% ttl = {'triptest_shocks_only', 'triptest_shocks_only_light'};
% ttl = {'triptest_both', 'triptest_both_light'};
g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);

g.cell_metrics = BAfc_match_celltypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 13;
g.fontSize2 = 12;
g.xlinewidth = 2;
g.pre_time = 5;
g.post_time = 5;
g.test_time = 0.5;
g.testvalue = 3;
g.bin_time = 0.01;
g.smoothvalue = 5;
g.plotwin = [1 1];
g.spxwin = [0.02 0.1];
g.spxbin = 0.002;
g.timeaxis = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.clim = [-5.5 13];
g.clustnum = 3;
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
g.roi_pca = round((g.pre_time-0.5)/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);
g.bR = 'LA';
g.cluster_data = 3; % 1=left heatmap (no-light), 2=right heatmap (light), 3=both heatmaps

%idx_PN = strcmp(g.cell_metrics.brainRegion,'LA') & strcmp(g.cell_metrics.putativeCellType,'PN');
idx_PN = strcmp(g.cell_metrics.brainRegion, g.bR);

%% Initiate figure
fig = figure('Position', [400, 100, 1400, 800]);
t = tiledlayout(fig,2,3,'TileSpacing', 'compact', 'Padding', 'none');

%% Prepare data and identify light-inhibited neurons
PSTHall = [];

for hmp = 1:size(ttl,2)
    psth_spx_og =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    psth_spx = zscore(psth_spx_og,0,2);
    psth_spx = smoothdata(psth_spx,2,'sgolay', g.smoothvalue);
    psth_spx_og_PN = psth_spx_og(idx_PN,:);
    psth_spx_PN = psth_spx(idx_PN,:);
    hmpdata.(['psth_' num2str(hmp)]) = psth_spx_PN(:,g.roi_pca);
end

% Select data for clustering based on g.cluster_data
switch g.cluster_data
    case 1
        % Use left heatmap (no-light) only
        PSTHall = hmpdata.psth_1;
        fprintf('Clustering based on: No-light condition (left heatmap)\n');
    case 2
        % Use right heatmap (light) only
        PSTHall = hmpdata.psth_2;
        fprintf('Clustering based on: Light condition (right heatmap)\n');
    case 3
        % Use both heatmaps concatenated
        PSTHall = [hmpdata.psth_1, hmpdata.psth_2];
        fprintf('Clustering based on: Both conditions (left + right heatmaps)\n');
    otherwise
        error('g.cluster_data must be 1 (no-light), 2 (light), or 3 (both)');
end

%% Identify light-inhibited neurons
% Use light trial (hmp=2) to find neurons with decreased firing in -0.5 to 0s vs -5 to -0.5s
fprintf('Identifying light-inhibited neurons...\n');
psth_light_raw = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{2}, ...
    'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
psth_light_raw_PN = psth_light_raw(idx_PN,:);

% Define windows: pre-baseline (-5 to -0.5s) and pre-stim (-0.5 to 0s)
pre_baseline_win = round((g.pre_time-5)/g.bin_time+1:(g.pre_time-0.5)/g.bin_time);
pre_stim_win = round((g.pre_time-0.5)/g.bin_time+1:g.pre_time/g.bin_time);

% Calculate mean firing rates in each window
mean_pre_baseline = mean(psth_light_raw_PN(:, pre_baseline_win), 2);
mean_pre_stim = mean(psth_light_raw_PN(:, pre_stim_win), 2);

% Identify light-inhibited neurons: significantly decreased firing
% Using two-sample t-test for each neuron and checking if pre_stim < pre_baseline
light_inhibited_idx = false(size(psth_light_raw_PN, 1), 1);
for i = 1:size(psth_light_raw_PN, 1)
    baseline_data = psth_light_raw_PN(i, pre_baseline_win)';
    prestim_data = psth_light_raw_PN(i, pre_stim_win)';

    % Check if mean firing decreased and if it's statistically significant
    if mean_pre_stim(i) < mean_pre_baseline(i)
        % Use two-sample t-test since window sizes are different
        [~, p_val] = ttest2(baseline_data, prestim_data, 'Tail', 'right');
        if p_val < 0.05
            light_inhibited_idx(i) = true;
        end
    end
end

fprintf('Found %d light-inhibited neurons (%.1f%% of total)\n', ...
    sum(light_inhibited_idx), 100*sum(light_inhibited_idx)/length(light_inhibited_idx));

%% Identify stimulus-inhibited and non-responsive neurons (from non-light-inhibited neurons)
% Use light trials to identify neurons inhibited or non-responsive to stimulus
fprintf('Identifying stimulus-inhibited and non-responsive neurons...\n');
psth_light_stim_raw = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{2}, ...
    'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
psth_light_stim_raw_PN = psth_light_stim_raw(idx_PN,:);

% Define windows: baseline (-0.5 to 0s) and stimulus (0 to 0.5s)
baseline_win = round((g.pre_time-0.5)/g.bin_time+1:g.pre_time/g.bin_time);
stim_win = round(g.pre_time/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);

% Calculate mean firing rates in each window
mean_baseline = mean(psth_light_stim_raw_PN(:, baseline_win), 2);
mean_stim = mean(psth_light_stim_raw_PN(:, stim_win), 2);

% Initialize indices
stim_inhibited_idx = false(size(psth_light_stim_raw_PN, 1), 1);
non_responsive_idx = false(size(psth_light_stim_raw_PN, 1), 1);

% Only test non-light-inhibited neurons
non_light_inhibited = find(~light_inhibited_idx);

for i = 1:length(non_light_inhibited)
    neuron_idx = non_light_inhibited(i);
    baseline_data = psth_light_stim_raw_PN(neuron_idx, baseline_win)';
    stim_data = psth_light_stim_raw_PN(neuron_idx, stim_win)';

    % Two-sample t-test
    [~, p_val] = ttest2(baseline_data, stim_data);

    if p_val < 0.05
        % Significant difference - check direction
        if mean_stim(neuron_idx) < mean_baseline(neuron_idx)
            % Stimulus-inhibited
            stim_inhibited_idx(neuron_idx) = true;
        end
        % If mean_stim > mean_baseline, it's excited (keep for clustering)
    else
        % Not significant - non-responsive
        non_responsive_idx(neuron_idx) = true;
    end
end

fprintf('Found %d stimulus-inhibited neurons (%.1f%% of non-light-inhibited)\n', ...
    sum(stim_inhibited_idx), 100*sum(stim_inhibited_idx)/sum(~light_inhibited_idx));
fprintf('Found %d non-responsive neurons (%.1f%% of non-light-inhibited)\n', ...
    sum(non_responsive_idx), 100*sum(non_responsive_idx)/sum(~light_inhibited_idx));

% Separate data for clustering (only excited neurons)
excited_idx = ~light_inhibited_idx & ~stim_inhibited_idx & ~non_responsive_idx;
fprintf('Using %d excited neurons for clustering (%.1f%% of total)\n', ...
    sum(excited_idx), 100*sum(excited_idx)/length(excited_idx));

PSTHall_for_clustering = PSTHall(excited_idx, :);

%% K-means clustering (only excited neurons)
fprintf('Running k-means clustering with k=%d (only excited neurons)...\n', g.clustnum);

% Remove any rows with NaN from all data first
nan_rows_all = any(isnan(PSTHall), 2);

% Remove any rows with NaN from clustering data
nan_rows_clustering = any(isnan(PSTHall_for_clustering), 2);
PSTHall_clean = PSTHall_for_clustering(~nan_rows_clustering, :);

% Run k-means with multiple replicates to find best solution
[Clusters_clean, Centroids, sumd, D] = kmeans(PSTHall_clean, g.clustnum, ...
    'Distance', 'sqeuclidean', ...
    'Replicates', 20, ...
    'MaxIter', 500, ...
    'Display', 'final');

% Map back to original indexing (all PN neurons)
n_neurons = size(PSTHall, 1);
Clusters = nan(n_neurons, 1);

% Assign cluster numbers to excited neurons
excited_indices = find(excited_idx);
excited_clean_indices = excited_indices(~nan_rows_clustering);
Clusters(excited_clean_indices) = Clusters_clean;

% Assign special groups to separate cluster numbers
Clusters(light_inhibited_idx) = g.clustnum + 1;  % Light-inhibited
Clusters(stim_inhibited_idx) = g.clustnum + 2;    % Stimulus-inhibited
Clusters(non_responsive_idx) = g.clustnum + 3;    % Non-responsive

% For visualization, create hierarchical ordering based on k-means results
% Use cluster centroids to determine dendrogram
centroid_dist = pdist(Centroids, 'euclidean');
Dend_centroids = linkage(centroid_dist, 'average');

% Within each cluster, order neurons by similarity to centroid
leafOrder = [];
for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    if ~isempty(clust_idx)
        % Get distances to centroid for this cluster
        dists_to_cent = pdist2(PSTHall(clust_idx, :), Centroids(c, :), 'euclidean');
        [~, sort_idx] = sort(dists_to_cent);
        leafOrder = [leafOrder; clust_idx(sort_idx)];
    end
end

% Get no-light PSTH for sorting special groups
psth_nolight_sort = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{1}, ...
    'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
psth_nolight_sort_PN = psth_nolight_sort(idx_PN,:);
post_stim_win = round(g.pre_time/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);

% Add light-inhibited neurons
light_inhib_idx = find(Clusters == g.clustnum + 1);
if ~isempty(light_inhib_idx)
    % Sort by post-stimulus activity (descending = most active on top)
    post_stim_activity = mean(psth_nolight_sort_PN(light_inhib_idx, post_stim_win), 2);
    [~, sort_idx] = sort(post_stim_activity, 'descend');
    leafOrder = [leafOrder; light_inhib_idx(sort_idx)];
end

% Add stimulus-inhibited neurons
stim_inhib_idx = find(Clusters == g.clustnum + 2);
if ~isempty(stim_inhib_idx)
    % Sort by inhibition strength (most inhibited on top)
    inhibition_strength = mean_baseline(stim_inhib_idx) - mean_stim(stim_inhib_idx);
    [~, sort_idx] = sort(inhibition_strength, 'descend');
    leafOrder = [leafOrder; stim_inhib_idx(sort_idx)];
end

% Add non-responsive neurons
non_resp_idx = find(Clusters == g.clustnum + 3);
if ~isempty(non_resp_idx)
    % Sort by baseline firing rate
    baseline_activity = mean(psth_nolight_sort_PN(non_resp_idx, baseline_win), 2);
    [~, sort_idx] = sort(baseline_activity, 'descend');
    leafOrder = [leafOrder; non_resp_idx(sort_idx)];
end

leafOrderfl = leafOrder;

%% Plot heatmaps
mploc = [1 2];
for hmp = 1:size(ttl,2)
    ax = nexttile(t,mploc(hmp),[2 1]);
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    psth_spx = zscore(psth_spx,0,2);   
    psth_spx = smoothdata(psth_spx,2,'sgolay', g.smoothvalue);
    psth_spx_PN = psth_spx(idx_PN,:);
    psth_spx_sorted = [psth_spx_PN(leafOrderfl,:)];
    matrix = psth_spx_sorted(:,(g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    imagesc(g.timeaxis,1:size(matrix,1),matrix);
    clim(g.clim)
    colormap(g.colors.Heatmap); 
    yline(sum(idx_PN)+0.5, 'k', 'LineWidth', 1);
    xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
    ylabel('Cell number')
    set(gca, 'FontSize', g.fontSize2);
    title(ttl{hmp}, 'FontSize', g.fontSize1)
    if hmp == 1
        cb = colorbar('westoutside', 'FontSize', g.fontSize2);
    end
    hold on
    n_clu = find(diff(Clusters(leafOrderfl))~=0);
    yline(n_clu+0.5, 'Color', 'k', 'LineWidth', 1);
    hold off
end

%% Plot response magnitude boxplots (Peak and AUC) - separate plot per cluster
% Get raw PSTH data for both conditions
psth_nolight_raw = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{1}, ...
    'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
psth_light_raw = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{2}, ...
    'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

psth_nolight_raw_PN = psth_nolight_raw(idx_PN,:);
psth_light_raw_PN = psth_light_raw(idx_PN,:);

% Define windows
baseline_win = round((g.pre_time-0.5)/g.bin_time+1:g.pre_time/g.bin_time);
stim_win = round(g.pre_time/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);

% Prepare data structure for plotting
total_clusters = g.clustnum + 3;
colors_cond = [0 0 0; 0 0 1]; % black for no-light, blue for light

% Create a new figure for cluster-specific boxplots
fig_boxplots = figure('Position', [1450, 100, 400, 800]);
t_box = tiledlayout(fig_boxplots, total_clusters, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for c = 1:total_clusters
    clust_idx = find(Clusters == c);

    if ~isempty(clust_idx)
        % No-light condition
        baseline_nolight = psth_nolight_raw_PN(clust_idx, baseline_win);
        stim_nolight = psth_nolight_raw_PN(clust_idx, stim_win);

        % Peak response (max - baseline mean)
        peak_nolight = max(stim_nolight, [], 2) - mean(baseline_nolight, 2);

        % AUC (area under curve: sum of response - baseline)
        auc_nolight = sum(stim_nolight, 2) - sum(baseline_nolight, 2);

        % Light condition
        baseline_light = psth_light_raw_PN(clust_idx, baseline_win);
        stim_light = psth_light_raw_PN(clust_idx, stim_win);

        peak_light = max(stim_light, [], 2) - mean(baseline_light, 2);
        auc_light = sum(stim_light, 2) - sum(baseline_light, 2);

        % Plot for this cluster
        ax1 = nexttile(t_box, c);

        % Plot Peak on left y-axis
        yyaxis left
        hold on

        % Peak - no-light
        boxplot(peak_nolight, 'Positions', 1, 'Width', 0.3, ...
            'Colors', colors_cond(1,:), 'Symbol', '');
        scatter(ones(length(peak_nolight), 1), peak_nolight, ...
            15, colors_cond(1,:), 'filled', 'MarkerFaceAlpha', 0.4);

        % Peak - light
        boxplot(peak_light, 'Positions', 1.4, 'Width', 0.3, ...
            'Colors', colors_cond(2,:), 'Symbol', '');
        scatter(ones(length(peak_light), 1)*1.4, peak_light, ...
            15, colors_cond(2,:), 'filled', 'MarkerFaceAlpha', 0.4);

        % Connect paired data points for Peak
        for n = 1:length(peak_nolight)
            plot([1, 1.4], [peak_nolight(n), peak_light(n)], ...
                '-', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 0.5);
        end

        ylabel('Peak (Hz)')
        ax1.YColor = 'k';

        % Statistical test for Peak (paired Wilcoxon signed-rank test)
        p_peak = signrank(peak_nolight, peak_light);

        % Add significance marker for Peak
        y_max_peak = max([peak_nolight; peak_light]);
        y_min_peak = min([peak_nolight; peak_light]);
        y_range_peak = y_max_peak - y_min_peak;
        sig_y_peak = y_max_peak + 0.1*y_range_peak;

        if p_peak < 0.001
            text(1.2, sig_y_peak, '***', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        elseif p_peak < 0.01
            text(1.2, sig_y_peak, '**', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        elseif p_peak < 0.05
            text(1.2, sig_y_peak, '*', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        else
            text(1.2, sig_y_peak, 'n.s.', 'HorizontalAlignment', 'center', 'FontSize', 10);
        end

        % Draw line for significance
        plot([1, 1.4], [sig_y_peak*0.98, sig_y_peak*0.98], 'k-', 'LineWidth', 1);

        hold off

        % Plot AUC on right y-axis
        yyaxis right
        hold on

        % AUC - no-light
        boxplot(auc_nolight, 'Positions', 2.2, 'Width', 0.3, ...
            'Colors', colors_cond(1,:), 'Symbol', '');
        scatter(ones(length(auc_nolight), 1)*2.2, auc_nolight, ...
            15, colors_cond(1,:), 'filled', 'MarkerFaceAlpha', 0.4);

        % AUC - light
        boxplot(auc_light, 'Positions', 2.6, 'Width', 0.3, ...
            'Colors', colors_cond(2,:), 'Symbol', '');
        scatter(ones(length(auc_light), 1)*2.6, auc_light, ...
            15, colors_cond(2,:), 'filled', 'MarkerFaceAlpha', 0.4);

        % Connect paired data points for AUC
        for n = 1:length(auc_nolight)
            plot([2.2, 2.6], [auc_nolight(n), auc_light(n)], ...
                '-', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 0.5);
        end

        ylabel('AUC (Hz × time)')
        ax1.YColor = 'k';

        % Statistical test for AUC (paired Wilcoxon signed-rank test)
        p_auc = signrank(auc_nolight, auc_light);

        % Add significance marker for AUC
        y_max_auc = max([auc_nolight; auc_light]);
        y_min_auc = min([auc_nolight; auc_light]);
        y_range_auc = y_max_auc - y_min_auc;
        sig_y_auc = y_max_auc + 0.1*y_range_auc;

        if p_auc < 0.001
            text(2.4, sig_y_auc, '***', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        elseif p_auc < 0.01
            text(2.4, sig_y_auc, '**', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        elseif p_auc < 0.05
            text(2.4, sig_y_auc, '*', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        else
            text(2.4, sig_y_auc, 'n.s.', 'HorizontalAlignment', 'center', 'FontSize', 10);
        end

        % Draw line for significance
        plot([2.2, 2.6], [sig_y_auc*0.98, sig_y_auc*0.98], 'k-', 'LineWidth', 1);

        hold off

        % Formatting
        xlim([0.5 3.1])
        xticks([1.2 2.4])
        xticklabels({'Peak', 'AUC'})

        % Create cluster label
        if c <= g.clustnum
            title(sprintf('Cluster %d (n=%d)', c, length(clust_idx)), 'FontSize', g.fontSize2);
        elseif c == g.clustnum + 1
            title(sprintf('Light-Inhibited (n=%d)', length(clust_idx)), 'FontSize', g.fontSize2);
        elseif c == g.clustnum + 2
            title(sprintf('Stim-Inhibited (n=%d)', length(clust_idx)), 'FontSize', g.fontSize2);
        else
            title(sprintf('Non-Responsive (n=%d)', length(clust_idx)), 'FontSize', g.fontSize2);
        end

        set(gca, 'FontSize', g.fontSize2-2)
        grid on
    end
end

%% Display results
fprintf('\n=== K-means Clustering Results ===\n');
fprintf('Total sum of distances: %.2f\n', sum(sumd));
fprintf('Cluster sizes: ');
for c = 1:g.clustnum
    fprintf('C%d:%d ', c, sum(Clusters == c));
end
fprintf('Light-inhib:%d ', sum(Clusters == g.clustnum + 1));
fprintf('Stim-inhib:%d ', sum(Clusters == g.clustnum + 2));
fprintf('Non-resp:%d ', sum(Clusters == g.clustnum + 3));
fprintf('\n');

% Silhouette analysis
if ~any(nan_rows_clustering)
    silh_vals = silhouette(PSTHall_clean, Clusters_clean);
    fprintf('Mean Silhouette Score: %.3f\n', mean(silh_vals));
end

%% ===== CLUSTER-SPECIFIC METRICS CALCULATION (Z-SCORE BASED) =====

fprintf('\n========== CLUSTER-SPECIFIC ANALYSIS (Z-SCORE) ==========\n');

% Initialize storage for cluster metrics
cluster_metrics = struct();

% Include all clusters (excited + 3 special groups)
total_clusters = g.clustnum + 3;

for c = 1:total_clusters
    % Get neurons in this cluster
    clust_idx = find(Clusters == c);
    n_neurons = length(clust_idx);

    if c <= g.clustnum
        fprintf('\n--- Cluster %d (%d neurons) ---\n', c, n_neurons);
    elseif c == g.clustnum + 1
        fprintf('\n--- Light-Inhibited Cluster (%d neurons) ---\n', n_neurons);
    elseif c == g.clustnum + 2
        fprintf('\n--- Stimulus-Inhibited Cluster (%d neurons) ---\n', n_neurons);
    else
        fprintf('\n--- Non-Responsive Cluster (%d neurons) ---\n', n_neurons);
    end

    % Extract PSTH data for each condition
    for hmp = 1:size(ttl,2)
        % Get raw PSTH
        psth_spx_og = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, ...
            'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
        psth_cluster_raw = psth_spx_og(idx_PN,:);
        psth_cluster_raw = psth_cluster_raw(clust_idx,:);

        % Z-score each neuron relative to its own baseline (-0.5 to 0s)
        baseline_window = round((g.pre_time-0.5)/g.bin_time+1:g.pre_time/g.bin_time);
        baseline_mean = mean(psth_cluster_raw(:, baseline_window), 2);
        baseline_std = std(psth_cluster_raw(:, baseline_window), 0, 2);

        % Compute z-score for entire trace
        psth_cluster_zscore = (psth_cluster_raw - baseline_mean) ./ baseline_std;

        % Store data
        cluster_metrics(c).(['psth_raw_' num2str(hmp)]) = psth_cluster_raw;
        cluster_metrics(c).(['psth_zscore_' num2str(hmp)]) = psth_cluster_zscore;

        % Calculate mean response across neurons
        cluster_metrics(c).(['mean_psth_raw_' num2str(hmp)]) = mean(psth_cluster_raw, 1, 'omitnan');
        cluster_metrics(c).(['sem_psth_raw_' num2str(hmp)]) = std(psth_cluster_raw, 0, 1, 'omitnan') / sqrt(n_neurons);
        cluster_metrics(c).(['mean_psth_zscore_' num2str(hmp)]) = mean(psth_cluster_zscore, 1, 'omitnan');
        cluster_metrics(c).(['sem_psth_zscore_' num2str(hmp)]) = std(psth_cluster_zscore, 0, 1, 'omitnan') / sqrt(n_neurons);

        % Calculate z-score response in stimulus window (0 to 0.5s)
        stim_window = round(g.pre_time/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);
        baseline_window_resp = round((g.pre_time-0.5)/g.bin_time+1:g.pre_time/g.bin_time);

        % Mean z-score in each window
        stim_zscore = mean(psth_cluster_zscore(:, stim_window), 2);
        baseline_zscore = mean(psth_cluster_zscore(:, baseline_window_resp), 2);

        cluster_metrics(c).(['stim_zscore_' num2str(hmp)]) = stim_zscore;
        cluster_metrics(c).(['baseline_zscore_' num2str(hmp)]) = baseline_zscore;
        cluster_metrics(c).(['mean_stim_zscore_' num2str(hmp)]) = mean(stim_zscore, 'omitnan');
        cluster_metrics(c).(['mean_baseline_zscore_' num2str(hmp)]) = mean(baseline_zscore, 'omitnan');
    end

    % Calculate difference metrics (Light - No Light) in z-score space
    cluster_metrics(c).delta_stim_zscore = cluster_metrics(c).stim_zscore_2 - cluster_metrics(c).stim_zscore_1;
    cluster_metrics(c).mean_delta_zscore = mean(cluster_metrics(c).delta_stim_zscore, 'omitnan');
    cluster_metrics(c).sem_delta_zscore = std(cluster_metrics(c).delta_stim_zscore, 0, 'omitnan') / sqrt(n_neurons);

    % Statistical test (paired t-test)
    [~, cluster_metrics(c).p_val, ~, stats] = ttest(cluster_metrics(c).stim_zscore_1, ...
                                                      cluster_metrics(c).stim_zscore_2);
    cluster_metrics(c).t_stat = stats.tstat;
    cluster_metrics(c).effect_size = mean(cluster_metrics(c).delta_stim_zscore) / std(cluster_metrics(c).delta_stim_zscore);

    fprintf('Mean stim z-score (no-light): %.3f\n', cluster_metrics(c).mean_stim_zscore_1);
    fprintf('Mean stim z-score (light): %.3f\n', cluster_metrics(c).mean_stim_zscore_2);
    fprintf('Mean z-score change: %.3f (p=%.4f, d=%.3f)\n', ...
        cluster_metrics(c).mean_delta_zscore, cluster_metrics(c).p_val, cluster_metrics(c).effect_size);

    cluster_metrics(c).n_neurons = n_neurons;
end

%% ===== COMPREHENSIVE VISUALIZATION =====

fig2 = figure('Position', [50, 50, 1800, 1200]);
t2 = tiledlayout(fig2, 4, total_clusters, 'TileSpacing', 'compact', 'Padding', 'compact');

for c = 1:total_clusters

    % Row 1: Average Z-scored PSTH comparison (no-light vs light)
    nexttile(t2, c);

    % Create time vector matching PSTH length
    n_bins = length(cluster_metrics(c).mean_psth_zscore_1);
    time_vec_full = linspace(-g.pre_time, g.post_time, n_bins);

    plot(time_vec_full, cluster_metrics(c).mean_psth_zscore_1, 'k', 'LineWidth', 2); hold on;
    plot(time_vec_full, cluster_metrics(c).mean_psth_zscore_2, 'b', 'LineWidth', 2);

    % Add SEM shading
    fill([time_vec_full, fliplr(time_vec_full)], ...
         [cluster_metrics(c).mean_psth_zscore_1 + cluster_metrics(c).sem_psth_zscore_1, ...
          fliplr(cluster_metrics(c).mean_psth_zscore_1 - cluster_metrics(c).sem_psth_zscore_1)], ...
         'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([time_vec_full, fliplr(time_vec_full)], ...
         [cluster_metrics(c).mean_psth_zscore_2 + cluster_metrics(c).sem_psth_zscore_2, ...
          fliplr(cluster_metrics(c).mean_psth_zscore_2 - cluster_metrics(c).sem_psth_zscore_2)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    xline(0, '--r', 'LineWidth', 1.5);
    xline(-0.5, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); % baseline window start
    xlabel('Time (s)'); ylabel('Z-score');
    if c <= g.clustnum
        title(sprintf('C%d: Avg Z-score Response', c));
    elseif c == g.clustnum + 1
        title('Light-Inhib: Avg Z-score Response');
    elseif c == g.clustnum + 2
        title('Stim-Inhib: Avg Z-score Response');
    else
        title('Non-Resp: Avg Z-score Response');
    end
    legend({'No Light', 'Light'}, 'Location', 'best');
    set(gca, 'FontSize', 10);
    xlim([-g.plotwin(1) g.plotwin(2)]);

    % Row 2: Heatmap of individual neurons (sorted by z-score change)
    nexttile(t2, total_clusters + c);
    [~, sort_idx] = sort(cluster_metrics(c).delta_stim_zscore, 'descend');
    psth_sorted = cluster_metrics(c).psth_zscore_1(sort_idx, :);
    imagesc(time_vec_full, 1:cluster_metrics(c).n_neurons, psth_sorted);
    colormap(gca, g.colors.Heatmap);
    clim([-3 5]);
    xline(0, '--w', 'LineWidth', 1.5);
    xline(-0.5, ':', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Neuron #');
    if c <= g.clustnum
        title(sprintf('C%d: No-Light Z-scored', c));
    elseif c == g.clustnum + 1
        title('Light-Inhib: No-Light Z-scored');
    elseif c == g.clustnum + 2
        title('Stim-Inhib: No-Light Z-scored');
    else
        title('Non-Resp: No-Light Z-scored');
    end
    set(gca, 'FontSize', 10);
    xlim([-g.plotwin(1) g.plotwin(2)]);

    % Row 3: Delta Z-score PSTH (Light - No Light) heatmap
    nexttile(t2, 2*total_clusters + c);
    delta_psth_z = cluster_metrics(c).psth_zscore_2(sort_idx,:) - cluster_metrics(c).psth_zscore_1(sort_idx,:);
    imagesc(time_vec_full, 1:cluster_metrics(c).n_neurons, delta_psth_z);

    % Create custom red-white-blue colormap
    n = 128;
    red_to_white = [linspace(0.5,1,n/2)', linspace(0,1,n/2)', linspace(0,1,n/2)'];
    white_to_blue = [linspace(1,0,n/2)', linspace(1,0,n/2)', linspace(1,1,n/2)'];
    redblue_map = [red_to_white; white_to_blue];
    colormap(gca, redblue_map);
    clim([-3 3]);
    xline(0, '--k', 'LineWidth', 1.5);
    xline(-0.5, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Neuron #');
    if c <= g.clustnum
        title(sprintf('C%d: Δ Z-score (Light-NoLight)', c));
    elseif c == g.clustnum + 1
        title('Light-Inhib: Δ Z-score (Light-NoLight)');
    elseif c == g.clustnum + 2
        title('Stim-Inhib: Δ Z-score (Light-NoLight)');
    else
        title('Non-Resp: Δ Z-score (Light-NoLight)');
    end
    set(gca, 'FontSize', 10);
    xlim([-g.plotwin(1) g.plotwin(2)]);

    % Row 4: Single-neuron scatter plot (z-scores)
    nexttile(t2, 3*total_clusters + c);
    scatter(cluster_metrics(c).stim_zscore_1, cluster_metrics(c).stim_zscore_2, ...
            30, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    max_val = max(abs([cluster_metrics(c).stim_zscore_1; cluster_metrics(c).stim_zscore_2]));
    min_val = min([cluster_metrics(c).stim_zscore_1; cluster_metrics(c).stim_zscore_2]);
    plot([min_val max_val], [min_val max_val], '--k', 'LineWidth', 1.5);
    xlabel('No-Light (z-score)'); ylabel('Light (z-score)');
    if c <= g.clustnum
        title(sprintf('C%d: p=%.4f', c, cluster_metrics(c).p_val));
    elseif c == g.clustnum + 1
        title(sprintf('Light-Inhib: p=%.4f', cluster_metrics(c).p_val));
    elseif c == g.clustnum + 2
        title(sprintf('Stim-Inhib: p=%.4f', cluster_metrics(c).p_val));
    else
        title(sprintf('Non-Resp: p=%.4f', cluster_metrics(c).p_val));
    end
    axis square;
    set(gca, 'FontSize', 10);
end

%% ===== COMPREHENSIVE VISUALIZATION - FIRING RATE (Hz) VERSION =====

fig3 = figure('Position', [100, 100, 1800, 1200]);
t3 = tiledlayout(fig3, 4, total_clusters, 'TileSpacing', 'compact', 'Padding', 'compact');

for c = 1:total_clusters

    % Row 1: Average Firing Rate PSTH comparison (no-light vs light)
    nexttile(t3, c);

    % Create time vector matching PSTH length
    n_bins = length(cluster_metrics(c).mean_psth_raw_1);
    time_vec_full = linspace(-g.pre_time, g.post_time, n_bins);

    plot(time_vec_full, cluster_metrics(c).mean_psth_raw_1, 'k', 'LineWidth', 2); hold on;
    plot(time_vec_full, cluster_metrics(c).mean_psth_raw_2, 'b', 'LineWidth', 2);

    % Add SEM shading
    fill([time_vec_full, fliplr(time_vec_full)], ...
         [cluster_metrics(c).mean_psth_raw_1 + cluster_metrics(c).sem_psth_raw_1, ...
          fliplr(cluster_metrics(c).mean_psth_raw_1 - cluster_metrics(c).sem_psth_raw_1)], ...
         'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([time_vec_full, fliplr(time_vec_full)], ...
         [cluster_metrics(c).mean_psth_raw_2 + cluster_metrics(c).sem_psth_raw_2, ...
          fliplr(cluster_metrics(c).mean_psth_raw_2 - cluster_metrics(c).sem_psth_raw_2)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    xline(0, '--r', 'LineWidth', 1.5);
    xline(-0.5, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); % baseline window start
    xlabel('Time (s)'); ylabel('Firing Rate (Hz)');
    if c <= g.clustnum
        title(sprintf('C%d: Avg Firing Rate', c));
    elseif c == g.clustnum + 1
        title('Light-Inhib: Avg Firing Rate');
    elseif c == g.clustnum + 2
        title('Stim-Inhib: Avg Firing Rate');
    else
        title('Non-Resp: Avg Firing Rate');
    end
    legend({'No Light', 'Light'}, 'Location', 'best');
    set(gca, 'FontSize', 10);
    xlim([-g.plotwin(1) g.plotwin(2)]);

    % Row 2: Heatmap of individual neurons (sorted by raw FR change)
    nexttile(t3, total_clusters + c);

    % Calculate firing rate change in stimulus window for sorting
    stim_window = round(g.pre_time/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);
    stim_fr_1 = mean(cluster_metrics(c).psth_raw_1(:, stim_window), 2);
    stim_fr_2 = mean(cluster_metrics(c).psth_raw_2(:, stim_window), 2);
    delta_stim_fr = stim_fr_2 - stim_fr_1;

    [~, sort_idx] = sort(delta_stim_fr, 'descend');
    psth_sorted = cluster_metrics(c).psth_raw_1(sort_idx, :);
    imagesc(time_vec_full, 1:cluster_metrics(c).n_neurons, psth_sorted);
    colormap(gca, g.colors.Heatmap);
    % Auto-scale color limits based on data
    clim([0 prctile(psth_sorted(:), 95)]);
    xline(0, '--w', 'LineWidth', 1.5);
    xline(-0.5, ':', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Neuron #');
    if c <= g.clustnum
        title(sprintf('C%d: No-Light FR (Hz)', c));
    elseif c == g.clustnum + 1
        title('Light-Inhib: No-Light FR (Hz)');
    elseif c == g.clustnum + 2
        title('Stim-Inhib: No-Light FR (Hz)');
    else
        title('Non-Resp: No-Light FR (Hz)');
    end
    set(gca, 'FontSize', 10);
    xlim([-g.plotwin(1) g.plotwin(2)]);

    % Row 3: Delta Firing Rate PSTH (Light - No Light) heatmap
    nexttile(t3, 2*total_clusters + c);
    delta_psth_fr = cluster_metrics(c).psth_raw_2(sort_idx,:) - cluster_metrics(c).psth_raw_1(sort_idx,:);
    imagesc(time_vec_full, 1:cluster_metrics(c).n_neurons, delta_psth_fr);

    % Create custom red-white-blue colormap
    n = 128;
    red_to_white = [linspace(0.5,1,n/2)', linspace(0,1,n/2)', linspace(0,1,n/2)'];
    white_to_blue = [linspace(1,0,n/2)', linspace(1,0,n/2)', linspace(1,1,n/2)'];
    redblue_map = [red_to_white; white_to_blue];
    colormap(gca, redblue_map);

    % Symmetric color limits around zero
    max_delta = max(abs(delta_psth_fr(:)));
    clim([-max_delta max_delta]);

    xline(0, '--k', 'LineWidth', 1.5);
    xline(-0.5, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Neuron #');
    if c <= g.clustnum
        title(sprintf('C%d: Δ FR (Light-NoLight, Hz)', c));
    elseif c == g.clustnum + 1
        title('Light-Inhib: Δ FR (Light-NoLight, Hz)');
    elseif c == g.clustnum + 2
        title('Stim-Inhib: Δ FR (Light-NoLight, Hz)');
    else
        title('Non-Resp: Δ FR (Light-NoLight, Hz)');
    end
    set(gca, 'FontSize', 10);
    xlim([-g.plotwin(1) g.plotwin(2)]);

    % Row 4: Single-neuron scatter plot (firing rates in Hz)
    nexttile(t3, 3*total_clusters + c);
    scatter(stim_fr_1, stim_fr_2, 30, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    max_val = max([stim_fr_1; stim_fr_2]);
    min_val = min([stim_fr_1; stim_fr_2]);
    plot([min_val max_val], [min_val max_val], '--k', 'LineWidth', 1.5);
    xlabel('No-Light (Hz)'); ylabel('Light (Hz)');

    % Calculate p-value for firing rate change
    [~, p_val_fr] = ttest(stim_fr_1, stim_fr_2);

    if c <= g.clustnum
        title(sprintf('C%d: p=%.4f', c, p_val_fr));
    elseif c == g.clustnum + 1
        title(sprintf('Light-Inhib: p=%.4f', p_val_fr));
    elseif c == g.clustnum + 2
        title(sprintf('Stim-Inhib: p=%.4f', p_val_fr));
    else
        title(sprintf('Non-Resp: p=%.4f', p_val_fr));
    end
    axis square;
    set(gca, 'FontSize', 10);
end

%% ===== POOLED PSTH OF ALL EXCITED NEURONS =====

fprintf('\n========== POOLED PSTH OF ALL EXCITED NEURONS ==========\n');

% Get indices of all excited neurons (clusters 1 through g.clustnum)
all_excited_idx = [];
for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    all_excited_idx = [all_excited_idx; clust_idx];
end

fprintf('Total excited neurons: %d\n', length(all_excited_idx));

% Get global neuron IDs (mapping from idx_PN to global cell_metrics)
global_neuron_ids = find(idx_PN);
excited_global_ids = global_neuron_ids(all_excited_idx);

% Pool spikes from all excited neurons for each stimulus condition
fig_pooled = figure('Position', [150, 150, 1800, 600]);
t_pooled = tiledlayout(fig_pooled, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% Define time bins for pooled PSTH
time_bin_edges_pooled = -g.pre_time:g.bin_time:g.post_time;
n_bins_pooled = length(time_bin_edges_pooled) - 1;
time_centers_pooled = time_bin_edges_pooled(1:end-1) + g.bin_time/2;

% Store pooled PSTHs for overlay plot
pooled_psth_smooth = cell(1, size(ttl, 2));
pooled_psth_raw = cell(1, size(ttl, 2));

for hmp = 1:size(ttl, 2)
    % Initialize pooled spike times
    all_pooled_spikes = [];
    total_trials = 0;

    % Loop through all excited neurons
    for i = 1:length(excited_global_ids)
        neuron_id = excited_global_ids(i);

        % Get spike times for this neuron
        spike_times = g.cell_metrics.spikes.times{neuron_id};

        % Get stimulus times for this neuron
        stim_times = g.cell_metrics.general.(ttl{hmp}){neuron_id};

        % For each trial, collect spikes relative to stimulus onset
        for trial = 1:length(stim_times)
            trial_spikes = spike_times(spike_times >= (stim_times(trial) - g.pre_time) & ...
                                       spike_times < (stim_times(trial) + g.post_time));
            trial_spikes_rel = trial_spikes - stim_times(trial);
            all_pooled_spikes = [all_pooled_spikes; trial_spikes_rel];
            total_trials = total_trials + 1;
        end
    end

    % Create PSTH from pooled spikes
    spike_counts_pooled = histcounts(all_pooled_spikes, time_bin_edges_pooled);
    firing_rate_pooled = spike_counts_pooled / (total_trials * g.bin_time);  % Hz per trial

    % Smooth the pooled PSTH
    firing_rate_pooled_smooth = smoothdata(firing_rate_pooled, 'sgolay', g.smoothvalue);

    % Store for overlay plot
    pooled_psth_smooth{hmp} = firing_rate_pooled_smooth;
    pooled_psth_raw{hmp} = firing_rate_pooled;

    fprintf('Stimulus %d (%s): Pooled %d trials from %d neurons\n', ...
        hmp, ttl{hmp}, total_trials, length(excited_global_ids));

    % Plot pooled PSTH
    nexttile(t_pooled, hmp);

    % Plot firing rate
    plot(time_centers_pooled, firing_rate_pooled_smooth, 'k', 'LineWidth', 2.5);
    hold on;

    % Add raw data as thin line for reference
    plot(time_centers_pooled, firing_rate_pooled, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);

    % Add stimulus onset line
    xline(0, '--r', 'LineWidth', 2, 'Alpha', 0.8);

    % Add baseline reference
    baseline_bins_pooled = time_centers_pooled < 0;
    baseline_mean_pooled = mean(firing_rate_pooled_smooth(baseline_bins_pooled));
    yline(baseline_mean_pooled, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

    % Formatting
    xlabel('Time from Stimulus Onset (s)', 'FontSize', g.fontSize1);
    ylabel('Firing Rate (Hz)', 'FontSize', g.fontSize1);
    title(sprintf('%s\n(%d neurons, %d trials)', ttl{hmp}, ...
        length(excited_global_ids), total_trials), 'FontSize', g.fontSize1);

    xlim([-g.plotwin(1) g.plotwin(2)]);
    grid on;
    set(gca, 'FontSize', g.fontSize2);

    legend({'Smoothed', 'Raw', 'Stimulus Onset', 'Baseline Mean'}, ...
        'Location', 'best', 'FontSize', g.fontSize2-2);
end

%% Panel 3: Overlay both conditions with different colors
nexttile(t_pooled, 3);

% Define colors for each condition
colors_stim = [0 0 0; 0 0.4 0.8]; % Black for no-light, Blue for light

% Plot both conditions
hold on;
h1 = plot(time_centers_pooled, pooled_psth_smooth{1}, 'Color', colors_stim(1,:), ...
    'LineWidth', 2.5, 'DisplayName', strrep(ttl{1}, '_', ' '));
h2 = plot(time_centers_pooled, pooled_psth_smooth{2}, 'Color', colors_stim(2,:), ...
    'LineWidth', 2.5, 'DisplayName', strrep(ttl{2}, '_', ' '));

% Add shaded error regions or raw data as thin lines
plot(time_centers_pooled, pooled_psth_raw{1}, 'Color', [colors_stim(1,:) 0.3], ...
    'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(time_centers_pooled, pooled_psth_raw{2}, 'Color', [colors_stim(2,:) 0.3], ...
    'LineWidth', 0.5, 'HandleVisibility', 'off');

% Add stimulus onset line
xline(0, '--r', 'LineWidth', 2, 'Alpha', 0.8, 'DisplayName', 'Stimulus Onset');

% Add baseline reference (using no-light baseline)
baseline_bins_pooled = time_centers_pooled < 0;
baseline_mean_pooled = mean(pooled_psth_smooth{1}(baseline_bins_pooled));
yline(baseline_mean_pooled, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
    'DisplayName', 'Baseline');

% Formatting
xlabel('Time from Stimulus Onset (s)', 'FontSize', g.fontSize1);
ylabel('Firing Rate (Hz)', 'FontSize', g.fontSize1);
title(sprintf('Overlay Comparison\n(%d neurons)', length(excited_global_ids)), ...
    'FontSize', g.fontSize1);

xlim([-g.plotwin(1) g.plotwin(2)]);
grid on;
set(gca, 'FontSize', g.fontSize2);

legend('show', 'Location', 'best', 'FontSize', g.fontSize2-2);
hold off;

% Overall title
t_pooled.Title.String = sprintf('Pooled PSTH: All Excited Neurons (%s, n=%d)', ...
    g.bR, length(excited_global_ids));
t_pooled.Title.FontSize = g.fontSize1 + 2;

fprintf('Pooled PSTH figure created.\n');
fprintf('========================================\n\n');


