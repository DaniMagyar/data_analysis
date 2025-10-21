% function BAfc_clustering_analysis
% K-means clustering on z-scored PSTH data - comparing responses within clusters
%
% This script:
% 1. Prepares Z-scored PSTH data for clustering (LA neurons, +/- 0.5s around stimulus onset).
% 2. Performs K-means clustering (k=5 by default).
% 3. Calculates response onset latency to determine sorting order within clusters.
% 4. Generates a figure comparing the mean cluster responses across the three stimuli (CS, US, CS+US)
% 5. **NEW:** Generates a separate figure showing the mathematical justification for the choice of k (Elbow and Silhouette methods).
clear all; close all
% 
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
    'MD307_001_kilosort',...
    'MD309_001_kilosort',...
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

ttl = {'triptest_sound_only','triptest_shocks_only', 'triptest_both'};

hmptitles = {'CS', 'US', 'CS + US'}; % Titles for plots

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

% Clustering Parameters
g.clustnum = 5;                 % The chosen number of clusters (k)
g.k_max = 10;                   % Maximum k to test for optimal k determination

% Onset Detection Parameters
g.testvalue = 3;                % significance threshold (z-score)
g.onset_threshold = g.testvalue; % z-score threshold for onset
g.bin_time = 0.01;
g.min_consec_bins = max(1, round(0.02 / g.bin_time)); % minimum consecutive bins (20 ms)

g.smoothvalue = 5;
g.plotwin = [0.5 1];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.clim = [-5.5 13];

% ROI for k-means/PCA (time axis for cluster comparison plot: +/- 0.5 s)
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time; % 0 to 0.5s from stimulus onset (for latency check)
g.roi_pca_idx = round((g.pre_time-0.5)/g.bin_time:(g.pre_time+0.5)/g.bin_time); % -0.5s to 0.5s (for clustering)
g.timeaxis_pca = -0.5:g.bin_time:0.5;

g.bR = 'LA';
idx_PN = strcmp(g.cell_metrics.brainRegion, g.bR) & strcmp(g.cell_metrics.putativeCellType, 'PN');
% idx_PN = strcmp(g.cell_metrics.brainRegion, g.bR);
%% Prepare data

PSTHall = [];
PSTH_ROI_per_stim = cell(1, numel(ttl)); % Store PSTH data in ROI for each stimulus
psthZ_full = cell(1, numel(ttl)); % Store full PSTH for latency detection

for hmp = 1:numel(ttl)
    % 1. Get raw PSTH
    psth_spx_og =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    % 2. Z-score and smooth
    psth_spx = zscore(psth_spx_og,0,2);
    psth_spx = smoothdata(psth_spx,2,'sgolay', g.smoothvalue);
    % 3. Filter by brain region
    
    psth_spx_PN = psth_spx(idx_PN,:);

    % save full (time x) PSTH for latency detection later
    psthZ_full{hmp} = psth_spx_PN;

    % data used for k-means (only ±0.5 s around CS/US onset, concatenated across stims)
    PSTHall = [PSTHall psth_spx_PN(:,g.roi_pca_idx)]; %#ok<AGROW>
    
    % Store the ROI data for later cluster-wise plotting
    PSTH_ROI_per_stim{hmp} = psth_spx_PN(:,g.roi_pca_idx);
end

n_neurons = size(PSTHall, 1);

% Remove any rows with NaN for clustering (this is the data matrix used below)
nan_rows = any(isnan(PSTHall), 2);
PSTHall_clean = PSTHall(~nan_rows, :);

%% K-means clustering (using the chosen g.clustnum)
fprintf('Running final k-means clustering with k=%d...\n', g.clustnum);

% Run k-means with multiple replicates to find best solution
rng(1); % for reproducibility
[Clusters_clean, Centroids, sumd, D] = kmeans(PSTHall_clean, g.clustnum, ...
    'Distance', 'sqeuclidean', ...
    'Replicates', 20, ...
    'MaxIter', 500, ...
    'Display', 'final'); %#ok<ASGLU>

% Map back to original indexing
Clusters = nan(n_neurons, 1);
Clusters(~nan_rows) = Clusters_clean;

%% === Compute response onset latencies and determine sort order ===
event_inds = g.roi; 

lat_per_stim = nan(n_neurons, numel(ttl));
for s = 1:numel(ttl)
    Z = psthZ_full{s}; % [n_neurons x time]
    for n = 1:n_neurons
        if all(isnan(Z(n,event_inds)))
            lat_per_stim(n,s) = NaN;
        else
            lat_per_stim(n,s) = first_onset_latency(Z(n,:), event_inds, g.onset_threshold, g.min_consec_bins, g.bin_time);
        end
    end
end

min_lat = nan(n_neurons,1);
for n = 1:n_neurons
    l = lat_per_stim(n,:);
    l = l(~isnan(l));
    if ~isempty(l)
        min_lat(n) = min(l);
    else
        min_lat(n) = NaN; % no onset detected in any stim
    end
end

%% === Order neurons WITHIN each cluster by onset latency (earliest first) ===
leafOrder = [];
for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    if ~isempty(clust_idx)
        % Only use units that were part of the clustering process for distance calculation
        clust_idx_clean_mask = ~nan_rows(clust_idx);
        clust_idx_clean = clust_idx(clust_idx_clean_mask);
        
        dists_to_cent = ones(numel(clust_idx),1) * inf;
        if ~isempty(clust_idx_clean)
            % distance to centroid only for the non-NaN units
            dists_to_cent(clust_idx_clean_mask) = pdist2(PSTHall(clust_idx_clean, :), Centroids(c, :), 'euclidean');
        end

        lats_c = min_lat(clust_idx);
        % Sort by: hasLatency (non-NaN first), latency (ascending), distance to centroid (ascending)
        sort_matrix = [isnan(lats_c), lats_c, dists_to_cent];
        [~, sort_idx] = sortrows(sort_matrix, [1 2 3]);
        leafOrder = [leafOrder; clust_idx(sort_idx)]; %#ok<AGROW>

        fprintf('Cluster %d: %d units, %d with detected onset (any stimulus).\n', c, numel(clust_idx), sum(~isnan(lats_c)));
    end
end

leafOrderfl = leafOrder;
colors_clusters = lines(g.clustnum);

%% 2. Feature Comparison: Extract and Compare Key Response Features
features = struct();
feature_names = {'Peak_Response', 'Onset_Latency', 'Response_Duration', 'AUC', 'Time_to_Peak'};

for s = 1:numel(ttl)
    Z = psthZ_full{s}; % Full z-scored PSTH

    for n = 1:n_neurons
        % Peak response (max z-score in ROI)
        features(n).Peak_Response(s) = max(Z(n, event_inds));

        % Onset latency (already computed)
        features(n).Onset_Latency(s) = lat_per_stim(n, s);

        % Response duration (time above threshold)
        above_thresh = Z(n, event_inds) >= g.onset_threshold;
        features(n).Response_Duration(s) = sum(above_thresh) * g.bin_time;

        % Area under curve (sum of positive z-scores)
        features(n).AUC(s) = sum(max(Z(n, event_inds), 0)) * g.bin_time;

        % Time to peak (latency of maximum response)
        [~, peak_idx] = max(Z(n, event_inds));
        features(n).Time_to_Peak(s) = (peak_idx - 1) * g.bin_time;

        % Store cluster assignment
        features(n).Cluster = Clusters(n);
    end
end

% Create violin/box plots for each feature


fig_features = figure('Position', [150, 150, 1600, 900]);
t_feat = tiledlayout(fig_features, 3, numel(feature_names), 'Padding', 'compact', 'TileSpacing', 'compact');

for f = 1:numel(feature_names)
    for s = 1:numel(ttl)
        ax = nexttile(t_feat, (s-1)*numel(feature_names) + f);

        % Collect data for each cluster
        data_by_cluster = cell(1, g.clustnum);
        for c = 1:g.clustnum
            clust_idx = find(Clusters == c);
            vals = [features(clust_idx).(feature_names{f})];
            vals = vals(s:numel(ttl):end); % Get values for stimulus s
            data_by_cluster{c} = vals(~isnan(vals)); % Remove NaNs
        end

        % Create box plot
        hold on;
        positions = 1:g.clustnum;
        for c = 1:g.clustnum
            if ~isempty(data_by_cluster{c})
                boxplot_data = data_by_cluster{c}(:);
                bp = boxplot(boxplot_data, 'Positions', c, 'Width', 0.6, 'Colors', colors_clusters(c,:), ...
                    'Symbol', '', 'Whisker', 1.5);
                set(bp, 'LineWidth', 1.5);

                % Add individual points with jitter
                x_jitter = c + 0.15 * (rand(size(boxplot_data)) - 0.5);
                scatter(x_jitter, boxplot_data, 20, colors_clusters(c,:), 'filled', 'MarkerFaceAlpha', 0.3);
            end
        end

        xlim([0.5, g.clustnum + 0.5]);
        xticks(1:g.clustnum);
        xticklabels(arrayfun(@(x) sprintf('C%d', x), 1:g.clustnum, 'UniformOutput', false));

        if s == 1
            title(strrep(feature_names{f}, '_', ' '), 'FontSize', g.fontSize1);
        end

        if f == 1
            ylabel(hmptitles{s}, 'FontSize', g.fontSize2, 'FontWeight', 'bold');
        end

        if s == numel(ttl)
            xlabel('Cluster', 'FontSize', g.fontSize2);
        end

        grid on;
        set(gca, 'FontSize', g.fontSize2-1);
    end
end

t_feat.Title.String = 'Response Feature Comparison Across Clusters';
t_feat.Title.FontSize = g.fontSize1 + 2;

%% 4. Cluster Centroid Comparison
fig_centroid = figure('Position', [250, 250, 1200, 400]);
t_cent = tiledlayout(fig_centroid, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for s = 1:numel(ttl)
    nexttile(t_cent, s);
    hold on;
    
    % Plot centroids for each cluster
    stim_offset = (s-1) * length(g.roi_pca_idx);
    for c = 1:g.clustnum
        centroid_segment = Centroids(c, stim_offset + (1:length(g.roi_pca_idx)));
        plot(g.timeaxis_pca, centroid_segment, 'Color', colors_clusters(c,:), ...
            'LineWidth', 2.5, 'DisplayName', sprintf('C%d', c));
    end
    
    xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha', 0.5);
    yline(0, ':', 'k', 'LineWidth', 1);
    
    xlabel('Time from Onset (s)', 'FontSize', g.fontSize2);
    ylabel('Z-score (Centroid)', 'FontSize', g.fontSize2);
    title(hmptitles{s}, 'FontSize', g.fontSize1);
    legend('show', 'Location', 'best', 'FontSize', g.fontSize2-2);
    grid on;
    xlim([-0.5, 0.5]);
    set(gca, 'FontSize', g.fontSize2);
end

t_cent.Title.String = 'K-means Cluster Centroids (Prototype Responses)';
t_cent.Title.FontSize = g.fontSize1 + 2;

%% Plot heatmaps and mean PSTH per cluster

fig = figure('Position', [400, 100, 1200, 800]); 
t = tiledlayout(fig,g.clustnum,3,'TileSpacing', 'compact', 'Padding', 'none'); 

hmp_axs = gobjects(1, numel(ttl)); 

for hmp = 1:numel(ttl)
    % Get the full PSTH for heatmap plotting window (+/- 0.5 to +1s around onset)
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    psth_spx = zscore(psth_spx,0,2);
    psth_spx = smoothdata(psth_spx,2,'sgolay', g.smoothvalue);
    psth_spx_PN = psth_spx(idx_PN,:);
    psth_spx_sorted = [psth_spx_PN(leafOrderfl,:)];
    matrix = psth_spx_sorted(:,(g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    
    % The heatmap will span all rows of the tiled layout
    ax = nexttile(t,hmp,[g.clustnum 1]);
    hmp_axs(hmp) = ax; 
    
    imagesc(g.timeaxis_hmp,1:size(matrix,1),matrix);
    clim(g.clim)
    colormap(g.colors.Heatmap);
    xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
    
    if hmp == 1
        ylabel('Cell number')
    else
        set(gca, 'YTickLabel', []);
    end
    
    if hmp == 3
        xlabel('Time from Stimulus Onset (s)')
    end
    
    set(gca, 'FontSize', g.fontSize2);
    title(hmptitles{hmp}, 'FontSize', g.fontSize1)
    
    hold on
    % Draw lines separating the clusters on the heatmap
    n_clu = find(diff(Clusters(leafOrderfl))~=0);
    yline(n_clu+0.5, 'Color', 'k', 'LineWidth', 1);
    hold off
end

% Add colorbar to the last heatmap plot
cb = colorbar(hmp_axs(end), 'westoutside', 'FontSize', g.fontSize2);
ylabel(cb, 'Z-score Firing Rate', 'FontSize', g.fontSize2);



t.Title.String = sprintf('K-means Clustering of LA Neurons (k=%d), Sorted by Earliest Onset Latency', g.clustnum);
t.Title.FontSize = g.fontSize1 + 2;
%% === FIGURE 10: CS vs US comparison (Peak and AUC) with cluster sizes ===

% Calculate mean peak responses and AUC for CS and US for each cluster
cs_response = zeros(g.clustnum, 1);
us_response = zeros(g.clustnum, 1);
cs_auc = zeros(g.clustnum, 1);
us_auc = zeros(g.clustnum, 1);
cluster_sizes = zeros(g.clustnum, 1);

for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    cluster_sizes(c) = length(clust_idx);
    
    if ~isempty(clust_idx)
        % CS response (stimulus 1)
        cs_data = PSTH_ROI_per_stim{1}(clust_idx, :);
        cs_response(c) = nanmean(max(cs_data, [], 2));
        % CS AUC (area under curve - sum of positive z-scores)
        cs_auc(c) = nanmean(sum(max(cs_data, 0), 2) * g.bin_time);
        
        % US response (stimulus 2)
        us_data = PSTH_ROI_per_stim{2}(clust_idx, :);
        us_response(c) = nanmean(max(us_data, [], 2));
        % US AUC
        us_auc(c) = nanmean(sum(max(us_data, 0), 2) * g.bin_time);
    end
end

fig_comparison = figure('Position', [100, 100, 1400, 600]);
t_comp = tiledlayout(fig_comparison, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Normalize cluster sizes for marker sizes (50 to 500)
min_marker_size = 100;
max_marker_size = 1500;
if max(cluster_sizes) > min(cluster_sizes)
    normalized_sizes = min_marker_size + (cluster_sizes - min(cluster_sizes)) * ...
        (max_marker_size - min_marker_size) / (max(cluster_sizes) - min(cluster_sizes));
else
    normalized_sizes = ones(g.clustnum, 1) * mean([min_marker_size, max_marker_size]);
end

% ---- Left Panel: Peak Response ----
nexttile(t_comp, 1);
hold on;

% Plot each cluster with size proportional to number of neurons
for c = 1:g.clustnum
    scatter(cs_response(c), us_response(c), normalized_sizes(c), colors_clusters(c,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
    
    % Add cluster label
    text(cs_response(c), us_response(c), sprintf('  C%d\n  (n=%d)', c, cluster_sizes(c)), ...
        'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
end

% Add unity line (equal response to CS and US)
max_val = max([cs_response; us_response]);
min_val = min([cs_response; us_response]);
plot([min_val max_val], [min_val max_val], '--k', 'LineWidth', 1.5);

% Add quadrant lines at zero
xline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

% Labels and formatting
xlabel('CS Peak Response (Z-score)', 'FontSize', g.fontSize1);
ylabel('US Peak Response (Z-score)', 'FontSize', g.fontSize1);
title('Peak Response', 'FontSize', g.fontSize1);

% Add region labels
text(0.02, 0.98, 'CS > US', 'Units', 'normalized', 'FontSize', 10, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', [0.4 0.4 0.4]);
text(0.98, 0.02, 'US > CS', 'Units', 'normalized', 'FontSize', 10, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', [0.4 0.4 0.4]);

axis equal;
grid on;
set(gca, 'FontSize', g.fontSize2);

% Add legend for cluster sizes
legend_text = sprintf('Marker size indicates\ncluster size (n=%d-%d)', min(cluster_sizes), max(cluster_sizes));
text(0.98, 0.98, legend_text, 'Units', 'normalized', 'FontSize', 9, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

% ---- Right Panel: AUC (Area Under Curve) ----
nexttile(t_comp, 2);
hold on;

% Plot each cluster with size proportional to number of neurons
for c = 1:g.clustnum
    scatter(cs_auc(c), us_auc(c), normalized_sizes(c), colors_clusters(c,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
    
    % Add cluster label
    text(cs_auc(c), us_auc(c), sprintf('  C%d\n  (n=%d)', c, cluster_sizes(c)), ...
        'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
end

% Add unity line
max_val_auc = max([cs_auc; us_auc]);
min_val_auc = min([cs_auc; us_auc]);
plot([min_val_auc max_val_auc], [min_val_auc max_val_auc], '--k', 'LineWidth', 1.5);

% Add quadrant lines at zero
xline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

% Labels and formatting
xlabel('CS AUC (Z-score × s)', 'FontSize', g.fontSize1);
ylabel('US AUC (Z-score × s)', 'FontSize', g.fontSize1);
title('Area Under Curve', 'FontSize', g.fontSize1);

% Add region labels
text(0.02, 0.98, 'CS > US', 'Units', 'normalized', 'FontSize', 10, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', [0.4 0.4 0.4]);
text(0.98, 0.02, 'US > CS', 'Units', 'normalized', 'FontSize', 10, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', [0.4 0.4 0.4]);

axis equal;
grid on;
set(gca, 'FontSize', g.fontSize2);

% Add legend for cluster sizes
text(0.98, 0.98, legend_text, 'Units', 'normalized', 'FontSize', 9, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);

% Overall title
t_comp.Title.String = sprintf('%s: Cluster Responses to CS vs US (k=%d) - Marker size = cluster size', g.bR, g.clustnum);
t_comp.Title.FontSize = g.fontSize1 + 1;

hold off;
fprintf('CS vs US comparison plot created for %s\n', g.bR);

% fprintf('Detailed CS vs US comparison plot created for %s\n', g.bR);
%% === Local function: onset latency detector ===
function lat = first_onset_latency(z_trace, event_inds, threshold, min_consec, bin_time)
% z_trace: 1 x T vector (z-scored, smoothed PSTH)
% event_inds: indices for [0, test_time]
% threshold: z-score threshold (e.g., 3)
% min_consec: minimum consecutive bins above threshold to accept as onset
% bin_time: bin width in seconds
    seg = z_trace(event_inds);
    if any(isnan(seg))
        lat = NaN;
        return
    end
    isAbove = seg >= threshold;
    % moving window: require min_consec consecutive bins above threshold
    winsum = conv(double(isAbove), ones(1,min_consec), 'same');
    idx = find(winsum >= min_consec, 1, 'first');
    if isempty(idx)
        lat = NaN;
    else
        lat = (idx-1) * bin_time; % seconds from t=0
    end
end



%% Stimulus Integration Analysis
% This analysis examines how neurons integrate CS and US information
% by comparing the observed CS+US response to the predicted linear sum
% of CS and US responses separately.
%
% Key metrics:
% 1. Integration Index: (CS+US) - (CS + US) 
%    - Positive = Supra-additive (synergistic)
%    - Near zero = Linear/additive
%    - Negative = Sub-additive (suppressive)
% 2. Peak-based and AUC-based integration measures
% 3. Temporal evolution of integration effects

fprintf('\n=== STIMULUS INTEGRATION ANALYSIS ===\n\n');

%% Calculate Integration Metrics for Each Neuron

integration_metrics = struct();

for n = 1:n_neurons
    % Get responses for each stimulus condition (ROI only)
    cs_trace = PSTH_ROI_per_stim{1}(n, :);      % CS alone
    us_trace = PSTH_ROI_per_stim{2}(n, :);      % US alone
    csus_trace = PSTH_ROI_per_stim{3}(n, :);    % CS+US combined
    
    % Linear prediction: CS + US
    predicted_linear = cs_trace + us_trace;
    
    % --- 1. Peak-based Integration ---
    cs_peak = max(cs_trace);
    us_peak = max(us_trace);
    csus_peak = max(csus_trace);
    predicted_peak = cs_peak + us_peak;
    
    % Integration index (peak-based)
    integration_metrics(n).peak_integration = csus_peak - predicted_peak;
    integration_metrics(n).peak_integration_ratio = csus_peak / max(predicted_peak, eps); % avoid division by zero
    
    % Store individual peaks
    integration_metrics(n).cs_peak = cs_peak;
    integration_metrics(n).us_peak = us_peak;
    integration_metrics(n).csus_peak = csus_peak;
    integration_metrics(n).predicted_peak = predicted_peak;
    
    % --- 2. AUC-based Integration ---
    cs_auc = sum(max(cs_trace, 0)) * g.bin_time;
    us_auc = sum(max(us_trace, 0)) * g.bin_time;
    csus_auc = sum(max(csus_trace, 0)) * g.bin_time;
    predicted_auc = cs_auc + us_auc;
    
    integration_metrics(n).auc_integration = csus_auc - predicted_auc;
    integration_metrics(n).auc_integration_ratio = csus_auc / max(predicted_auc, eps);
    
    % Store individual AUCs
    integration_metrics(n).cs_auc = cs_auc;
    integration_metrics(n).us_auc = us_auc;
    integration_metrics(n).csus_auc = csus_auc;
    integration_metrics(n).predicted_auc = predicted_auc;
    
    % --- 3. Time-resolved Integration ---
    % Calculate integration at each time point
    integration_metrics(n).temporal_integration = csus_trace - predicted_linear;
    integration_metrics(n).temporal_observed = csus_trace;
    integration_metrics(n).temporal_predicted = predicted_linear;
    
    % --- 4. Classify Integration Type ---
    % Using peak-based metric with threshold
    threshold = 1.0; % z-score units (can be adjusted)
    if integration_metrics(n).peak_integration > threshold
        integration_metrics(n).type = 'Supra-additive';
        integration_metrics(n).type_code = 1;
    elseif integration_metrics(n).peak_integration < -threshold
        integration_metrics(n).type = 'Sub-additive';
        integration_metrics(n).type_code = -1;
    else
        integration_metrics(n).type = 'Linear';
        integration_metrics(n).type_code = 0;
    end
    
    % Store cluster assignment
    integration_metrics(n).cluster = Clusters(n);
end

%% Summarize Integration by Cluster

fprintf('Integration Summary by Cluster:\n');
fprintf('%-10s %-15s %-15s %-15s %-15s\n', 'Cluster', 'N Neurons', 'Mean Peak Int', 'Mean AUC Int', 'Integration Type');
fprintf('%s\n', repmat('-', 1, 80));

cluster_integration = struct();

for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    
    if ~isempty(clust_idx)
        % Extract metrics
        peak_int = [integration_metrics(clust_idx).peak_integration];
        auc_int = [integration_metrics(clust_idx).auc_integration];
        types = [integration_metrics(clust_idx).type_code];
        
        % Calculate statistics
        cluster_integration(c).mean_peak_int = nanmean(peak_int);
        cluster_integration(c).sem_peak_int = nanstd(peak_int) / sqrt(sum(~isnan(peak_int)));
        cluster_integration(c).mean_auc_int = nanmean(auc_int);
        cluster_integration(c).sem_auc_int = nanstd(auc_int) / sqrt(sum(~isnan(auc_int)));
        
        % Count integration types
        n_supra = sum(types == 1);
        n_linear = sum(types == 0);
        n_sub = sum(types == -1);
        
        cluster_integration(c).n_supra = n_supra;
        cluster_integration(c).n_linear = n_linear;
        cluster_integration(c).n_sub = n_sub;
        
        % Dominant type
        [~, max_idx] = max([n_supra, n_linear, n_sub]);
        type_labels = {'Supra-add', 'Linear', 'Sub-add'};
        dominant_type = type_labels{max_idx};
        
        fprintf('C%-9d %-15d %-15.2f %-15.2f %-15s\n', ...
            c, length(clust_idx), cluster_integration(c).mean_peak_int, ...
            cluster_integration(c).mean_auc_int, dominant_type);
    end
end

%% FIGURE 1: Integration Overview - Observed vs Predicted Responses

fig_overview = figure('Position', [50, 50, 1400, 600]);
t_overview = tiledlayout(fig_overview, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

colors_clusters = lines(g.clustnum);

% --- Panel A: Peak Response ---
nexttile(t_overview, 1);
hold on;

% Unity line
max_val = max([[integration_metrics.csus_peak], [integration_metrics.predicted_peak]]);
min_val = min([[integration_metrics.csus_peak], [integration_metrics.predicted_peak]]);
plot([min_val max_val], [min_val max_val], '--k', 'LineWidth', 2, 'DisplayName', 'Linear (Unity)');

% Plot each cluster
for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    observed = [integration_metrics(clust_idx).csus_peak];
    predicted = [integration_metrics(clust_idx).predicted_peak];
    
    scatter(predicted, observed, 80, colors_clusters(c,:), 'filled', ...
        'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('C%d', c));
end

% Add region shading
xl = xlim; yl = ylim;
patch([xl(1) xl(2) xl(2) xl(1)], [xl(1) xl(2)+10 yl(2) yl(1)], [0.9 0.9 1], ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.98, 0.02, 'Supra-additive', 'Units', 'normalized', 'FontSize', 10, ...
    'Color', [0 0 0.5], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

patch([xl(1) xl(2) xl(2) xl(1)], [yl(1) xl(1) xl(2) xl(2)], [1 0.9 0.9], ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.98, 0.98, 'Sub-additive', 'Units', 'normalized', 'FontSize', 10, ...
    'Color', [0.5 0 0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

xlabel('Predicted Peak (CS + US, Z-score)', 'FontSize', g.fontSize1);
ylabel('Observed Peak (CS+US, Z-score)', 'FontSize', g.fontSize1);
title('Peak Response Integration', 'FontSize', g.fontSize1);
legend('show', 'Location', 'northwest', 'FontSize', g.fontSize2-2);
axis equal;
grid on;
set(gca, 'FontSize', g.fontSize2);

% --- Panel B: AUC ---
nexttile(t_overview, 2);
hold on;

% Unity line
max_val_auc = max([[integration_metrics.csus_auc], [integration_metrics.predicted_auc]]);
min_val_auc = min([[integration_metrics.csus_auc], [integration_metrics.predicted_auc]]);
plot([min_val_auc max_val_auc], [min_val_auc max_val_auc], '--k', 'LineWidth', 2, 'DisplayName', 'Linear (Unity)');

% Plot each cluster
for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    observed_auc = [integration_metrics(clust_idx).csus_auc];
    predicted_auc = [integration_metrics(clust_idx).predicted_auc];
    
    scatter(predicted_auc, observed_auc, 80, colors_clusters(c,:), 'filled', ...
        'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('C%d', c));
end

% Add region shading
xl = xlim; yl = ylim;
patch([xl(1) xl(2) xl(2) xl(1)], [xl(1) xl(2)+10 yl(2) yl(1)], [0.9 0.9 1], ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.98, 0.02, 'Supra-additive', 'Units', 'normalized', 'FontSize', 10, ...
    'Color', [0 0 0.5], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

patch([xl(1) xl(2) xl(2) xl(1)], [yl(1) xl(1) xl(2) xl(2)], [1 0.9 0.9], ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(0.98, 0.98, 'Sub-additive', 'Units', 'normalized', 'FontSize', 10, ...
    'Color', [0.5 0 0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

xlabel('Predicted AUC (CS + US, Z×s)', 'FontSize', g.fontSize1);
ylabel('Observed AUC (CS+US, Z×s)', 'FontSize', g.fontSize1);
title('Area Under Curve Integration', 'FontSize', g.fontSize1);
legend('show', 'Location', 'northwest', 'FontSize', g.fontSize2-2);
axis equal;
grid on;
set(gca, 'FontSize', g.fontSize2);

t_overview.Title.String = sprintf('%s: Stimulus Integration - Observed vs Predicted Responses', g.bR);
t_overview.Title.FontSize = g.fontSize1 + 2;

%% FIGURE 2: Integration Index Distribution by Cluster

fig_distribution = figure('Position', [100, 100, 1400, 800]);
t_dist = tiledlayout(fig_distribution, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- Panel A: Peak Integration Index Distribution ---
nexttile(t_dist, 1);
hold on;

for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    peak_int = [integration_metrics(clust_idx).peak_integration];
    peak_int = peak_int(~isnan(peak_int));
    
    if ~isempty(peak_int)
        % Violin plot / histogram
        [counts, edges] = histcounts(peak_int, 20);
        centers = edges(1:end-1) + diff(edges)/2;
        
        % Normalize for visualization
        counts_norm = counts / max(counts) * 0.4; % scale to 0.4 cluster width
        
        % Plot as filled area
        x_fill = [c - counts_norm, fliplr(c + counts_norm)];
        y_fill = [centers, fliplr(centers)];
        fill(x_fill, y_fill, colors_clusters(c,:), 'FaceAlpha', 0.5, 'EdgeColor', colors_clusters(c,:), 'LineWidth', 1.5);
        
        % Add mean line
        mean_int = nanmean(peak_int);
        plot([c-0.4, c+0.4], [mean_int, mean_int], 'k', 'LineWidth', 2);
    end
end

yline(0, '--k', 'LineWidth', 1.5, 'Label', 'Linear', 'LabelHorizontalAlignment', 'left');
yline(1, ':', 'Color', [0 0 0.8], 'LineWidth', 1, 'Label', 'Supra threshold');
yline(-1, ':', 'Color', [0.8 0 0], 'LineWidth', 1, 'Label', 'Sub threshold');

xlim([0.5, g.clustnum + 0.5]);
xticks(1:g.clustnum);
xticklabels(arrayfun(@(x) sprintf('C%d', x), 1:g.clustnum, 'UniformOutput', false));
xlabel('Cluster', 'FontSize', g.fontSize2);
ylabel('Peak Integration Index (Z-score)', 'FontSize', g.fontSize2);
title('Peak-based Integration', 'FontSize', g.fontSize1);
grid on;
set(gca, 'FontSize', g.fontSize2);

% --- Panel B: AUC Integration Index Distribution ---
nexttile(t_dist, 2);
hold on;

for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    auc_int = [integration_metrics(clust_idx).auc_integration];
    auc_int = auc_int(~isnan(auc_int));
    
    if ~isempty(auc_int)
        [counts, edges] = histcounts(auc_int, 20);
        centers = edges(1:end-1) + diff(edges)/2;
        counts_norm = counts / max(counts) * 0.4;
        
        x_fill = [c - counts_norm, fliplr(c + counts_norm)];
        y_fill = [centers, fliplr(centers)];
        fill(x_fill, y_fill, colors_clusters(c,:), 'FaceAlpha', 0.5, 'EdgeColor', colors_clusters(c,:), 'LineWidth', 1.5);
        
        mean_int = nanmean(auc_int);
        plot([c-0.4, c+0.4], [mean_int, mean_int], 'k', 'LineWidth', 2);
    end
end

yline(0, '--k', 'LineWidth', 1.5, 'Label', 'Linear', 'LabelHorizontalAlignment', 'left');

xlim([0.5, g.clustnum + 0.5]);
xticks(1:g.clustnum);
xticklabels(arrayfun(@(x) sprintf('C%d', x), 1:g.clustnum, 'UniformOutput', false));
xlabel('Cluster', 'FontSize', g.fontSize2);
ylabel('AUC Integration Index (Z×s)', 'FontSize', g.fontSize2);
title('AUC-based Integration', 'FontSize', g.fontSize1);
grid on;
set(gca, 'FontSize', g.fontSize2);

% --- Panel C: Integration Type Proportions ---
nexttile(t_dist, 3);
hold on;

type_matrix = zeros(g.clustnum, 3); % columns: supra, linear, sub
for c = 1:g.clustnum
    type_matrix(c, 1) = cluster_integration(c).n_supra;
    type_matrix(c, 2) = cluster_integration(c).n_linear;
    type_matrix(c, 3) = cluster_integration(c).n_sub;
end

% Normalize to proportions
type_proportions = type_matrix ./ sum(type_matrix, 2);

% Stacked bar chart
b = bar(1:g.clustnum, type_proportions, 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.7);
b(1).CData = repmat([0.2 0.4 0.9], g.clustnum, 1); % Supra = blue
b(2).CData = repmat([0.7 0.7 0.7], g.clustnum, 1); % Linear = gray
b(3).CData = repmat([0.9 0.4 0.2], g.clustnum, 1); % Sub = red

xlabel('Cluster', 'FontSize', g.fontSize2);
ylabel('Proportion of Neurons', 'FontSize', g.fontSize2);
title('Integration Type Distribution', 'FontSize', g.fontSize1);
legend({'Supra-additive', 'Linear', 'Sub-additive'}, 'Location', 'best', 'FontSize', g.fontSize2-2);
xticks(1:g.clustnum);
xticklabels(arrayfun(@(x) sprintf('C%d', x), 1:g.clustnum, 'UniformOutput', false));
ylim([0 1]);
grid on;
set(gca, 'FontSize', g.fontSize2);

% --- Panel D: Mean Integration by Cluster (with error bars) ---
nexttile(t_dist, 4);
hold on;

mean_peak_ints = zeros(g.clustnum, 1);
sem_peak_ints = zeros(g.clustnum, 1);

for c = 1:g.clustnum
    mean_peak_ints(c) = cluster_integration(c).mean_peak_int;
    sem_peak_ints(c) = cluster_integration(c).sem_peak_int;
end

% Bar plot with error bars
b = bar(1:g.clustnum, mean_peak_ints, 'FaceColor', 'flat', 'BarWidth', 0.6);
for c = 1:g.clustnum
    b.CData(c,:) = colors_clusters(c,:);
end

errorbar(1:g.clustnum, mean_peak_ints, sem_peak_ints, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

yline(0, '--k', 'LineWidth', 1.5);

xlabel('Cluster', 'FontSize', g.fontSize2);
ylabel('Mean Peak Integration Index (Z)', 'FontSize', g.fontSize2);
title('Mean Integration by Cluster', 'FontSize', g.fontSize1);
xticks(1:g.clustnum);
xticklabels(arrayfun(@(x) sprintf('C%d', x), 1:g.clustnum, 'UniformOutput', false));
grid on;
set(gca, 'FontSize', g.fontSize2);

t_dist.Title.String = sprintf('%s: Integration Index Distribution Across Clusters', g.bR);
t_dist.Title.FontSize = g.fontSize1 + 2;

%% FIGURE 3: Temporal Evolution of Integration

fig_temporal = figure('Position', [150, 150, 1600, 800]);
t_temp = tiledlayout(fig_temporal, g.clustnum, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

for c = 1:g.clustnum
    ax = nexttile(t_temp, c);
    hold on;
    
    clust_idx = find(Clusters == c);
    
    if ~isempty(clust_idx)
        % Calculate mean temporal traces
        observed_traces = zeros(length(clust_idx), length(g.timeaxis_pca));
        predicted_traces = zeros(length(clust_idx), length(g.timeaxis_pca));
        integration_traces = zeros(length(clust_idx), length(g.timeaxis_pca));
        
        for i = 1:length(clust_idx)
            n = clust_idx(i);
            observed_traces(i,:) = integration_metrics(n).temporal_observed;
            predicted_traces(i,:) = integration_metrics(n).temporal_predicted;
            integration_traces(i,:) = integration_metrics(n).temporal_integration;
        end
        
        % Calculate mean and SEM
        mean_observed = nanmean(observed_traces, 1);
        sem_observed = nanstd(observed_traces, 0, 1) / sqrt(size(observed_traces, 1));
        
        mean_predicted = nanmean(predicted_traces, 1);
        sem_predicted = nanstd(predicted_traces, 0, 1) / sqrt(size(predicted_traces, 1));
        
        mean_integration = nanmean(integration_traces, 1);
        sem_integration = nanstd(integration_traces, 0, 1) / sqrt(size(integration_traces, 1));
        
        % Plot predicted (CS + US)
        fill([g.timeaxis_pca, fliplr(g.timeaxis_pca)], ...
             [mean_predicted + sem_predicted, fliplr(mean_predicted - sem_predicted)], ...
             [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(g.timeaxis_pca, mean_predicted, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 2, 'DisplayName', 'Predicted (CS+US)');
        
        % Plot observed (CS+US)
        fill([g.timeaxis_pca, fliplr(g.timeaxis_pca)], ...
             [mean_observed + sem_observed, fliplr(mean_observed - sem_observed)], ...
             colors_clusters(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(g.timeaxis_pca, mean_observed, '-', 'Color', colors_clusters(c,:), 'LineWidth', 2, 'DisplayName', 'Observed (CS+US)');
        
        % Add integration difference in separate subplot or as shaded area
        % Plot integration as area below zero line
        integration_color = mean_integration;
        pos_int = mean_integration; pos_int(pos_int < 0) = 0;
        neg_int = mean_integration; neg_int(neg_int > 0) = 0;
        
        area(g.timeaxis_pca, pos_int, 'FaceColor', [0.2 0.4 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Supra-add');
        area(g.timeaxis_pca, neg_int, 'FaceColor', [0.9 0.4 0.2], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Sub-add');
        
        % Reference lines
        xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.7);
        yline(0, ':', 'k', 'LineWidth', 1);
        
        ylabel(sprintf('C%d (n=%d)\nZ-score', c, length(clust_idx)), 'FontSize', g.fontSize2);
        xlim([-0.5, 0.5]);
        
        if c == 1
            legend('show', 'Location', 'best', 'FontSize', g.fontSize2-2);
            title('Temporal Evolution: Observed vs Predicted', 'FontSize', g.fontSize1);
        end
        
        if c == g.clustnum
            xlabel('Time from Stimulus Onset (s)', 'FontSize', g.fontSize2);
        else
            set(gca, 'XTickLabel', []);
        end
        
        grid on;
        set(gca, 'FontSize', g.fontSize2);
    end
end

t_temp.Title.String = sprintf('%s: Time-Resolved Integration Analysis by Cluster', g.bR);
t_temp.Title.FontSize = g.fontSize1 + 2;

%% FIGURE 4: Integration Scatter Matrix

fig_scatter = figure('Position', [200, 200, 1000, 900]);
t_scatter = tiledlayout(fig_scatter, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- Panel A: CS vs US Response (colored by integration) ---
nexttile(t_scatter, 1);
hold on;

cs_peaks = [integration_metrics.cs_peak];
us_peaks = [integration_metrics.us_peak];
peak_ints = [integration_metrics.peak_integration];

% Color by integration index
scatter(cs_peaks, us_peaks, 80, peak_ints, 'filled', 'MarkerFaceAlpha', 0.7);
colormap(gca, bluewhitered(256));
cb = colorbar;
ylabel(cb, 'Integration Index', 'FontSize', g.fontSize2-1);
caxis([-max(abs(peak_ints)), max(abs(peak_ints))]);

xlabel('CS Peak Response (Z)', 'FontSize', g.fontSize2);
ylabel('US Peak Response (Z)', 'FontSize', g.fontSize2);
title('CS vs US Colored by Integration', 'FontSize', g.fontSize1);
grid on;
set(gca, 'FontSize', g.fontSize2);

% --- Panel B: Response Magnitude vs Integration ---
nexttile(t_scatter, 2);
hold on;

response_magnitude = sqrt(cs_peaks.^2 + us_peaks.^2); % Euclidean norm

for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    scatter(response_magnitude(clust_idx), peak_ints(clust_idx), 80, colors_clusters(c,:), ...
        'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('C%d', c));
end

yline(0, '--k', 'LineWidth', 1.5);

xlabel('Response Magnitude (√(CS² + US²))', 'FontSize', g.fontSize2);
ylabel('Peak Integration Index (Z)', 'FontSize', g.fontSize2);
title('Magnitude vs Integration', 'FontSize', g.fontSize1);
legend('show', 'Location', 'best', 'FontSize', g.fontSize2-3);
grid on;
set(gca, 'FontSize', g.fontSize2);

% --- Panel C: CS Selectivity vs Integration ---
nexttile(t_scatter, 3);
hold on;

% Calculate selectivity: (CS - US) / (CS + US)
cs_selectivity = (cs_peaks - us_peaks) ./ (cs_peaks + us_peaks + eps);

for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    scatter(cs_selectivity(clust_idx), peak_ints(clust_idx), 80, colors_clusters(c,:), ...
        'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('C%d', c));
end

yline(0, '--k', 'LineWidth', 1.5);
xline(0, '--k', 'LineWidth', 1);

xlabel('CS Selectivity [(CS-US)/(CS+US)]', 'FontSize', g.fontSize2);
ylabel('Peak Integration Index (Z)', 'FontSize', g.fontSize2);
title('Selectivity vs Integration', 'FontSize', g.fontSize1);
text(0.98, 0.02, 'US-selective', 'Units', 'normalized', 'FontSize', 9, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
text(0.02, 0.02, 'CS-selective', 'Units', 'normalized', 'FontSize', 9, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
grid on;
set(gca, 'FontSize', g.fontSize2);

% --- Panel D: Integration Type by Cluster (pie charts) ---
nexttile(t_scatter, 4);
axis off;

% Create mini pie charts for each cluster
n_rows = ceil(sqrt(g.clustnum));
n_cols = ceil(g.clustnum / n_rows);

for c = 1:g.clustnum
    % Calculate subplot position
    row = ceil(c / n_cols);
    col = mod(c-1, n_cols) + 1;
    
    % Position within the tile
    width = 1 / n_cols;
    height = 1 / n_rows;
    left = (col - 1) * width;
    bottom = 1 - row * height;
    
    ax_pie = axes('Parent', fig_scatter, 'Position', [0.565 + left*0.4, 0.125 + bottom*0.4, width*0.35, height*0.35]);
    
    pie_data = [cluster_integration(c).n_supra, cluster_integration(c).n_linear, cluster_integration(c).n_sub];
    pie(ax_pie, pie_data);
    colormap(ax_pie, [0.2 0.4 0.9; 0.7 0.7 0.7; 0.9 0.4 0.2]);
    title(ax_pie, sprintf('C%d', c), 'FontSize', g.fontSize2-2);
end

% Add legend for pie charts
axes('Parent', fig_scatter, 'Position', [0.565, 0.05, 0.4, 0.05]);
axis off;
text(0.1, 0.5, '■ Supra-additive', 'FontSize', g.fontSize2-2, 'Color', [0.2 0.4 0.9]);
text(0.4, 0.5, '■ Linear', 'FontSize', g.fontSize2-2, 'Color', [0.7 0.7 0.7]);
text(0.7, 0.5, '■ Sub-additive', 'FontSize', g.fontSize2-2, 'Color', [0.9 0.4 0.2]);

t_scatter.Title.String = sprintf('%s: Integration Analysis - Multi-dimensional Relationships', g.bR);
t_scatter.Title.FontSize = g.fontSize1 + 2;

%% Statistical Testing

fprintf('\n=== STATISTICAL TESTS ===\n\n');

% Test 1: Is integration significantly different from zero for each cluster?
fprintf('One-sample t-test (Integration vs Zero) by Cluster:\n');
fprintf('%-10s %-10s %-10s %-15s %-15s\n', 'Cluster', 'N', 't-stat', 'p-value', 'Sig?');
fprintf('%s\n', repmat('-', 1, 70));

for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    peak_int = [integration_metrics(clust_idx).peak_integration];
    peak_int = peak_int(~isnan(peak_int));
    
    if length(peak_int) > 2
        [h, p, ~, stats] = ttest(peak_int, 0);
        sig_str = '';
        if p < 0.001
            sig_str = '***';
        elseif p < 0.01
            sig_str = '**';
        elseif p < 0.05
            sig_str = '*';
        else
            sig_str = 'n.s.';
        end
        
        fprintf('C%-9d %-10d %-10.3f %-15.4f %-15s\n', c, length(peak_int), stats.tstat, p, sig_str);
    else
        fprintf('C%-9d %-10d %-10s %-15s %-15s\n', c, length(peak_int), 'N/A', 'N/A', 'Too few');
    end
end

% Test 2: Compare integration between clusters (ANOVA)
fprintf('\n\nOne-way ANOVA: Comparing Integration Across Clusters\n');

all_integrations = [];
all_cluster_labels = [];

for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    peak_int = [integration_metrics(clust_idx).peak_integration];
    peak_int = peak_int(~isnan(peak_int));
    
    all_integrations = [all_integrations, peak_int]; %#ok<AGROW>
    all_cluster_labels = [all_cluster_labels, repmat(c, 1, length(peak_int))]; %#ok<AGROW>
end

if length(unique(all_cluster_labels)) > 1
    [p_anova, tbl, stats_anova] = anova1(all_integrations, all_cluster_labels, 'off');
    fprintf('F(%d,%d) = %.3f, p = %.4f\n', tbl{2,3}, tbl{3,3}, tbl{2,5}, p_anova);
    
    if p_anova < 0.05
        fprintf('Significant difference detected. Running post-hoc tests...\n');
        % Perform multiple comparisons
        [c_posthoc, ~, ~, gnames] = multcompare(stats_anova, 'Display', 'off');
        
        fprintf('\nPost-hoc comparisons (Tukey-Kramer):\n');
        fprintf('%-15s %-15s %-15s\n', 'Comparison', 'Mean Diff', 'p-value');
        fprintf('%s\n', repmat('-', 1, 50));
        
        for i = 1:size(c_posthoc, 1)
            if c_posthoc(i, 6) < 0.05
                fprintf('C%d vs C%d      %-15.3f %-15.4f *\n', ...
                    c_posthoc(i,1), c_posthoc(i,2), c_posthoc(i,4), c_posthoc(i,6));
            end
        end
    else
        fprintf('No significant differences between clusters.\n');
    end
else
    fprintf('Not enough clusters for ANOVA.\n');
end

% Test 3: Correlation between CS/US response and integration
fprintf('\n\nCorrelation Analysis:\n');

% CS response vs integration
cs_peaks_all = [integration_metrics.cs_peak];
peak_ints_all = [integration_metrics.peak_integration];
valid_idx = ~isnan(cs_peaks_all) & ~isnan(peak_ints_all);

if sum(valid_idx) > 10
    [r_cs, p_cs] = corrcoef(cs_peaks_all(valid_idx), peak_ints_all(valid_idx));
    fprintf('CS Peak vs Integration: r = %.3f, p = %.4f\n', r_cs(1,2), p_cs(1,2));
end

% US response vs integration
us_peaks_all = [integration_metrics.us_peak];
valid_idx = ~isnan(us_peaks_all) & ~isnan(peak_ints_all);

if sum(valid_idx) > 10
    [r_us, p_us] = corrcoef(us_peaks_all(valid_idx), peak_ints_all(valid_idx));
    fprintf('US Peak vs Integration: r = %.3f, p = %.4f\n', r_us(1,2), p_us(1,2));
end

%% Summary Table

fprintf('\n\n=== INTEGRATION SUMMARY TABLE ===\n\n');
fprintf('%-8s | %-8s | %-12s | %-12s | %-10s | %-10s | %-10s\n', ...
    'Cluster', 'N Units', 'Mean Peak', 'Mean AUC', 'N Supra', 'N Linear', 'N Sub');
fprintf('%s\n', repmat('-', 1, 90));

for c = 1:g.clustnum
    clust_idx = find(Clusters == c);
    fprintf('C%-7d | %-8d | %+11.2f | %+11.2f | %-10d | %-10d | %-10d\n', ...
        c, length(clust_idx), ...
        cluster_integration(c).mean_peak_int, ...
        cluster_integration(c).mean_auc_int, ...
        cluster_integration(c).n_supra, ...
        cluster_integration(c).n_linear, ...
        cluster_integration(c).n_sub);
end

fprintf('\n');

%% Export Results (Optional)

% Save integration metrics to file
integration_table = struct2table(integration_metrics);
output_filename = sprintf('%s_integration_analysis.csv', g.bR);
% writetable(integration_table, fullfile(g.mainFolder, output_filename));
% fprintf('Integration metrics saved to: %s\n', output_filename);

fprintf('\n=== STIMULUS INTEGRATION ANALYSIS COMPLETE ===\n\n');

%% Helper function for blue-white-red colormap
function cmap = bluewhitered(n)
    % Creates a blue-white-red colormap for diverging data
    if nargin < 1
        n = 256;
    end
    
    % Create color gradient
    r = [0:1/(n/2-1):1, ones(1,n/2)]';
    g = [0:1/(n/2-1):1, 1:-1/(n/2-1):0]';
    b = [ones(1,n/2), 1:-1/(n/2-1):0]';
    
    cmap = [r, g, b];
end