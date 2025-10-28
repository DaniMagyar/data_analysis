% comparison_heatmap_004
% Clean heatmap-only version: K-means clustering comparing CS and US responses
% PNs and INs clustered separately, plotted together (PNs on top, INs on bottom)

clear all; close all

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

ttl = {'triptest_sound_only','triptest_shocks_only'}; % CS and US only

hmptitles = {'CS', 'US'};

g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
g.cell_metrics = BAfc_match_celltypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 13;
g.fontSize2 = 12;
g.xlinewidth = 2;
g.pre_time = 5;
g.post_time = 1;
g.test_time = 0.5;

% Clustering Parameters
g.clustnum = 5;

% Manual clustering thresholds
g.excitation_threshold = 2;   % Z-score threshold for excitatory response
g.inhibition_threshold = -2;  % Z-score threshold for inhibitory response

% Colormap limits: use percentile to avoid extreme values
g.use_percentile = true;  % Set to false to use min/max
g.clim_percentile = 99;   % Use 2.5th to 97.5th percentile (95% range)

% Onset Detection Parameters
g.testvalue = 3;
g.onset_threshold = g.testvalue;
g.bin_time = 0.001;
g.min_consec_bins = max(1, round(0.02 / g.bin_time));

g.smoothvalue = 101;
g.plotwin = [0.5 1];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);

% ROI for k-means
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
g.roi_pca_idx = round((g.pre_time-0.5)/g.bin_time:(g.pre_time+0.5)/g.bin_time);
g.timeaxis_pca = -0.5:g.bin_time:0.5;

g.bR = 'LA';

% Get PN and IN indices
idx_PN = strcmp(g.cell_metrics.brainRegion, g.bR) & strcmp(g.cell_metrics.putativeCellType, 'PN');
idx_IN = strcmp(g.cell_metrics.brainRegion, g.bR) & strcmp(g.cell_metrics.putativeCellType, 'IN');

%% CALCULATE ALL PSTHs ONCE AT THE BEGINNING

fprintf('Calculating PSTHs once for all stimuli...\n');
psthZ_full = cell(1, numel(ttl));

for hmp = 1:numel(ttl)
    fprintf('  Processing %s...\n', hmptitles{hmp});
    psth_spx_og = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Z-score using only baseline period (before stimulus)
    baseline_idx = 1:(g.pre_time / g.bin_time);
    baseline_mean = mean(psth_spx_og(:, baseline_idx), 2);
    baseline_std = std(psth_spx_og(:, baseline_idx), 0, 2);
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;

    psth_spx = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);

    % Store for all neurons
    psthZ_full{hmp} = psth_spx;
end

fprintf('PSTHs calculated.\n\n');

%% K-MEANS CLUSTERING

cell_types = {'PN', 'IN'};
idx_celltypes = {idx_PN, idx_IN};
results = struct();

for ct = 1:2
    fprintf('=== K-means Clustering for %s ===\n', cell_types{ct});

    idx_curr = idx_celltypes{ct};

    % Prepare data using pre-calculated PSTHs
    PSTHall = [];
    for hmp = 1:numel(ttl)
        psth_spx_curr = psthZ_full{hmp}(idx_curr, :);
        PSTHall = [PSTHall psth_spx_curr(:, g.roi_pca_idx)]; %#ok<AGROW>
    end

    n_neurons = size(PSTHall, 1);

    % Remove NaN rows
    nan_rows = any(isnan(PSTHall), 2);
    PSTHall_clean = PSTHall(~nan_rows, :);

    % K-means clustering
    fprintf('Running k-means clustering with k=%d...\n', g.clustnum);
    rng(1);
    [Clusters_clean, Centroids, ~, ~] = kmeans(PSTHall_clean, g.clustnum, ...
        'Distance', 'sqeuclidean', ...
        'Replicates', 20, ...
        'MaxIter', 500, ...
        'Display', 'final');

    % Map back to original indexing
    Clusters = nan(n_neurons, 1);
    Clusters(~nan_rows) = Clusters_clean;

    % Compute onset and offset latencies using pre-calculated PSTHs
    event_inds = g.roi;
    onset_lat_per_stim = nan(n_neurons, numel(ttl));
    offset_lat_per_stim = nan(n_neurons, numel(ttl));

    for s = 1:numel(ttl)
        Z = psthZ_full{s}(idx_curr, :);
        for n = 1:n_neurons
            if all(isnan(Z(n, event_inds)))
                onset_lat_per_stim(n, s) = NaN;
                offset_lat_per_stim(n, s) = NaN;
            else
                [onset_lat_per_stim(n, s), offset_lat_per_stim(n, s)] = compute_onset_offset_latency(Z(n, :), event_inds, g.onset_threshold, g.min_consec_bins, g.bin_time);
            end
        end
    end

    min_onset_lat = nan(n_neurons, 1);
    min_offset_lat = nan(n_neurons, 1);
    for n = 1:n_neurons
        l_onset = onset_lat_per_stim(n, :);
        l_onset = l_onset(~isnan(l_onset));
        if ~isempty(l_onset)
            min_onset_lat(n) = min(l_onset);
        else
            min_onset_lat(n) = NaN;
        end

        l_offset = offset_lat_per_stim(n, :);
        l_offset = l_offset(~isnan(l_offset));
        if ~isempty(l_offset)
            min_offset_lat(n) = min(l_offset);
        else
            min_offset_lat(n) = NaN;
        end
    end

    % Order neurons within each cluster by onset latency
    leafOrder = [];
    for c = 1:g.clustnum
        clust_idx = find(Clusters == c);
        if ~isempty(clust_idx)
            clust_idx_clean_mask = ~nan_rows(clust_idx);
            clust_idx_clean = clust_idx(clust_idx_clean_mask);

            dists_to_cent = ones(numel(clust_idx), 1) * inf;
            if ~isempty(clust_idx_clean)
                dists_to_cent(clust_idx_clean_mask) = pdist2(PSTHall(clust_idx_clean, :), Centroids(c, :), 'euclidean');
            end

            onset_lats_c = min_onset_lat(clust_idx);
            offset_lats_c = min_offset_lat(clust_idx);
            % Sort by: hasOnset (non-NaN first), onset latency (ascending), offset latency (ascending), distance to centroid
            sort_matrix = [isnan(onset_lats_c), onset_lats_c, offset_lats_c, dists_to_cent];
            [~, sort_idx] = sortrows(sort_matrix, [1 2 3 4]);
            leafOrder = [leafOrder; clust_idx(sort_idx)]; %#ok<AGROW>

            fprintf('Cluster %d: %d units, %d with detected onset.\n', c, numel(clust_idx), sum(~isnan(onset_lats_c)));
        end
    end

    % Store results
    results(ct).leafOrder = leafOrder;
    results(ct).Clusters = Clusters;
    results(ct).idx_curr = idx_curr;
    results(ct).n_neurons = n_neurons;
end

%% Calculate matrices for k-means figure using pre-calculated PSTHs

fprintf('\nPreparing k-means heatmap matrices...\n');
matrices_all = cell(1, numel(ttl));

for hmp = 1:numel(ttl)
    % Use pre-calculated PSTH
    psth_spx = psthZ_full{hmp};

    % Extract and sort PNs
    psth_PN = psth_spx(results(1).idx_curr, :);
    psth_PN_sorted = psth_PN(results(1).leafOrder, :);
    matrix_PN = psth_PN_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);

    % Extract and sort INs
    psth_IN = psth_spx(results(2).idx_curr, :);
    psth_IN_sorted = psth_IN(results(2).leafOrder, :);
    matrix_IN = psth_IN_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);

    % Combine: PNs on top, INs on bottom
    matrix_combined = [matrix_PN; matrix_IN];
    matrices_all{hmp} = matrix_combined;
end

% Determine global color limits
all_values = [matrices_all{:}];
if g.use_percentile
    clim_min = prctile(all_values(:), (100 - g.clim_percentile) / 2);
    clim_max = prctile(all_values(:), 100 - (100 - g.clim_percentile) / 2);
else
    clim_min = min(all_values(:));
    clim_max = max(all_values(:));
end

%% Plot k-means heatmap

fprintf('Plotting k-means clustering figure...\n');
fig = figure('Position', [400, 100, 800, 800]);
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

colors_clusters = lines(g.clustnum);

for hmp = 1:numel(ttl)
    matrix_combined = matrices_all{hmp};

    ax = nexttile(t, hmp);
    imagesc(g.timeaxis_hmp, 1:size(matrix_combined, 1), matrix_combined);
    clim([clim_min clim_max]);
    colormap(g.colors.Heatmap);
    xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha', 1);

    if hmp == 1
        ylabel('Cell number', 'FontSize', g.fontSize2);
    else
        set(gca, 'YTickLabel', []);
    end

    xlabel('Time from Stimulus Onset (s)', 'FontSize', g.fontSize2);
    set(gca, 'FontSize', g.fontSize2);
    title(hmptitles{hmp}, 'FontSize', g.fontSize1);

    hold on;

    % Draw lines separating clusters within PNs
    Clusters_PN = results(1).Clusters(results(1).leafOrder);
    n_clu_PN = find(diff(Clusters_PN) ~= 0);
    yline(n_clu_PN + 0.5, 'Color', colors_clusters(1, :), 'LineWidth', 1);

    % Draw line separating PNs from INs
    n_PN = results(1).n_neurons;
    yline(n_PN + 0.5, 'Color', 'k', 'LineWidth', 3);

    % Draw lines separating clusters within INs
    Clusters_IN = results(2).Clusters(results(2).leafOrder);
    n_clu_IN = find(diff(Clusters_IN) ~= 0);
    yline(n_clu_IN + n_PN + 0.5, 'Color', colors_clusters(2, :), 'LineWidth', 1);

    hold off;
end

% Add colorbar
cb = colorbar(nexttile(t, 2), 'eastoutside', 'FontSize', g.fontSize2);
ylabel(cb, 'Z-score Firing Rate', 'FontSize', g.fontSize2);

title(t, sprintf('%s: K-means Clustering (k=%d) - PNs (top) and INs (bottom)', g.bR, g.clustnum), 'FontSize', g.fontSize1 + 2);

%% MANUAL CLUSTERING using pre-calculated PSTHs

fprintf('\n=== Manual Clustering ===\n');
results_manual = struct();

for ct = 1:2
    fprintf('Manual clustering for %s...\n', cell_types{ct});

    idx_curr = idx_celltypes{ct};
    n_neurons = sum(idx_curr);

    % Use pre-calculated PSTHs
    psth_CS_curr = psthZ_full{1}(idx_curr, :);
    psth_US_curr = psthZ_full{2}(idx_curr, :);

    % Calculate peak responses and latencies in ROI (0 to 0.5s)
    CS_peak = max(psth_CS_curr(:, g.roi), [], 2);
    US_peak = max(psth_US_curr(:, g.roi), [], 2);
    CS_min = min(psth_CS_curr(:, g.roi), [], 2);
    US_min = min(psth_US_curr(:, g.roi), [], 2);

    % Compute onset and offset latencies for manual clustering
    CS_onset_lat = nan(n_neurons, 1);
    CS_offset_lat = nan(n_neurons, 1);
    US_onset_lat = nan(n_neurons, 1);
    US_offset_lat = nan(n_neurons, 1);

    for n = 1:n_neurons
        [CS_onset_lat(n), CS_offset_lat(n)] = compute_onset_offset_latency(psth_CS_curr(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [US_onset_lat(n), US_offset_lat(n)] = compute_onset_offset_latency(psth_US_curr(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
    end

    % Classify neurons
    CS_excited = CS_peak >= g.excitation_threshold;
    US_excited = US_peak >= g.excitation_threshold;
    CS_inhibited = CS_min <= g.inhibition_threshold;
    US_inhibited = US_min <= g.inhibition_threshold;

    % Initialize cluster assignments
    Clusters_manual = zeros(n_neurons, 1);

    % Cluster 1: CS-selective (CS excited, US not excited)
    idx_CS_selective = CS_excited & ~US_excited;
    Clusters_manual(idx_CS_selective) = 1;

    % Cluster 2: US-selective (US excited, CS not excited)
    idx_US_selective = US_excited & ~CS_excited;
    Clusters_manual(idx_US_selective) = 2;

    % Cluster 3: Multisensory (both CS and US excited)
    idx_multisensory = CS_excited & US_excited;
    Clusters_manual(idx_multisensory) = 3;

    % Cluster 4: Inhibited (not excited by either, but inhibited by at least one)
    idx_inhibited = ~CS_excited & ~US_excited & (CS_inhibited | US_inhibited);
    Clusters_manual(idx_inhibited) = 4;

    % Cluster 5: Non-responsive (neither excited nor inhibited)
    idx_nonresponsive = ~CS_excited & ~US_excited & ~CS_inhibited & ~US_inhibited;
    Clusters_manual(idx_nonresponsive) = 5;

    % Sort within each cluster
    leafOrder_manual = [];

    for c = 1:5
        clust_idx = find(Clusters_manual == c);
        if ~isempty(clust_idx)
            % Sort by onset latency (primary), offset latency (secondary), then peak response
            if c == 1  % CS-selective: sort by CS latencies and peak
                onset_c = CS_onset_lat(clust_idx);
                offset_c = CS_offset_lat(clust_idx);
                peak_c = CS_peak(clust_idx);
                sort_matrix = [isnan(onset_c), onset_c, offset_c, -peak_c];  % negative peak for descending
                [~, sort_idx] = sortrows(sort_matrix, [1 2 3 4]);
            elseif c == 2  % US-selective: sort by US latencies and peak
                onset_c = US_onset_lat(clust_idx);
                offset_c = US_offset_lat(clust_idx);
                peak_c = US_peak(clust_idx);
                sort_matrix = [isnan(onset_c), onset_c, offset_c, -peak_c];
                [~, sort_idx] = sortrows(sort_matrix, [1 2 3 4]);
            elseif c == 3  % Multisensory: sort by earliest onset, then offset, then combined peak
                onset_c = min(CS_onset_lat(clust_idx), US_onset_lat(clust_idx));
                offset_c = min(CS_offset_lat(clust_idx), US_offset_lat(clust_idx));
                combined_peak = CS_peak(clust_idx) + US_peak(clust_idx);
                sort_matrix = [isnan(onset_c), onset_c, offset_c, -combined_peak];
                [~, sort_idx] = sortrows(sort_matrix, [1 2 3 4]);
            elseif c == 4  % Inhibited: sort by minimum (most inhibited first)
                combined_min = min(CS_min(clust_idx), US_min(clust_idx));
                [~, sort_idx] = sort(combined_min, 'ascend');
            else  % Non-responsive: sort by firing rate
                % Get original indices in full cell_metrics
                full_idx = find(idx_curr);
                fr = g.cell_metrics.firingRate(full_idx(clust_idx));
                [~, sort_idx] = sort(fr, 'descend');
            end
            leafOrder_manual = [leafOrder_manual; clust_idx(sort_idx)]; %#ok<AGROW>

            fprintf('  Cluster %d: %d units\n', c, numel(clust_idx));
        end
    end

    % Store results
    results_manual(ct).leafOrder = leafOrder_manual;
    results_manual(ct).Clusters = Clusters_manual;
    results_manual(ct).idx_curr = idx_curr;
    results_manual(ct).n_neurons = n_neurons;
end

%% Calculate matrices for manual clustering figure using pre-calculated PSTHs

fprintf('\nPreparing manual clustering heatmap matrices...\n');
matrices_manual = cell(1, numel(ttl));

for hmp = 1:numel(ttl)
    % Use pre-calculated PSTH
    psth_spx = psthZ_full{hmp};

    % Extract and sort PNs
    psth_PN = psth_spx(results_manual(1).idx_curr, :);
    psth_PN_sorted = psth_PN(results_manual(1).leafOrder, :);
    matrix_PN = psth_PN_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);

    % Extract and sort INs
    psth_IN = psth_spx(results_manual(2).idx_curr, :);
    psth_IN_sorted = psth_IN(results_manual(2).leafOrder, :);
    matrix_IN = psth_IN_sorted(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);

    % Combine: PNs on top, INs on bottom
    matrix_combined = [matrix_PN; matrix_IN];
    matrices_manual{hmp} = matrix_combined;
end

% Determine color limits (use same as k-means for consistency)
all_values_manual = [matrices_manual{:}];
if g.use_percentile
    clim_min_manual = prctile(all_values_manual(:), (100 - g.clim_percentile) / 2);
    clim_max_manual = prctile(all_values_manual(:), 100 - (100 - g.clim_percentile) / 2);
else
    clim_min_manual = min(all_values_manual(:));
    clim_max_manual = max(all_values_manual(:));
end

%% Plot manual clustering heatmap

fprintf('Plotting manual clustering figure...\n');
fig2 = figure('Position', [500, 100, 800, 800]);
t2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

cluster_colors_manual = [
    1 0 0;      % Red: CS-selective
    0 0 1;      % Blue: US-selective
    0.5 0 0.5;  % Purple: Multisensory
    0 0.7 0;    % Green: Inhibited
    0.5 0.5 0.5 % Gray: Non-responsive
];

for hmp = 1:numel(ttl)
    matrix_combined = matrices_manual{hmp};

    ax = nexttile(t2, hmp);
    imagesc(g.timeaxis_hmp, 1:size(matrix_combined, 1), matrix_combined);
    clim([clim_min_manual clim_max_manual]);
    colormap(g.colors.Heatmap);
    xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha', 1);

    if hmp == 1
        ylabel('Cell number', 'FontSize', g.fontSize2);
    else
        set(gca, 'YTickLabel', []);
    end

    xlabel('Time from Stimulus Onset (s)', 'FontSize', g.fontSize2);
    set(gca, 'FontSize', g.fontSize2);
    title(hmptitles{hmp}, 'FontSize', g.fontSize1);

    hold on;

    % Draw lines separating clusters within PNs
    Clusters_PN = results_manual(1).Clusters(results_manual(1).leafOrder);
    n_clu_PN = find(diff(Clusters_PN) ~= 0);
    for i = 1:length(n_clu_PN)
        cluster_before = Clusters_PN(n_clu_PN(i));
        yline(n_clu_PN(i) + 0.5, 'Color', cluster_colors_manual(cluster_before, :), 'LineWidth', 1);
    end

    % Draw line separating PNs from INs
    n_PN = results_manual(1).n_neurons;
    yline(n_PN + 0.5, 'Color', 'k', 'LineWidth', 3);

    % Draw lines separating clusters within INs
    Clusters_IN = results_manual(2).Clusters(results_manual(2).leafOrder);
    n_clu_IN = find(diff(Clusters_IN) ~= 0);
    for i = 1:length(n_clu_IN)
        cluster_before = Clusters_IN(n_clu_IN(i));
        yline(n_clu_IN(i) + n_PN + 0.5, 'Color', cluster_colors_manual(cluster_before, :), 'LineWidth', 1);
    end

    hold off;
end

% Add colorbar
cb2 = colorbar(nexttile(t2, 2), 'eastoutside', 'FontSize', g.fontSize2);
ylabel(cb2, 'Z-score Firing Rate', 'FontSize', g.fontSize2);

title(t2, sprintf('%s: Manual Clustering - PNs (top) and INs (bottom)\nCS-sel | US-sel | Multi | Inhibited | Non-resp', g.bR), 'FontSize', g.fontSize1 + 2);

fprintf('\nDone.\n');

%% Local function: onset and offset latency detector
function [onset_lat, offset_lat] = compute_onset_offset_latency(z_trace, event_inds, threshold, min_consec, bin_time)
    seg = z_trace(event_inds);
    if any(isnan(seg))
        onset_lat = NaN;
        offset_lat = NaN;
        return
    end

    % Detect onset
    isAbove = seg >= threshold;
    winsum = conv(double(isAbove), ones(1, min_consec), 'same');
    onset_idx = find(winsum >= min_consec, 1, 'first');

    if isempty(onset_idx)
        onset_lat = NaN;
        offset_lat = NaN;
    else
        onset_lat = (onset_idx - 1) * bin_time;

        % Detect offset: first time after onset when activity goes below threshold
        seg_after_onset = seg(onset_idx:end);
        isBelow = seg_after_onset < threshold;
        winsum_below = conv(double(isBelow), ones(1, min_consec), 'same');
        offset_idx_relative = find(winsum_below >= min_consec, 1, 'first');

        if isempty(offset_idx_relative)
            offset_lat = NaN;  % No offset detected (activity stays above threshold)
        else
            offset_idx = onset_idx + offset_idx_relative - 1;
            offset_lat = (offset_idx - 1) * bin_time;
        end
    end
end
