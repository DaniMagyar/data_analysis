% comparison_heatmap_005
% Third figure: Sort neurons within each cluster by peak latency
% Manual clustering (CS-selective, US-selective, Multisensory, Inhibited, Non-responsive)

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
g.post_time = 2;
g.test_time = 0.5;

% Clustering Parameters
g.clustnum = 5;
g.alpha = 0.8;  % Weight for duration in rank score

% Manual clustering thresholds
g.excitation_threshold = 2;   % Z-score threshold for excitatory response
g.inhibition_threshold = -2;  % Z-score threshold for inhibitory response (not used)
g.inhibition_fr_drop = 0.30;  % Firing rate drop threshold (30%)

% Colormap limits: use percentile to avoid extreme values
g.use_percentile = true;
g.clim_percentile = 95;

% Onset Detection Parameters
g.testvalue = g.excitation_threshold;
g.onset_threshold = g.testvalue;
g.bin_time = 0.001;
g.min_consec_bins = max(1, round(0.020 / g.bin_time));

g.smoothvalue = 101;
g.plotwin = [1 1];
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);

% ROI for clustering
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;

g.bR = 'Astria';

% Get PN and IN indices
idx_PN = strcmp(g.cell_metrics.brainRegion, g.bR) & strcmp(g.cell_metrics.putativeCellType, 'PN');
idx_IN = strcmp(g.cell_metrics.brainRegion, g.bR) & strcmp(g.cell_metrics.putativeCellType, 'IN');

%% CALCULATE ALL PSTHs ONCE AT THE BEGINNING

fprintf('Calculating PSTHs once for all stimuli...\n');
psthZ_full = cell(1, numel(ttl));
psthHz_full = cell(1, numel(ttl));

for hmp = 1:numel(ttl)
    fprintf('  Processing %s...\n', hmptitles{hmp});
    psth_spx_og = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Store raw Hz
    psthHz_full{hmp} = psth_spx_og;

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

%% MANUAL CLUSTERING using pre-calculated PSTHs

fprintf('\n=== Manual Clustering ===\n');
results_manual = struct();

cell_types = {'PN', 'IN'};
idx_celltypes = {idx_PN, idx_IN};

for ct = 1:2
    fprintf('Processing %s...\n', cell_types{ct});

    idx_curr = idx_celltypes{ct};
    n_neurons = sum(idx_curr);

    % Use pre-calculated PSTHs
    psth_CS_curr = psthZ_full{1}(idx_curr, :);
    psth_US_curr = psthZ_full{2}(idx_curr, :);
    psth_CS_Hz = psthHz_full{1}(idx_curr, :);
    psth_US_Hz = psthHz_full{2}(idx_curr, :);

    % Calculate peak responses and latencies in ROI (0 to 0.5s)
    CS_peak = max(psth_CS_curr(:, g.roi), [], 2);
    US_peak = max(psth_US_curr(:, g.roi), [], 2);

    % Calculate firing rate drop for inhibition detection (using Hz)
    baseline_idx = 1:(g.pre_time / g.bin_time);
    CS_baseline_fr = mean(psth_CS_Hz(:, baseline_idx), 2);
    US_baseline_fr = mean(psth_US_Hz(:, baseline_idx), 2);
    CS_test_fr = mean(psth_CS_Hz(:, g.roi), 2);
    US_test_fr = mean(psth_US_Hz(:, g.roi), 2);
    CS_fr_drop = (CS_baseline_fr - CS_test_fr) ./ (CS_baseline_fr + eps);
    US_fr_drop = (US_baseline_fr - US_test_fr) ./ (US_baseline_fr + eps);

    % Compute onset and offset latencies for manual clustering
    CS_onset_lat = nan(n_neurons, 1);
    CS_offset_lat = nan(n_neurons, 1);
    US_onset_lat = nan(n_neurons, 1);
    US_offset_lat = nan(n_neurons, 1);
    CS_peak_lat = nan(n_neurons, 1);
    US_peak_lat = nan(n_neurons, 1);

    for n = 1:n_neurons
        [CS_onset_lat(n), CS_offset_lat(n)] = compute_onset_offset_latency(psth_CS_curr(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);
        [US_onset_lat(n), US_offset_lat(n)] = compute_onset_offset_latency(psth_US_curr(n, :), g.roi, g.onset_threshold, g.min_consec_bins, g.bin_time);

        % Compute peak latency (time of maximum response in ROI)
        [~, CS_peak_idx] = max(psth_CS_curr(n, g.roi));
        [~, US_peak_idx] = max(psth_US_curr(n, g.roi));
        CS_peak_lat(n) = (CS_peak_idx - 1) * g.bin_time;
        US_peak_lat(n) = (US_peak_idx - 1) * g.bin_time;
    end

    % Initialize cluster assignments
    Clusters_manual = ones(n_neurons, 1);  % All INs in cluster 1

    % For PNs, classify into clusters
    if ct == 1
        % Classify neurons
        CS_excited = CS_peak >= g.excitation_threshold;
        US_excited = US_peak >= g.excitation_threshold;
        CS_inhibited = CS_fr_drop >= g.inhibition_fr_drop;
        US_inhibited = US_fr_drop >= g.inhibition_fr_drop;

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
    end

    % Sort within each cluster by rank score: onset + alpha * duration
    leafOrder_manual = [];

    if ct == 2  % INs: separate responsive/non-responsive, sort responsives
        % Classify as responsive or non-responsive
        CS_excited = CS_peak >= g.excitation_threshold;
        US_excited = US_peak >= g.excitation_threshold;
        CS_inhibited = CS_fr_drop >= g.inhibition_fr_drop;
        US_inhibited = US_fr_drop >= g.inhibition_fr_drop;

        is_responsive = CS_excited | US_excited | CS_inhibited | US_inhibited;
        Clusters_manual(is_responsive) = 1;  % Responsive
        Clusters_manual(~is_responsive) = 2;  % Non-responsive

        % Sort responsive INs by rank score
        resp_idx = find(is_responsive);
        if ~isempty(resp_idx)
            CS_duration = CS_offset_lat(resp_idx) - CS_onset_lat(resp_idx);
            US_duration = US_offset_lat(resp_idx) - US_onset_lat(resp_idx);
            CS_rank_score = CS_onset_lat(resp_idx) + g.alpha * CS_duration;
            US_rank_score = US_onset_lat(resp_idx) + g.alpha * US_duration;
            rank_score = mean([CS_rank_score, US_rank_score], 2, 'omitnan');

            sort_matrix = [isnan(rank_score), rank_score];
            [~, sort_idx] = sortrows(sort_matrix, [1 2]);
            leafOrder_manual = resp_idx(sort_idx);
            fprintf('  INs responsive: %d units, %d with valid rank score\n', numel(resp_idx), sum(~isnan(rank_score)));
        end

        % Sort non-responsive INs by firing rate
        nonresp_idx = find(~is_responsive);
        if ~isempty(nonresp_idx)
            full_idx = find(idx_curr);
            fr = g.cell_metrics.firingRate(full_idx(nonresp_idx));
            [~, sort_idx] = sort(fr, 'descend');
            leafOrder_manual = [leafOrder_manual; nonresp_idx(sort_idx)];
            fprintf('  INs non-responsive: %d units\n', numel(nonresp_idx));
        end
    else  % PNs: cluster and sort
        % Process clusters in order: 1, 2, 3, 5, 4 (inhibited last)
        cluster_order = [1, 2, 3, 5, 4];
        for c = cluster_order
            clust_idx = find(Clusters_manual == c);
            if ~isempty(clust_idx)
                % Determine sorting criteria based on cluster type
                if c == 1  % CS-selective: use CS latencies
                    onset_c = CS_onset_lat(clust_idx);
                    offset_c = CS_offset_lat(clust_idx);
                elseif c == 2  % US-selective: use US latencies
                    onset_c = US_onset_lat(clust_idx);
                    offset_c = US_offset_lat(clust_idx);
                elseif c == 3  % Multisensory: calculate rank score for each, then average
                    CS_duration = CS_offset_lat(clust_idx) - CS_onset_lat(clust_idx);
                    US_duration = US_offset_lat(clust_idx) - US_onset_lat(clust_idx);
                    CS_rank_score = CS_onset_lat(clust_idx) + g.alpha * CS_duration;
                    US_rank_score = US_onset_lat(clust_idx) + g.alpha * US_duration;
                    rank_score = mean([CS_rank_score, US_rank_score], 2, 'omitnan');

                    % Sort by rank score (ascending)
                    sort_matrix = [isnan(rank_score), rank_score];
                    [~, sort_idx] = sortrows(sort_matrix, [1 2]);
                    leafOrder_manual = [leafOrder_manual; clust_idx(sort_idx)]; %#ok<AGROW>
                    fprintf('  Cluster %d: %d units, %d with valid rank score\n', c, numel(clust_idx), sum(~isnan(rank_score)));
                    continue;
                elseif c == 4  % Inhibited: sort by mean z-score during test time
                    mean_zscore = mean([mean(psth_CS_curr(clust_idx, g.roi), 2), mean(psth_US_curr(clust_idx, g.roi), 2)], 2);
                    [~, sort_idx] = sort(mean_zscore, 'descend');
                    leafOrder_manual = [leafOrder_manual; clust_idx(sort_idx)]; %#ok<AGROW>
                    fprintf('  Cluster %d (Inhibited): %d units\n', c, numel(clust_idx));
                    continue;
                else  % Non-responsive: sort by mean z-score during test time
                    mean_zscore = mean([mean(psth_CS_curr(clust_idx, g.roi), 2), mean(psth_US_curr(clust_idx, g.roi), 2)], 2);
                    [~, sort_idx] = sort(mean_zscore, 'descend');
                    leafOrder_manual = [leafOrder_manual; clust_idx(sort_idx)]; %#ok<AGROW>
                    fprintf('  Cluster %d (Non-responsive): %d units\n', c, numel(clust_idx));
                    continue;
                end

                % Calculate rank score: onset + alpha * duration
                duration_c = offset_c - onset_c;
                rank_score = onset_c + g.alpha * duration_c;

                % Sort by rank score (ascending)
                sort_matrix = [isnan(rank_score), rank_score];
                [~, sort_idx] = sortrows(sort_matrix, [1 2]);

                leafOrder_manual = [leafOrder_manual; clust_idx(sort_idx)]; %#ok<AGROW>

                fprintf('  Cluster %d: %d units, %d with valid rank score\n', c, numel(clust_idx), sum(~isnan(rank_score)));
            end
        end
    end

    % Store results
    results_manual(ct).leafOrder = leafOrder_manual;
    results_manual(ct).Clusters = Clusters_manual;
    results_manual(ct).idx_curr = idx_curr;
    results_manual(ct).n_neurons = n_neurons;
end

%% Calculate matrices for peak latency-sorted manual clustering figure using pre-calculated PSTHs

fprintf('\nPreparing peak latency-sorted manual clustering heatmap matrices...\n');
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

% Determine color limits (use same as 004 for consistency)
all_values_manual = [matrices_manual{:}];
if g.use_percentile
    clim_min_manual = prctile(all_values_manual(:), (100 - g.clim_percentile) / 2);
    clim_max_manual = prctile(all_values_manual(:), 100 - (100 - g.clim_percentile) / 2);
else
    clim_min_manual = min(all_values_manual(:));
    clim_max_manual = max(all_values_manual(:));
end

%% Plot peak latency-sorted manual clustering heatmap

fprintf('Plotting peak latency-sorted manual clustering figure...\n');
fig = figure('Position', [500, 100, 800, 800]);
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

cluster_colors_manual = [
    0.8 0.2 0.2;    % Soft red: CS-selective
    0.2 0.4 0.8;    % Soft blue: US-selective
    0.6 0.2 0.6;    % Soft purple: Multisensory
    0.2 0.6 0.3;    % Soft green: Inhibited
    0.6 0.6 0.6     % Medium gray: Non-responsive
];

for hmp = 1:numel(ttl)
    matrix_combined = matrices_manual{hmp};

    ax = nexttile(t, hmp);
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
        yline(n_clu_PN(i) + 0.5, 'Color', cluster_colors_manual(cluster_before, :), 'LineWidth', 2);
    end

    % Draw line separating PNs from INs
    n_PN = results_manual(1).n_neurons;
    yline(n_PN + 0.5, 'Color', 'k', 'LineWidth', 4);

    % Draw line separating responsive/non-responsive INs
    Clusters_IN = results_manual(2).Clusters(results_manual(2).leafOrder);
    n_clu_IN = find(diff(Clusters_IN) ~= 0);
    if ~isempty(n_clu_IN)
        yline(n_clu_IN(1) + n_PN + 0.5, 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
    end

    hold off;
end

% Add colorbar
cb = colorbar(nexttile(t, 2), 'eastoutside', 'FontSize', g.fontSize2);
ylabel(cb, 'Z-score Firing Rate', 'FontSize', g.fontSize2);

title(t, sprintf('%s: Rank Score-Sorted Manual Clustering (alpha=%.1f) - PNs (top) and INs (bottom)\nCS-sel | US-sel | Multi | Inhibited | Non-resp', g.bR, g.alpha), 'FontSize', g.fontSize1 + 2);

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

