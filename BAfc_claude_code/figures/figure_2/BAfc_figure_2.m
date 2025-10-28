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
g.fontSize1 = 14;
g.fontSize2 = 12;
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
brain_regions = {'LA', 'BA', 'Astria', 'CeA'};
cell_type_filter = {'all', 'all', 'all', 'all'};  % LA/BA: PN only, Astria/CeA: all

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
    CS_test_fr = mean(psth_CS_Hz(:, g.roi), 2);
    US_test_fr = mean(psth_US_Hz(:, g.roi), 2);
    CS_fr_drop = (CS_baseline_fr - CS_test_fr) ./ (CS_baseline_fr + eps);
    US_fr_drop = (US_baseline_fr - US_test_fr) ./ (US_baseline_fr + eps);

    % Classify neurons
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

    % Store results
    results_all{br}.Clusters = Clusters;
    results_all{br}.leafOrder = leafOrder;
    results_all{br}.psth_CS = psth_CS;
    results_all{br}.psth_US = psth_US;
    results_all{br}.psth_CS_Hz = psth_CS_Hz;
    results_all{br}.psth_US_Hz = psth_US_Hz;
    results_all{br}.n_neurons = n_neurons;
end

%% Create figure
fig = figure('Position', [-100, -100, 1500, 3000], 'Units', 'pixels');
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
    if br == 1
        ylabel('Neuron #', 'FontSize', g.fontSize2);
    end

    % Set yticks to first and last
    n_neurons_plot = size(matrix, 1);
    yticks([1, n_neurons_plot]);

    if br == 1
        title('CS', 'FontSize', g.fontSize1);
    end
    if br == 4
        xlabel('Time (s)', 'FontSize', g.fontSize2);
        xticks([-1 0 1]);
    else
        set(gca, 'XTickLabel', []);
    end

    % Add cluster lines (black)
    hold on;
    Clusters_sorted = res.Clusters(res.leafOrder);
    n_clu = find(diff(Clusters_sorted) ~= 0);
    for i = 1:length(n_clu)
        yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
    end
    hold off;

    % Add brain region label
    text(-0.9, res.n_neurons/2, brain_regions{br}, 'FontSize', g.fontSize1+2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Rotation', 90);

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

    % Add cluster lines (black)
    hold on;
    for i = 1:length(n_clu)
        yline(n_clu(i) + 0.5, 'Color', 'k', 'LineWidth', 1);
    end
    hold off;

    % Stacked bar chart - column 3 (middle, half width)
    ax = nexttile(t, (br-1)*5 + 3);

    % Count neurons in each cluster (in order they appear in heatmap)
    cluster_order = [1, 2, 3, 4, 5];  % Same as heatmap top to bottom
    cluster_counts = zeros(5, 1);
    for i = 1:5
        c = cluster_order(i);
        cluster_counts(i) = sum(res.Clusters == c);
    end

    % Calculate proportions
    cluster_props = cluster_counts / res.n_neurons * 100;

    % Create vertical stacked bar (following Gergo style)
    % Reverse order: bar stacks bottom-up, we want cluster 5 at bottom, cluster 1 at top
    b = bar(0.5, flipud(cluster_props)', 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.5);

    % Color each segment (reversed order to match)
    for i = 1:5
        c = cluster_order(6-i);  % Map: i=1→c=5, i=2→c=4, i=3→c=3, i=4→c=2, i=5→c=1
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
    for i = 1:5
        c = cluster_order(6-i);  % Reversed cluster index
        prop = cluster_props(c);
        if round(prop) > 0
            text(0.5, cumulative + prop/2, sprintf('%.0f%%', prop), ...
                'HorizontalAlignment', 'center', 'FontSize', g.fontSize2-2, 'FontWeight', 'bold', 'Color', 'w');
        end
        cumulative = cumulative + prop;
    end

    % Add colored text labels as legend (all rows, matching lineplot positions)
    % Position text labels to the right of bar
    y_positions = [85, 70, 50, 30, 10];  % Approximate vertical positions for each cluster
    for i = 1:5
        c = cluster_order(i);
        % Only show label if cluster exists in this brain region
        if sum(res.Clusters == c) > 0
            text(1.3, y_positions(i), cluster_names{c}, ...
                'FontSize', g.fontSize2, 'FontWeight', 'bold', ...
                'Color', cluster_colors(c, :), 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle');
        end
    end

    % Add total n below legend
    text(0.5, 0.05, sprintf('(n=%d)', res.n_neurons), 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', g.fontSize2-1, 'FontWeight', 'normal');

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
    % Separate limits for responsive clusters (1-4) and inhibited (5)
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

    for c = [1 2 3 4 5]
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
                'FontSize', g.fontSize2-2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            hold off;
        end

        if c == 5  % Add separate scalebar to cluster 5 (inhibited) with different scale
            curr_ylim = ylim;
            y_range = curr_ylim(2) - curr_ylim(1);
            scalebar_size = 5;  % Smaller scalebar for inhibited neurons
            x_pos = time_vec(end) - 0.3;
            y_pos = curr_ylim(2) - scalebar_size - 0.05*y_range;  % Position near top

            hold on;
            plot([x_pos x_pos], [y_pos y_pos+scalebar_size], 'k-', 'LineWidth', 3);
            text(x_pos+0.2, y_pos+scalebar_size/2, sprintf('%d Hz', scalebar_size), ...
                'FontSize', g.fontSize2-2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            hold off;
        end
    end
end

% Add colorbar - manually create with full control
drawnow;  % Ensure all positions are updated
cb_width = 0.010;
cb_left = 0.025;
cb_bottom = 0.05;
cb_height = 0.15;
cb_ax = axes('Position', [cb_left, cb_bottom, cb_width, cb_height]);
imagesc(cb_ax, [0 1], [clim_min clim_max], repmat(linspace(clim_min, clim_max, 256)', 1, 10));
colormap(cb_ax, g.colors.Heatmap);
set(cb_ax, 'YDir', 'normal');
set(cb_ax, 'XTick', [], 'YAxisLocation', 'left');
ylabel(cb_ax, 'Z-score', 'FontSize', g.fontSize2);
cb_ax.FontSize = g.fontSize2;

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
