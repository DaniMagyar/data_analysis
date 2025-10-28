% BAfc_figure_1
% Final publication figure with all panels organized as separate functions
% TO CREATE PRINT QUALITY: exportgraphics(gcf, 'figure.tiff', 'Resolution', 300);

clear all

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

g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'});
g.cell_metrics = BAfc_match_celltypes('cell_metrics', g.cell_metrics);

clearvars -except g
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 14;  % Title font size
g.fontSize2 = 12;  % Axis font size

%% Initialize figure (A4 dimensions: 210mm x 297mm = 8.27" x 11.69" at 96 DPI)
fig = figure('Position', [0, 0, 1500, 1500], 'Units', 'pixels');  % A4 size
t = tiledlayout(fig, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

%% ROW 1: Experimental setup and example traces
panel_experimental_setup(t, g, 1, [1 2]);  % Spans 2 columns
panel_example_traces(t, g, 3, 'Example 1');
panel_example_traces(t, g, 4, 'Example 2');  % Duplicate for now

%% ROW 2: ISI examples, spike features, and waveforms
panel_example_ISI(t, g, 5, 'PN');
panel_example_ISI(t, g, 6, 'IN');
panel_spike_features(t, g, 7);
panel_normalized_waveforms(t, g, 8);

%% ROW 3: Firing rate distributions
panel_firing_rate_distribution(t, g, 9, 'LA');
panel_firing_rate_distribution(t, g, 10, 'BA');
panel_firing_rate_distribution(t, g, 11, 'Astria');
panel_firing_rate_distribution(t, g, 12, 'CeA');

%% ROW 4: Burst index distributions
panel_burst_index_distribution(t, g, 13, 'LA');
panel_burst_index_distribution(t, g, 14, 'BA');
panel_burst_index_distribution(t, g, 15, 'Astria');
panel_burst_index_distribution(t, g, 16, 'CeA');

%% ========================================================================
%% PANEL FUNCTIONS
%% ========================================================================

%% Row 1 Functions
function panel_experimental_setup(t, g, tile_num, tile_span)
    % Panel: Experimental setup schematic
    ax = nexttile(t, tile_num, tile_span);
    [img, cmap] = imread([g.mainFolder '\drawed_mouse.png']);
    if ~isempty(cmap)
        img = ind2rgb(img, cmap);
    end
    imshow(img, 'Parent', ax);
    title(ax, 'Experimental setup', 'FontSize', g.fontSize1);
end

function panel_example_traces(t, g, tile_num, example_name)
    % Panel: Example raw traces and PSTH
    % Raw traces
    rawMD292_002 = load([g.mainFolder '\MD292_002_kilosort\kilosort25preprocess\temp_wh_ch32.mat']);
    rawMD292_002 = rawMD292_002.temp_wh_Ch32.values * 0.195; % convert to mV
    cellID = intersect(find(strcmp(g.cell_metrics.animal, 'MD292_002')), find(g.cell_metrics.cluID == 341));
    ttlMD292_002 = load([g.mainFolder '\MD292_002_kilosort\kilosort25preprocess\TTLsKS.mat']);
    ttlMD292_002 = ttlMD292_002.shocks([5 7 36 41 43 47]);
    twin = [-0.01 0.04];
    fs = 30000;

    % Subplot 1: Raw traces
    ax1 = nexttile(t, tile_num, [1 1]);
    hold(ax1, 'on');
    for ii = 1:size(ttlMD292_002, 1)
        sn_curr = ttlMD292_002(ii, 1) * fs;
        timeaxis = linspace(twin(1), twin(2), (twin(2)-twin(1))*fs) * 1000;
        plot(ax1, timeaxis, rawMD292_002(round(sn_curr+twin(1)*fs):round(sn_curr+twin(2)*fs-1)), 'Color', [0 0 0]);
    end
    ylabel(ax1, 'Voltage (mV)', 'FontSize', g.fontSize2);
    ylim(ax1, [-0.5 0.3]);
    yticks(ax1, [-0.5 -0.1 0.3]);
    set(ax1, 'FontSize', g.fontSize2);
    ax1.XColor = 'none';
    yLimits = get(ax1, 'YLim');
    fill(ax1, [0 10 10 0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ax1.Layer = 'top';
    ax1.Box = 'off';
    title(ax1, example_name, 'FontSize', g.fontSize1);
end

%% Row 2 Functions
function panel_example_ISI(t, g, tile_num, cell_type)
    % Panel: Example ISI with ACG and waveform insets
    ax = nexttile(t, tile_num);

    if strcmp(cell_type, 'PN')
        cellID = intersect(find(strcmp(g.cell_metrics.animal, 'MD292_002')), find(g.cell_metrics.cluID == 479));
        color_primary = g.colors.PN_primary;
        title_text = 'Example PN ISI';
    else % IN
        cellID = intersect(find(strcmp(g.cell_metrics.animal, 'MD295_001')), find(g.cell_metrics.cluID == 106));
        color_primary = g.colors.IN_primary;
        title_text = 'Example IN ISI';
    end

    spike_timestamps = g.cell_metrics.spikes.times{cellID};
    ISI = diff(spike_timestamps);
    num_bins = 50;
    bin_edges = logspace(log10(min(ISI)), log10(max(ISI)), num_bins);

    histogram(ax, ISI, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', color_primary);
    xlim(ax, [0.001 80]);
    set(ax, 'XScale', 'log');
    xticks(ax, [0.001 0.1 10]);
    xlabel(ax, 'Time (ms)', 'FontSize', g.fontSize2);
    ylim(ax, [0 0.1]);
    ylabel(ax, 'Probability', 'FontSize', g.fontSize2);
    yticks(ax, [0 0.05 0.1]);
    set(ax, 'FontSize', g.fontSize2);
    title(ax, title_text, 'FontSize', g.fontSize1);
    ax.Box = 'off';

    % ACG inset (top right)
    p = ax.Position;
    acg_ax = axes('Position', [p(1)+p(3)*0.65 p(2)+p(4)*0.65 p(3)*0.3 p(4)*0.3]);
    bar(acg_ax, g.cell_metrics.acg.narrow(:, cellID), 'FaceColor', color_primary, 'EdgeColor', color_primary);
    hold(acg_ax, 'on');
    ylim_acg = ylim(acg_ax);
    y = ylim_acg(1) - (ylim_acg(2)-ylim_acg(1))*0.1;
    plot(acg_ax, [160 200], [y y], 'k', 'LineWidth', 2);
    text(acg_ax, 180, y - (ylim_acg(2)-ylim_acg(1))*0.08, '20 ms', 'HorizontalAlignment', 'center', 'FontSize', g.fontSize2-2);
    acg_ax.XColor = 'none';
    acg_ax.YColor = 'none';

    % Waveform inset (top right, below ACG)
    wf_ax = axes('Position', [p(1)+p(3)*0.65 p(2)+p(4)*0.30 p(3)*0.3 p(4)*0.3]);
    plot(wf_ax, g.cell_metrics.waveforms.filt{cellID}, 'Color', color_primary, 'LineWidth', 2);
    hold(wf_ax, 'on');
    y = min(g.cell_metrics.waveforms.filt{cellID}) * 1.2;
    plot(wf_ax, [33 48], [y y], 'k', 'LineWidth', 2);
    text(wf_ax, 40, y*1.2, '0.5 ms', 'HorizontalAlignment', 'center', 'FontSize', g.fontSize2-2);
    wf_ax.XColor = 'none';
    wf_ax.YColor = 'none';
end

function panel_spike_features(t, g, tile_num)
    % Panel: Spike features across unit types
    ax = nexttile(t, tile_num);

    idx_PN = strcmp(g.cell_metrics.putativeCellType, 'PN') & (strcmp(g.cell_metrics.brainRegion, 'LA') | strcmp(g.cell_metrics.brainRegion, 'BA'));
    idx_IN = strcmp(g.cell_metrics.putativeCellType, 'IN') & (strcmp(g.cell_metrics.brainRegion, 'LA') | strcmp(g.cell_metrics.brainRegion, 'BA'));
    idx_unknown = strcmp(g.cell_metrics.putativeCellType, 'unknown') & (strcmp(g.cell_metrics.brainRegion, 'LA') | strcmp(g.cell_metrics.brainRegion, 'BA'));

    g.data_ttp = g.cell_metrics.spikes.ttp;
    g.data_firing_rate = g.cell_metrics.firingRate;

    hold(ax, 'on');
    plot(ax, g.data_ttp(idx_unknown), g.data_firing_rate(idx_unknown), 'o', 'Color', [.7 .7 .7], 'MarkerSize', 6);
    plot(ax, g.data_ttp(idx_IN), g.data_firing_rate(idx_IN), 'o', 'Color', g.colors.IN_primary, 'MarkerSize', 6);
    plot(ax, g.data_ttp(idx_PN), g.data_firing_rate(idx_PN), 'o', 'Color', g.colors.PN_primary, 'MarkerSize', 6);

    % Set limits first
    xlim(ax, [0.1 0.8]);
    ylim_curr = ylim(ax);

    % Draw partial dashed lines
    plot(ax, [0.1 0.4], [10 10], '--k', 'LineWidth', 1);  % Horizontal line from left to 0.4
    plot(ax, [0.4 0.4], [10 ylim_curr(2)], '--k', 'LineWidth', 1);  % Vertical line from 10 to top

    ylabel(ax, 'Firing rate (Hz)', 'FontSize', g.fontSize2);
    xlabel(ax, 'Trough to peak (ms)', 'FontSize', g.fontSize2);
    set(ax, 'FontSize', g.fontSize2);
    title(ax, 'Spike features', 'FontSize', g.fontSize1);
    xticks(ax, [0.1 0.45 0.8]);
    yticks(ax, [ylim_curr(1) mean(ylim_curr) ylim_curr(2)]);
    ax.Box = 'off';
end

function panel_normalized_waveforms(t, g, tile_num)
    % Panel: Normalized waveforms across unit types
    ax = nexttile(t, tile_num);

    idx_PN = strcmp(g.cell_metrics.putativeCellType, 'PN') & (strcmp(g.cell_metrics.brainRegion, 'LA') | strcmp(g.cell_metrics.brainRegion, 'BA'));
    idx_IN = strcmp(g.cell_metrics.putativeCellType, 'IN') & (strcmp(g.cell_metrics.brainRegion, 'LA') | strcmp(g.cell_metrics.brainRegion, 'BA'));

    timeaxis = g.cell_metrics.waveforms.time{1};
    waveforms = zscore(cell2mat(g.cell_metrics.waveforms.filt(:)), [], 2);
    waveforms(g.cell_metrics.polarity > 0, :) = -waveforms(g.cell_metrics.polarity > 0, :);

    % Remove outliers
    wide_idx = find(idx_PN == 1);
    mean_wide = mean(waveforms(wide_idx, :), 1);
    narrow_idx = find(idx_IN == 1);
    mean_narrow = mean(waveforms(narrow_idx, :), 1);

    distances = zeros(size(waveforms, 1), 1);
    for ii = wide_idx
        distances(ii) = sqrt(sum((waveforms(ii, :) - mean_wide).^2));
    end
    for ii = narrow_idx
        distances(ii) = sqrt(sum((waveforms(ii, :) - mean_narrow).^2));
    end
    waveforms(distances > 7, :) = NaN;

    % Compute mean and std
    mean_IN = mean(waveforms(idx_IN, 13:48), 1, 'omitnan');
    std_IN = std(waveforms(idx_IN, 13:48), 0, 1, 'omitnan');
    mean_PN = mean(waveforms(idx_PN, 13:48), 1, 'omitnan');
    std_PN = std(waveforms(idx_PN, 13:48), 0, 1, 'omitnan');

    hold(ax, 'on');
    fill(ax, [timeaxis(13:48), fliplr(timeaxis(13:48))], ...
        [mean_IN + std_IN, fliplr(mean_IN - std_IN)], ...
        g.colors.IN_secondary, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(ax, timeaxis(13:48), mean_IN, 'Color', g.colors.IN_primary, 'LineWidth', 3);
    fill(ax, [timeaxis(13:48), fliplr(timeaxis(13:48))], ...
        [mean_PN + std_PN, fliplr(mean_PN - std_PN)], ...
        g.colors.PN_secondary, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(ax, timeaxis(13:48), mean_PN, 'Color', g.colors.PN_primary, 'LineWidth', 3);
    text(ax, 0.2, -2, ['PN (n=' num2str(sum(idx_PN)) ')'], 'Color', g.colors.PN_primary, 'FontSize', g.fontSize2-1);
    text(ax, 0.2, -2.5, ['IN (n=' num2str(sum(idx_IN)) ')'], 'Color', g.colors.IN_primary, 'FontSize', g.fontSize2-1);
    xlabel(ax, 'Time (ms)', 'FontSize', g.fontSize2);
    ylabel(ax, 'Z-score', 'FontSize', g.fontSize2);
    xlim_curr = xlim(ax);
    xticks(ax, [xlim_curr(1) mean(xlim_curr) xlim_curr(2)]);
    yticks(ax, [-4 0 4]);
    set(ax, 'FontSize', g.fontSize2);
    title(ax, 'Normalized waveforms', 'FontSize', g.fontSize1);
    ax.Box = 'off';
end

%% Row 3 Functions
function panel_firing_rate_distribution(t, g, tile_num, brain_region)
    % Panel: Firing rate distribution by brain region
    ax = nexttile(t, tile_num);

    % Get indices for brain region
    idx_region = strcmp(g.cell_metrics.brainRegion, brain_region);

    if strcmp(brain_region, 'LA') || strcmp(brain_region, 'BA')
        % Separate PN and IN
        idx_PN = idx_region & strcmp(g.cell_metrics.putativeCellType, 'PN');
        idx_IN = idx_region & strcmp(g.cell_metrics.putativeCellType, 'IN');

        firing_rate_PN = g.cell_metrics.firingRate(idx_PN);
        firing_rate_IN = g.cell_metrics.firingRate(idx_IN);

        % Log-spaced bins
        bin_edges = logspace(log10(0.1), log10(30), 20);

        hold(ax, 'on');
        histogram(ax, firing_rate_PN, 'BinEdges', bin_edges, 'FaceColor', g.colors.PN_primary, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        histogram(ax, firing_rate_IN, 'BinEdges', bin_edges, 'FaceColor', g.colors.IN_primary, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        legend(ax, {['PN (n=' num2str(sum(idx_PN)) ')'], ['IN (n=' num2str(sum(idx_IN)) ')']}, 'FontSize', g.fontSize2-1, 'Location', 'northeast');
    else
        % All neurons together
        firing_rate = g.cell_metrics.firingRate(idx_region);
        bin_edges = logspace(log10(0.1), log10(30), 20);
        histogram(ax, firing_rate, 'BinEdges', bin_edges, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        text(ax, 0.6, 0.9, ['n=' num2str(sum(idx_region))], 'Units', 'normalized', 'FontSize', g.fontSize2-1);
    end

    xlabel(ax, 'Firing Rate (Hz)', 'FontSize', g.fontSize2);
    ylabel(ax, 'Count', 'FontSize', g.fontSize2);
    title(ax, [brain_region ' Firing Rate'], 'FontSize', g.fontSize1);
    xlim(ax, [0.1 30]);
    ylim(ax, [0 50]);
    set(ax, 'XScale', 'log');
    xticks(ax, [0.1 1 10]);
    yticks(ax, [0 25 50]);
    set(ax, 'FontSize', g.fontSize2);
    ax.Box = 'off';
end

%% Row 4 Functions
function panel_burst_index_distribution(t, g, tile_num, brain_region)
    % Panel: Burst index distribution by brain region
    ax = nexttile(t, tile_num);

    % Get indices for brain region
    idx_region = strcmp(g.cell_metrics.brainRegion, brain_region);

    if strcmp(brain_region, 'LA') || strcmp(brain_region, 'BA')
        % Separate PN and IN
        idx_PN = idx_region & strcmp(g.cell_metrics.putativeCellType, 'PN');
        idx_IN = idx_region & strcmp(g.cell_metrics.putativeCellType, 'IN');

        burst_index_PN = g.cell_metrics.burstIndex_Royer2012(idx_PN);
        burst_index_IN = g.cell_metrics.burstIndex_Royer2012(idx_IN);

        % Log-spaced bins
        bin_edges = logspace(log10(0.1), log10(100), 20);

        hold(ax, 'on');
        histogram(ax, burst_index_PN, 'BinEdges', bin_edges, 'FaceColor', g.colors.PN_primary, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        histogram(ax, burst_index_IN, 'BinEdges', bin_edges, 'FaceColor', g.colors.IN_primary, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        legend(ax, {['PN (n=' num2str(sum(idx_PN)) ')'], ['IN (n=' num2str(sum(idx_IN)) ')']}, 'FontSize', g.fontSize2-1, 'Location', 'northeast');
    else
        % All neurons together
        burst_index = g.cell_metrics.burstIndex_Royer2012(idx_region);
        bin_edges = logspace(log10(0.1), log10(100), 20);
        histogram(ax, burst_index, 'BinEdges', bin_edges, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        text(ax, 0.6, 0.9, ['n=' num2str(sum(idx_region))], 'Units', 'normalized', 'FontSize', g.fontSize2-1);
    end

    xlabel(ax, 'Burst Index', 'FontSize', g.fontSize2);
    ylabel(ax, 'Count', 'FontSize', g.fontSize2);
    title(ax, [brain_region ' Burst Index'], 'FontSize', g.fontSize1);
    xlim(ax, [0.1 100]);
    ylim(ax, [0 30]);
    set(ax, 'XScale', 'log');
    xticks(ax, [0.1 1 10]);
    yticks(ax, [0 15 30]);
    set(ax, 'FontSize', g.fontSize2);
    ax.Box = 'off';
end
