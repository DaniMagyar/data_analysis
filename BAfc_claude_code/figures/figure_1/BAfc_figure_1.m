% BAfc_figure_1
% Final publication figure with all panels organized as separate functions
% TO CREATE PRINT QUALITY: exportgraphics(gcf, 'figure.tiff', 'Resolution', 300);

close all; clear all

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
g.fontSize1 = 12;  % Title font size
g.fontSize2 = 10;  % Axis font size

%% Initialize figure (A4 dimensions: 210mm x 297mm = 8.27" x 11.69" at 96 DPI)
fig = figure('Position', [0, 0, 1000, 1000], 'Units', 'pixels');  % A4 size
t = tiledlayout(fig, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

%% ROW 1: Experimental setup and example traces
ax_A = panel_experimental_setup(t, g, 1, [1 2]);  % Spans 2 columns
ax_B1 = panel_example_traces(t, g, 3, 'Example CS response', 'MD312_001', 84);
add_panel_label(ax_B1, 'B');
ax_B2 = panel_example_traces(t, g, 4, 'Example US response', 'MD292_002', 32);

%% ROW 2: ISI examples, spike features, and waveforms
ax_C1 = panel_example_ISI(t, g, 5, 'PN');
add_panel_label(ax_C1, 'C');
ax_C2 = panel_example_ISI(t, g, 6, 'IN');
ax_D1 = panel_spike_features(t, g, 7);
add_panel_label(ax_D1, 'D');
ax_D2 = panel_normalized_waveforms(t, g, 8);

%% ROW 3: Firing rate distributions
t_fr = tiledlayout(t, 1, 4, 'TileSpacing', 'tight');
t_fr.Layout.Tile = 9;
t_fr.Layout.TileSpan = [1 4];
ax_E = panel_firing_rate_distribution(t_fr, g, 1, 'LA', true);
add_panel_label(ax_E, 'E');
panel_firing_rate_distribution(t_fr, g, 2, 'BA', false);
panel_firing_rate_distribution(t_fr, g, 3, 'Astria', false);
panel_firing_rate_distribution(t_fr, g, 4, 'CeA', false);

%% ROW 4: Burst index distributions
t_bi = tiledlayout(t, 1, 4, 'TileSpacing', 'tight');
t_bi.Layout.Tile = 13;
t_bi.Layout.TileSpan = [1 4];
ax_F = panel_burst_index_distribution(t_bi, g, 1, 'LA', true);
add_panel_label(ax_F, 'F');
panel_burst_index_distribution(t_bi, g, 2, 'BA', false);
panel_burst_index_distribution(t_bi, g, 3, 'Astria', false);
panel_burst_index_distribution(t_bi, g, 4, 'CeA', false);

%% Export statistics
export_figure1_to_excel(g);
exportgraphics(gcf, 'figure_1.png', 'Resolution', 300);
%% ========================================================================
%% PANEL FUNCTIONS
%% ========================================================================

%% Helper function for panel labels
function add_panel_label(ax, label)
    % Add panel label (A, B, C, etc.) outside top-left corner
    text(ax, -0.20, 1.15, label, 'Units', 'normalized', ...
        'FontSize', 14, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
end

%% Row 1 Functions
function ax = panel_experimental_setup(t, g, tile_num, tile_span)
    % Panel: Experimental setup schematic (placeholder - add image manually)
    ax = nexttile(t, tile_num, tile_span);
    xlim(ax, [0 1]);
    ylim(ax, [0 1]);
    set(ax, 'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', []);
    title(ax, 'Experimental setup', 'FontSize', g.fontSize1);
    box(ax, 'off');
    % Add label A directly
    text(ax, -0.10, 1.15, 'A', 'Units', 'normalized', ...
        'FontSize', 14, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
end

function ax = panel_example_traces(t, g, tile_num, example_name, animal, channel)
    % Panel: Example raw traces and PSTH
    twin = [-0.04 0.06];
    fs = 30000;

    % Load data based on example number
    if tile_num == 3
        % Example 1: MD312_001
        rawMD312_001 = load([g.mainFolder '\MD312_001_kilosort\kilosort25preprocess\temp_wh_ch_122.mat']);
        rawMD312_001 = rawMD312_001.temp_wh_Ch122.values * 0.195; % convert to mV
        ttlMD312_001 = load([g.mainFolder '\MD312_001_kilosort\kilosort25preprocess\TTLsKS.mat']);
        ttlMD312_001 = ttlMD312_001.triptest_sound_only([2 3 4 5 6], :);

        raw_data = rawMD312_001;
        ttl_data = ttlMD312_001;
    else
        % Example 2: MD292_002
        rawMD292_002 = load([g.mainFolder '\MD292_002_kilosort\kilosort25preprocess\temp_wh_ch32.mat']);
        rawMD292_002 = rawMD292_002.temp_wh_Ch32.values * 0.195; % convert to mV
        ttlMD292_002 = load([g.mainFolder '\MD292_002_kilosort\kilosort25preprocess\TTLsKS.mat']);
        ttlMD292_002 = ttlMD292_002.shocks([5 7 36 41 43 47], :);

        raw_data = rawMD292_002;
        ttl_data = ttlMD292_002;
    end

    % Subplot 1: Raw traces
    ax = nexttile(t, tile_num, [1 1]);
    hold(ax, 'on');
    for ii = 1:size(ttl_data, 1)
        sn_curr = ttl_data(ii, 1) * fs;
        timeaxis = linspace(twin(1), twin(2), (twin(2)-twin(1))*fs) * 1000;
        plot(ax, timeaxis, raw_data(round(sn_curr+twin(1)*fs):round(sn_curr+twin(2)*fs-1)), 'Color', [0 0 0]);
    end
    xlim(ax, [twin(1)*1000, twin(2)*1000]);
    ylabel(ax, 'Voltage (mV)', 'FontSize', g.fontSize2);

    % Set ylim based on example
    if tile_num == 3
        ylim(ax, [-0.12 0.1]);
        yticks(ax, [-0.1 0 0.1]);
    else
        ylim(ax, [-0.5 0.3]);
        yticks(ax, [-0.5 -0.1 0.3]);
    end

    set(ax, 'FontSize', g.fontSize2);
    yLimits = get(ax, 'YLim');

    % Add stimulus indicators
    if tile_num == 3
        % Example 1 (CS): Add speaker symbol (red)
        text(ax, 7, yLimits(2), '🔊', 'Color', 'r', 'FontSize', 16, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    else
        % Example 2 (US): Add gray shaded area and lightning bolt (red, bold)
        fill(ax, [0 10 10 0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        text(ax, 5, yLimits(2), '⚡', 'Color', 'r', 'FontSize', 20, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end

    % Add dashed line at time 0
    plot(ax, [0 0], yLimits, '--k', 'LineWidth', 1);

    ax.Layer = 'top';
    ax.Box = 'off';
    xlabel(ax, 'Time (ms)', 'FontSize', g.fontSize2);
    title(ax, example_name, 'FontSize', g.fontSize1);
end

%% Row 2 Functions
function ax = panel_example_ISI(t, g, tile_num, cell_type)
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
    text(acg_ax, 180, y - (ylim_acg(2)-ylim_acg(1))*0.08, '20 ms', 'HorizontalAlignment', 'center', 'FontSize', g.fontSize2);
    acg_ax.XColor = 'none';
    acg_ax.YColor = 'none';

    % Waveform inset (top right, below ACG)
    wf_ax = axes('Position', [p(1)+p(3)*0.65 p(2)+p(4)*0.20 p(3)*0.3 p(4)*0.3]);
    plot(wf_ax, g.cell_metrics.waveforms.filt{cellID}, 'Color', color_primary, 'LineWidth', 2);
    hold(wf_ax, 'on');
    y = min(g.cell_metrics.waveforms.filt{cellID}) * 1.2;
    plot(wf_ax, [33 48], [y y], 'k', 'LineWidth', 2);
    text(wf_ax, 40, y*1.2, '0.5 ms', 'HorizontalAlignment', 'center', 'FontSize', g.fontSize2);
    wf_ax.XColor = 'none';
    wf_ax.YColor = 'none';
end

function ax = panel_spike_features(t, g, tile_num)
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

function ax = panel_normalized_waveforms(t, g, tile_num)
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
    text(ax, 0.2, -2, ['PN (n=' num2str(sum(idx_PN)) ')'], 'Color', g.colors.PN_primary, 'FontSize', g.fontSize2);
    text(ax, 0.2, -2.5, ['IN (n=' num2str(sum(idx_IN)) ')'], 'Color', g.colors.IN_primary, 'FontSize', g.fontSize2);
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
function ax = panel_firing_rate_distribution(t, g, tile_num, brain_region, show_ylabel)
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
        legend(ax, {['PN (n=' num2str(sum(idx_PN)) ')'], ['IN (n=' num2str(sum(idx_IN)) ')']}, 'FontSize', g.fontSize2, 'Location', 'northeast', 'Box', 'off');
    else
        % All neurons together
        firing_rate = g.cell_metrics.firingRate(idx_region);
        bin_edges = logspace(log10(0.1), log10(30), 20);
        histogram(ax, firing_rate, 'BinEdges', bin_edges, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        text(ax, 0.6, 0.9, ['n=' num2str(sum(idx_region))], 'Units', 'normalized', 'FontSize', g.fontSize2);
    end

    xlabel(ax, 'Hz', 'FontSize', g.fontSize2);
    if show_ylabel
        ylabel(ax, 'Count', 'FontSize', g.fontSize2);
    end
    if strcmp(brain_region, 'Astria')
        title(ax, 'AStria FR', 'FontSize', g.fontSize1);
    else
        title(ax, [brain_region ' FR'], 'FontSize', g.fontSize1);
    end
    xlim(ax, [0.1 30]);
    ylim(ax, [0 50]);
    set(ax, 'XScale', 'log');
    xticks(ax, [0.1 1 10]);
    yticks(ax, [0 25 50]);
    set(ax, 'FontSize', g.fontSize2);
    ax.Box = 'off';
end

%% Row 4 Functions
function ax = panel_burst_index_distribution(t, g, tile_num, brain_region, show_ylabel)
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
        legend(ax, {['PN (n=' num2str(sum(idx_PN)) ')'], ['IN (n=' num2str(sum(idx_IN)) ')']}, 'FontSize', g.fontSize2, 'Location', 'northeast', 'Box', 'off');
    else
        % All neurons together
        burst_index = g.cell_metrics.burstIndex_Royer2012(idx_region);
        bin_edges = logspace(log10(0.1), log10(100), 20);
        histogram(ax, burst_index, 'BinEdges', bin_edges, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        text(ax, 0.6, 0.9, ['n=' num2str(sum(idx_region))], 'Units', 'normalized', 'FontSize', g.fontSize2);
    end

    if show_ylabel
        ylabel(ax, 'Count', 'FontSize', g.fontSize2);
    end
    if strcmp(brain_region, 'Astria')
        title(ax, 'AStria Burst Index', 'FontSize', g.fontSize1);
    else
        title(ax, [brain_region ' Burst Index'], 'FontSize', g.fontSize1);
    end
    xlim(ax, [0.1 100]);
    ylim(ax, [0 30]);
    set(ax, 'XScale', 'log');
    xticks(ax, [0.1 1 10]);
    yticks(ax, [0 15 30]);
    set(ax, 'FontSize', g.fontSize2);
    ax.Box = 'off';
end

%% Statistics Export Function
function export_figure1_to_excel(g)
    output_filename = 'figure_1_data.xlsx';

    if exist(output_filename, 'file')
        delete(output_filename);
    end

    % Get indices for all analyses
    idx_PN_all = strcmp(g.cell_metrics.putativeCellType, 'PN');
    idx_IN_all = strcmp(g.cell_metrics.putativeCellType, 'IN');
    idx_unknown_all = strcmp(g.cell_metrics.putativeCellType, 'unknown');
    idx_LABA = strcmp(g.cell_metrics.brainRegion, 'LA') | strcmp(g.cell_metrics.brainRegion, 'BA');

    regions = {'LA', 'BA', 'Astria', 'CeA'};

    %% Sheet 1: Summary - Sample sizes
    sheet1 = {};
    sheet1{1,1} = 'FIGURE 1: SAMPLE SIZES';
    sheet1{2,1} = ['Generated: ' datestr(now)];
    sheet1{3,1} = '';

    unique_animals = unique(g.cell_metrics.animal);
    n_animals = length(unique_animals);
    n_total = length(g.cell_metrics.cluID);

    sheet1{4,1} = 'Number of animals (N)';
    sheet1{4,2} = n_animals;
    sheet1{5,1} = 'Animal IDs';
    sheet1{5,2} = strjoin(unique_animals, ', ');
    sheet1{6,1} = '';

    sheet1{7,1} = 'Total neurons';
    sheet1{7,2} = n_total;
    sheet1{8,1} = 'Putative pyramidal neurons (PNs)';
    sheet1{8,2} = sum(idx_PN_all);
    sheet1{8,3} = sprintf('%.1f%%', 100*sum(idx_PN_all)/n_total);
    sheet1{9,1} = 'Putative interneurons (INs)';
    sheet1{9,2} = sum(idx_IN_all);
    sheet1{9,3} = sprintf('%.1f%%', 100*sum(idx_IN_all)/n_total);
    sheet1{10,1} = 'Unknown cell type';
    sheet1{10,2} = sum(idx_unknown_all);
    sheet1{10,3} = sprintf('%.1f%%', 100*sum(idx_unknown_all)/n_total);
    sheet1{11,1} = '';

    row_idx = 12;
    sheet1{row_idx,1} = 'By brain region:';
    row_idx = row_idx + 1;
    for r = 1:length(regions)
        idx_region = strcmp(g.cell_metrics.brainRegion, regions{r});
        n_region = sum(idx_region);
        sheet1{row_idx,1} = regions{r};
        sheet1{row_idx,2} = n_region;
        sheet1{row_idx,3} = sprintf('%.1f%%', 100*n_region/n_total);
        row_idx = row_idx + 1;

        if strcmp(regions{r}, 'LA') || strcmp(regions{r}, 'BA')
            idx_PN = idx_region & idx_PN_all;
            idx_IN = idx_region & idx_IN_all;
            idx_unk = idx_region & idx_unknown_all;
            sheet1{row_idx,1} = '  PNs';
            sheet1{row_idx,2} = sum(idx_PN);
            sheet1{row_idx,3} = sprintf('%.1f%% of %s', 100*sum(idx_PN)/n_region, regions{r});
            row_idx = row_idx + 1;
            sheet1{row_idx,1} = '  INs';
            sheet1{row_idx,2} = sum(idx_IN);
            sheet1{row_idx,3} = sprintf('%.1f%% of %s', 100*sum(idx_IN)/n_region, regions{r});
            row_idx = row_idx + 1;
            sheet1{row_idx,1} = '  Unknown';
            sheet1{row_idx,2} = sum(idx_unk);
            sheet1{row_idx,3} = sprintf('%.1f%% of %s', 100*sum(idx_unk)/n_region, regions{r});
            row_idx = row_idx + 1;
        end
    end

    writecell(sheet1, output_filename, 'Sheet', 'Summary_SampleSizes');

    %% Sheet 2: Panel D - Spike Features (LA+BA combined and separate)
    idx_PN_LABA = idx_PN_all & idx_LABA;
    idx_IN_LABA = idx_IN_all & idx_LABA;
    idx_unk_LABA = idx_unknown_all & idx_LABA;

    idx_LA = strcmp(g.cell_metrics.brainRegion, 'LA');
    idx_BA = strcmp(g.cell_metrics.brainRegion, 'BA');
    idx_PN_LA = idx_PN_all & idx_LA;
    idx_IN_LA = idx_IN_all & idx_LA;
    idx_PN_BA = idx_PN_all & idx_BA;
    idx_IN_BA = idx_IN_all & idx_BA;

    % Combined LA+BA
    ttp_PN = g.cell_metrics.spikes.ttp(idx_PN_LABA);
    ttp_IN = g.cell_metrics.spikes.ttp(idx_IN_LABA);
    fr_PN = g.cell_metrics.firingRate(idx_PN_LABA);
    fr_IN = g.cell_metrics.firingRate(idx_IN_LABA);

    [p_ttp, ~, stats_ttp] = ranksum(ttp_PN, ttp_IN);
    [p_fr, ~, stats_fr] = ranksum(fr_PN, fr_IN);

    % LA only
    ttp_PN_LA = g.cell_metrics.spikes.ttp(idx_PN_LA);
    ttp_IN_LA = g.cell_metrics.spikes.ttp(idx_IN_LA);
    fr_PN_LA = g.cell_metrics.firingRate(idx_PN_LA);
    fr_IN_LA = g.cell_metrics.firingRate(idx_IN_LA);

    [p_ttp_LA, ~, stats_ttp_LA] = ranksum(ttp_PN_LA, ttp_IN_LA);
    [p_fr_LA, ~, stats_fr_LA] = ranksum(fr_PN_LA, fr_IN_LA);

    % BA only
    ttp_PN_BA = g.cell_metrics.spikes.ttp(idx_PN_BA);
    ttp_IN_BA = g.cell_metrics.spikes.ttp(idx_IN_BA);
    fr_PN_BA = g.cell_metrics.firingRate(idx_PN_BA);
    fr_IN_BA = g.cell_metrics.firingRate(idx_IN_BA);

    [p_ttp_BA, ~, stats_ttp_BA] = ranksum(ttp_PN_BA, ttp_IN_BA);
    [p_fr_BA, ~, stats_fr_BA] = ranksum(fr_PN_BA, fr_IN_BA);

    sheet2 = {};
    sheet2{1,1} = 'PANEL D: SPIKE FEATURES';
    sheet2{2,1} = '';

    % LA+BA Combined
    sheet2{3,1} = '=== LA + BA COMBINED ===';
    sheet2{4,1} = 'Sample sizes';
    sheet2{4,2} = 'n';
    sheet2{5,1} = 'PNs';
    sheet2{5,2} = sum(idx_PN_LABA);
    sheet2{6,1} = 'INs';
    sheet2{6,2} = sum(idx_IN_LABA);
    sheet2{7,1} = 'Unknown';
    sheet2{7,2} = sum(idx_unk_LABA);
    sheet2{8,1} = '';

    sheet2{9,1} = 'Trough-to-peak (TTP) duration (ms)';
    sheet2{9,2} = '';
    sheet2{9,3} = 'Mean';
    sheet2{9,4} = 'SEM';
    sheet2{9,5} = 'Median';
    sheet2{9,6} = 'SD';
    sheet2{10,1} = 'PNs';
    sheet2{10,2} = sum(idx_PN_LABA);
    sheet2{10,3} = mean(ttp_PN);
    sheet2{10,4} = std(ttp_PN)/sqrt(length(ttp_PN));
    sheet2{10,5} = median(ttp_PN);
    sheet2{10,6} = std(ttp_PN);
    sheet2{11,1} = 'INs';
    sheet2{11,2} = sum(idx_IN_LABA);
    sheet2{11,3} = mean(ttp_IN);
    sheet2{11,4} = std(ttp_IN)/sqrt(length(ttp_IN));
    sheet2{11,5} = median(ttp_IN);
    sheet2{11,6} = std(ttp_IN);
    sheet2{12,1} = '';
    sheet2{13,1} = 'Wilcoxon rank-sum test (PN vs IN TTP)';
    sheet2{14,1} = 'Z-statistic';
    sheet2{14,2} = stats_ttp.zval;
    sheet2{15,1} = 'p-value';
    sheet2{15,2} = p_ttp;
    sheet2{15,3} = format_significance(p_ttp);
    sheet2{16,1} = '';

    sheet2{17,1} = 'Firing rate (Hz)';
    sheet2{17,2} = '';
    sheet2{17,3} = 'Mean';
    sheet2{17,4} = 'SEM';
    sheet2{17,5} = 'Median';
    sheet2{17,6} = 'SD';
    sheet2{18,1} = 'PNs';
    sheet2{18,2} = sum(idx_PN_LABA);
    sheet2{18,3} = mean(fr_PN);
    sheet2{18,4} = std(fr_PN)/sqrt(length(fr_PN));
    sheet2{18,5} = median(fr_PN);
    sheet2{18,6} = std(fr_PN);
    sheet2{19,1} = 'INs';
    sheet2{19,2} = sum(idx_IN_LABA);
    sheet2{19,3} = mean(fr_IN);
    sheet2{19,4} = std(fr_IN)/sqrt(length(fr_IN));
    sheet2{19,5} = median(fr_IN);
    sheet2{19,6} = std(fr_IN);
    sheet2{20,1} = '';
    sheet2{21,1} = 'Wilcoxon rank-sum test (PN vs IN firing rate)';
    sheet2{22,1} = 'Z-statistic';
    sheet2{22,2} = stats_fr.zval;
    sheet2{23,1} = 'p-value';
    sheet2{23,2} = p_fr;
    sheet2{23,3} = format_significance(p_fr);
    sheet2{24,1} = '';

    % LA only
    sheet2{25,1} = '=== LA ONLY ===';
    sheet2{26,1} = 'Sample sizes';
    sheet2{26,2} = 'n';
    sheet2{27,1} = 'PNs';
    sheet2{27,2} = sum(idx_PN_LA);
    sheet2{28,1} = 'INs';
    sheet2{28,2} = sum(idx_IN_LA);
    sheet2{29,1} = '';

    sheet2{30,1} = 'Trough-to-peak (TTP) duration (ms)';
    sheet2{30,2} = '';
    sheet2{30,3} = 'Mean';
    sheet2{30,4} = 'SEM';
    sheet2{30,5} = 'Median';
    sheet2{30,6} = 'SD';
    sheet2{31,1} = 'PNs';
    sheet2{31,2} = sum(idx_PN_LA);
    sheet2{31,3} = mean(ttp_PN_LA);
    sheet2{31,4} = std(ttp_PN_LA)/sqrt(length(ttp_PN_LA));
    sheet2{31,5} = median(ttp_PN_LA);
    sheet2{31,6} = std(ttp_PN_LA);
    sheet2{32,1} = 'INs';
    sheet2{32,2} = sum(idx_IN_LA);
    sheet2{32,3} = mean(ttp_IN_LA);
    sheet2{32,4} = std(ttp_IN_LA)/sqrt(length(ttp_IN_LA));
    sheet2{32,5} = median(ttp_IN_LA);
    sheet2{32,6} = std(ttp_IN_LA);
    sheet2{33,1} = '';
    sheet2{34,1} = 'Wilcoxon rank-sum test (PN vs IN TTP)';
    sheet2{35,1} = 'Z-statistic';
    sheet2{35,2} = stats_ttp_LA.zval;
    sheet2{36,1} = 'p-value';
    sheet2{36,2} = p_ttp_LA;
    sheet2{36,3} = format_significance(p_ttp_LA);
    sheet2{37,1} = '';

    sheet2{38,1} = 'Firing rate (Hz)';
    sheet2{38,2} = '';
    sheet2{38,3} = 'Mean';
    sheet2{38,4} = 'SEM';
    sheet2{38,5} = 'Median';
    sheet2{38,6} = 'SD';
    sheet2{39,1} = 'PNs';
    sheet2{39,2} = sum(idx_PN_LA);
    sheet2{39,3} = mean(fr_PN_LA);
    sheet2{39,4} = std(fr_PN_LA)/sqrt(length(fr_PN_LA));
    sheet2{39,5} = median(fr_PN_LA);
    sheet2{39,6} = std(fr_PN_LA);
    sheet2{40,1} = 'INs';
    sheet2{40,2} = sum(idx_IN_LA);
    sheet2{40,3} = mean(fr_IN_LA);
    sheet2{40,4} = std(fr_IN_LA)/sqrt(length(fr_IN_LA));
    sheet2{40,5} = median(fr_IN_LA);
    sheet2{40,6} = std(fr_IN_LA);
    sheet2{41,1} = '';
    sheet2{42,1} = 'Wilcoxon rank-sum test (PN vs IN firing rate)';
    sheet2{43,1} = 'Z-statistic';
    sheet2{43,2} = stats_fr_LA.zval;
    sheet2{44,1} = 'p-value';
    sheet2{44,2} = p_fr_LA;
    sheet2{44,3} = format_significance(p_fr_LA);
    sheet2{45,1} = '';

    % BA only
    sheet2{46,1} = '=== BA ONLY ===';
    sheet2{47,1} = 'Sample sizes';
    sheet2{47,2} = 'n';
    sheet2{48,1} = 'PNs';
    sheet2{48,2} = sum(idx_PN_BA);
    sheet2{49,1} = 'INs';
    sheet2{49,2} = sum(idx_IN_BA);
    sheet2{50,1} = '';

    sheet2{51,1} = 'Trough-to-peak (TTP) duration (ms)';
    sheet2{51,2} = '';
    sheet2{51,3} = 'Mean';
    sheet2{51,4} = 'SEM';
    sheet2{51,5} = 'Median';
    sheet2{51,6} = 'SD';
    sheet2{52,1} = 'PNs';
    sheet2{52,2} = sum(idx_PN_BA);
    sheet2{52,3} = mean(ttp_PN_BA);
    sheet2{52,4} = std(ttp_PN_BA)/sqrt(length(ttp_PN_BA));
    sheet2{52,5} = median(ttp_PN_BA);
    sheet2{52,6} = std(ttp_PN_BA);
    sheet2{53,1} = 'INs';
    sheet2{53,2} = sum(idx_IN_BA);
    sheet2{53,3} = mean(ttp_IN_BA);
    sheet2{53,4} = std(ttp_IN_BA)/sqrt(length(ttp_IN_BA));
    sheet2{53,5} = median(ttp_IN_BA);
    sheet2{53,6} = std(ttp_IN_BA);
    sheet2{54,1} = '';
    sheet2{55,1} = 'Wilcoxon rank-sum test (PN vs IN TTP)';
    sheet2{56,1} = 'Z-statistic';
    sheet2{56,2} = stats_ttp_BA.zval;
    sheet2{57,1} = 'p-value';
    sheet2{57,2} = p_ttp_BA;
    sheet2{57,3} = format_significance(p_ttp_BA);
    sheet2{58,1} = '';

    sheet2{59,1} = 'Firing rate (Hz)';
    sheet2{59,2} = '';
    sheet2{59,3} = 'Mean';
    sheet2{59,4} = 'SEM';
    sheet2{59,5} = 'Median';
    sheet2{59,6} = 'SD';
    sheet2{60,1} = 'PNs';
    sheet2{60,2} = sum(idx_PN_BA);
    sheet2{60,3} = mean(fr_PN_BA);
    sheet2{60,4} = std(fr_PN_BA)/sqrt(length(fr_PN_BA));
    sheet2{60,5} = median(fr_PN_BA);
    sheet2{60,6} = std(fr_PN_BA);
    sheet2{61,1} = 'INs';
    sheet2{61,2} = sum(idx_IN_BA);
    sheet2{61,3} = mean(fr_IN_BA);
    sheet2{61,4} = std(fr_IN_BA)/sqrt(length(fr_IN_BA));
    sheet2{61,5} = median(fr_IN_BA);
    sheet2{61,6} = std(fr_IN_BA);
    sheet2{62,1} = '';
    sheet2{63,1} = 'Wilcoxon rank-sum test (PN vs IN firing rate)';
    sheet2{64,1} = 'Z-statistic';
    sheet2{64,2} = stats_fr_BA.zval;
    sheet2{65,1} = 'p-value';
    sheet2{65,2} = p_fr_BA;
    sheet2{65,3} = format_significance(p_fr_BA);

    writecell(sheet2, output_filename, 'Sheet', 'PanelD_SpikeFeatures');

    %% Sheet 3: Panel E - Firing Rate by Region
    sheet3 = {};
    sheet3{1,1} = 'PANEL E: FIRING RATE DISTRIBUTIONS BY REGION';
    sheet3{2,1} = '';

    row_idx = 3;
    for r = 1:length(regions)
        sheet3{row_idx,1} = ['--- ' regions{r} ' ---'];
        row_idx = row_idx + 1;
        idx_region = strcmp(g.cell_metrics.brainRegion, regions{r});

        if strcmp(regions{r}, 'LA') || strcmp(regions{r}, 'BA')
            idx_PN = idx_region & idx_PN_all;
            idx_IN = idx_region & idx_IN_all;

            fr_PN_region = g.cell_metrics.firingRate(idx_PN);
            fr_IN_region = g.cell_metrics.firingRate(idx_IN);

            sheet3{row_idx,1} = 'Cell Type';
            sheet3{row_idx,2} = 'n';
            sheet3{row_idx,3} = 'Mean (Hz)';
            sheet3{row_idx,4} = 'SEM';
            sheet3{row_idx,5} = 'Median';
            sheet3{row_idx,6} = 'SD';
            row_idx = row_idx + 1;

            sheet3{row_idx,1} = 'PNs';
            sheet3{row_idx,2} = sum(idx_PN);
            sheet3{row_idx,3} = mean(fr_PN_region);
            sheet3{row_idx,4} = std(fr_PN_region)/sqrt(length(fr_PN_region));
            sheet3{row_idx,5} = median(fr_PN_region);
            sheet3{row_idx,6} = std(fr_PN_region);
            row_idx = row_idx + 1;

            sheet3{row_idx,1} = 'INs';
            sheet3{row_idx,2} = sum(idx_IN);
            sheet3{row_idx,3} = mean(fr_IN_region);
            sheet3{row_idx,4} = std(fr_IN_region)/sqrt(length(fr_IN_region));
            sheet3{row_idx,5} = median(fr_IN_region);
            sheet3{row_idx,6} = std(fr_IN_region);
            row_idx = row_idx + 1;

            [p_region, ~, stats_region] = ranksum(fr_PN_region, fr_IN_region);
            sheet3{row_idx,1} = 'Wilcoxon rank-sum (PN vs IN)';
            row_idx = row_idx + 1;
            sheet3{row_idx,1} = 'Z-statistic';
            sheet3{row_idx,2} = stats_region.zval;
            row_idx = row_idx + 1;
            sheet3{row_idx,1} = 'p-value';
            sheet3{row_idx,2} = p_region;
            sheet3{row_idx,3} = format_significance(p_region);
            row_idx = row_idx + 1;
        else
            fr_region = g.cell_metrics.firingRate(idx_region);

            sheet3{row_idx,1} = 'All neurons';
            sheet3{row_idx,2} = 'n';
            sheet3{row_idx,3} = 'Mean (Hz)';
            sheet3{row_idx,4} = 'SEM';
            sheet3{row_idx,5} = 'Median';
            sheet3{row_idx,6} = 'SD';
            row_idx = row_idx + 1;

            sheet3{row_idx,1} = regions{r};
            sheet3{row_idx,2} = sum(idx_region);
            sheet3{row_idx,3} = mean(fr_region);
            sheet3{row_idx,4} = std(fr_region)/sqrt(length(fr_region));
            sheet3{row_idx,5} = median(fr_region);
            sheet3{row_idx,6} = std(fr_region);
            row_idx = row_idx + 1;
        end

        sheet3{row_idx,1} = '';
        row_idx = row_idx + 1;
    end

    writecell(sheet3, output_filename, 'Sheet', 'PanelE_FiringRate');

    %% Sheet 4: Panel F - Burst Index by Region
    sheet4 = {};
    sheet4{1,1} = 'PANEL F: BURST INDEX DISTRIBUTIONS BY REGION';
    sheet4{2,1} = '';

    row_idx = 3;
    for r = 1:length(regions)
        sheet4{row_idx,1} = ['--- ' regions{r} ' ---'];
        row_idx = row_idx + 1;
        idx_region = strcmp(g.cell_metrics.brainRegion, regions{r});

        if strcmp(regions{r}, 'LA') || strcmp(regions{r}, 'BA')
            idx_PN = idx_region & idx_PN_all;
            idx_IN = idx_region & idx_IN_all;

            bi_PN_region = g.cell_metrics.burstIndex_Royer2012(idx_PN);
            bi_IN_region = g.cell_metrics.burstIndex_Royer2012(idx_IN);

            sheet4{row_idx,1} = 'Cell Type';
            sheet4{row_idx,2} = 'n';
            sheet4{row_idx,3} = 'Mean';
            sheet4{row_idx,4} = 'SEM';
            sheet4{row_idx,5} = 'Median';
            sheet4{row_idx,6} = 'SD';
            row_idx = row_idx + 1;

            sheet4{row_idx,1} = 'PNs';
            sheet4{row_idx,2} = sum(idx_PN);
            sheet4{row_idx,3} = mean(bi_PN_region, 'omitnan');
            sheet4{row_idx,4} = std(bi_PN_region, 'omitnan')/sqrt(sum(~isnan(bi_PN_region)));
            sheet4{row_idx,5} = median(bi_PN_region, 'omitnan');
            sheet4{row_idx,6} = std(bi_PN_region, 'omitnan');
            row_idx = row_idx + 1;

            sheet4{row_idx,1} = 'INs';
            sheet4{row_idx,2} = sum(idx_IN);
            sheet4{row_idx,3} = mean(bi_IN_region, 'omitnan');
            sheet4{row_idx,4} = std(bi_IN_region, 'omitnan')/sqrt(sum(~isnan(bi_IN_region)));
            sheet4{row_idx,5} = median(bi_IN_region, 'omitnan');
            sheet4{row_idx,6} = std(bi_IN_region, 'omitnan');
            row_idx = row_idx + 1;

            [p_region, ~, stats_region] = ranksum(bi_PN_region(~isnan(bi_PN_region)), bi_IN_region(~isnan(bi_IN_region)));
            sheet4{row_idx,1} = 'Wilcoxon rank-sum (PN vs IN)';
            row_idx = row_idx + 1;
            sheet4{row_idx,1} = 'Z-statistic';
            sheet4{row_idx,2} = stats_region.zval;
            row_idx = row_idx + 1;
            sheet4{row_idx,1} = 'p-value';
            sheet4{row_idx,2} = p_region;
            sheet4{row_idx,3} = format_significance(p_region);
            row_idx = row_idx + 1;
        else
            bi_region = g.cell_metrics.burstIndex_Royer2012(idx_region);

            sheet4{row_idx,1} = 'All neurons';
            sheet4{row_idx,2} = 'n';
            sheet4{row_idx,3} = 'Mean';
            sheet4{row_idx,4} = 'SEM';
            sheet4{row_idx,5} = 'Median';
            sheet4{row_idx,6} = 'SD';
            row_idx = row_idx + 1;

            sheet4{row_idx,1} = regions{r};
            sheet4{row_idx,2} = sum(idx_region);
            sheet4{row_idx,3} = mean(bi_region, 'omitnan');
            sheet4{row_idx,4} = std(bi_region, 'omitnan')/sqrt(sum(~isnan(bi_region)));
            sheet4{row_idx,5} = median(bi_region, 'omitnan');
            sheet4{row_idx,6} = std(bi_region, 'omitnan');
            row_idx = row_idx + 1;
        end

        sheet4{row_idx,1} = '';
        row_idx = row_idx + 1;
    end

    writecell(sheet4, output_filename, 'Sheet', 'PanelF_BurstIndex');

    %% Sheet 5: Cross-Region Comparisons
    idx_LA_PN = strcmp(g.cell_metrics.brainRegion, 'LA') & idx_PN_all;
    idx_LA_IN = strcmp(g.cell_metrics.brainRegion, 'LA') & idx_IN_all;
    idx_BA_PN = strcmp(g.cell_metrics.brainRegion, 'BA') & idx_PN_all;
    idx_BA_IN = strcmp(g.cell_metrics.brainRegion, 'BA') & idx_IN_all;
    idx_Astria = strcmp(g.cell_metrics.brainRegion, 'Astria');
    idx_CeA = strcmp(g.cell_metrics.brainRegion, 'CeA');

    fr_LA_PN = g.cell_metrics.firingRate(idx_LA_PN);
    fr_LA_IN = g.cell_metrics.firingRate(idx_LA_IN);
    fr_BA_PN = g.cell_metrics.firingRate(idx_BA_PN);
    fr_BA_IN = g.cell_metrics.firingRate(idx_BA_IN);
    fr_Astria = g.cell_metrics.firingRate(idx_Astria);
    fr_CeA = g.cell_metrics.firingRate(idx_CeA);

    % Ensure column vectors
    fr_LA_PN = fr_LA_PN(:);
    fr_LA_IN = fr_LA_IN(:);
    fr_BA_PN = fr_BA_PN(:);
    fr_BA_IN = fr_BA_IN(:);
    fr_Astria = fr_Astria(:);
    fr_CeA = fr_CeA(:);

    group_labels = [ones(sum(idx_LA_PN), 1); 2*ones(sum(idx_LA_IN), 1); ...
        3*ones(sum(idx_BA_PN), 1); 4*ones(sum(idx_BA_IN), 1); ...
        5*ones(sum(idx_Astria), 1); 6*ones(sum(idx_CeA), 1)];
    fr_combined = [fr_LA_PN; fr_LA_IN; fr_BA_PN; fr_BA_IN; fr_Astria; fr_CeA];

    [p_kw, tbl_kw, stats_kw] = kruskalwallis(fr_combined, group_labels, 'off');

    sheet5 = {};
    sheet5{1,1} = 'CROSS-REGION COMPARISONS: FIRING RATE';
    sheet5{2,1} = '';
    sheet5{3,1} = 'Sample sizes';
    sheet5{3,2} = 'n';
    sheet5{4,1} = 'LA PNs';
    sheet5{4,2} = sum(idx_LA_PN);
    sheet5{5,1} = 'LA INs';
    sheet5{5,2} = sum(idx_LA_IN);
    sheet5{6,1} = 'BA PNs';
    sheet5{6,2} = sum(idx_BA_PN);
    sheet5{7,1} = 'BA INs';
    sheet5{7,2} = sum(idx_BA_IN);
    sheet5{8,1} = 'AStria';
    sheet5{8,2} = sum(idx_Astria);
    sheet5{9,1} = 'CeA';
    sheet5{9,2} = sum(idx_CeA);
    sheet5{10,1} = '';

    sheet5{11,1} = 'Kruskal-Wallis Test: Firing rate across regions';
    sheet5{12,1} = 'Chi-square statistic';
    sheet5{12,2} = tbl_kw{2,5};
    sheet5{13,1} = 'Degrees of freedom';
    sheet5{13,2} = tbl_kw{2,3};
    sheet5{14,1} = 'p-value';
    sheet5{14,2} = p_kw;
    sheet5{14,3} = format_significance(p_kw);
    sheet5{15,1} = '';

    row_idx = 16;
    if p_kw < 0.05
        sheet5{row_idx,1} = 'Post-hoc pairwise comparisons (Dunn-Sidak)';
        row_idx = row_idx + 1;
        sheet5{row_idx,1} = 'Comparison';
        sheet5{row_idx,2} = 'p-value';
        sheet5{row_idx,3} = 'Significance';
        row_idx = row_idx + 1;

        c = multcompare(stats_kw, 'Display', 'off', 'CType', 'dunn-sidak');
        region_names = {'LA_PN', 'LA_IN', 'BA_PN', 'BA_IN', 'AStria', 'CeA'};
        for i = 1:size(c, 1)
            sheet5{row_idx,1} = [region_names{c(i,1)} ' vs ' region_names{c(i,2)}];
            sheet5{row_idx,2} = c(i,6);
            sheet5{row_idx,3} = format_significance(c(i,6));
            row_idx = row_idx + 1;
        end
    else
        sheet5{row_idx,1} = 'No post-hoc tests (K-W p >= 0.05)';
    end

    writecell(sheet5, output_filename, 'Sheet', 'CrossRegion_FiringRate');

    %% Sheet 6: Cross-Region Comparisons - Burst Index
    bi_LA_PN = g.cell_metrics.burstIndex_Royer2012(idx_LA_PN);
    bi_LA_IN = g.cell_metrics.burstIndex_Royer2012(idx_LA_IN);
    bi_BA_PN = g.cell_metrics.burstIndex_Royer2012(idx_BA_PN);
    bi_BA_IN = g.cell_metrics.burstIndex_Royer2012(idx_BA_IN);
    bi_Astria = g.cell_metrics.burstIndex_Royer2012(idx_Astria);
    bi_CeA = g.cell_metrics.burstIndex_Royer2012(idx_CeA);

    % Remove NaNs before analysis
    bi_LA_PN_clean = bi_LA_PN(~isnan(bi_LA_PN));
    bi_LA_IN_clean = bi_LA_IN(~isnan(bi_LA_IN));
    bi_BA_PN_clean = bi_BA_PN(~isnan(bi_BA_PN));
    bi_BA_IN_clean = bi_BA_IN(~isnan(bi_BA_IN));
    bi_Astria_clean = bi_Astria(~isnan(bi_Astria));
    bi_CeA_clean = bi_CeA(~isnan(bi_CeA));

    % Ensure column vectors
    bi_LA_PN_clean = bi_LA_PN_clean(:);
    bi_LA_IN_clean = bi_LA_IN_clean(:);
    bi_BA_PN_clean = bi_BA_PN_clean(:);
    bi_BA_IN_clean = bi_BA_IN_clean(:);
    bi_Astria_clean = bi_Astria_clean(:);
    bi_CeA_clean = bi_CeA_clean(:);

    group_labels_bi = [ones(length(bi_LA_PN_clean), 1); 2*ones(length(bi_LA_IN_clean), 1); ...
        3*ones(length(bi_BA_PN_clean), 1); 4*ones(length(bi_BA_IN_clean), 1); ...
        5*ones(length(bi_Astria_clean), 1); 6*ones(length(bi_CeA_clean), 1)];
    bi_combined = [bi_LA_PN_clean; bi_LA_IN_clean; bi_BA_PN_clean; bi_BA_IN_clean; bi_Astria_clean; bi_CeA_clean];

    [p_kw_bi, tbl_kw_bi, stats_kw_bi] = kruskalwallis(bi_combined, group_labels_bi, 'off');

    sheet6 = {};
    sheet6{1,1} = 'CROSS-REGION COMPARISONS: BURST INDEX';
    sheet6{2,1} = '';
    sheet6{3,1} = 'Sample sizes (excluding NaN)';
    sheet6{3,2} = 'n';
    sheet6{4,1} = 'LA PNs';
    sheet6{4,2} = length(bi_LA_PN_clean);
    sheet6{5,1} = 'LA INs';
    sheet6{5,2} = length(bi_LA_IN_clean);
    sheet6{6,1} = 'BA PNs';
    sheet6{6,2} = length(bi_BA_PN_clean);
    sheet6{7,1} = 'BA INs';
    sheet6{7,2} = length(bi_BA_IN_clean);
    sheet6{8,1} = 'AStria';
    sheet6{8,2} = length(bi_Astria_clean);
    sheet6{9,1} = 'CeA';
    sheet6{9,2} = length(bi_CeA_clean);
    sheet6{10,1} = '';

    sheet6{11,1} = 'Kruskal-Wallis Test: Burst index across regions';
    sheet6{12,1} = 'Chi-square statistic';
    sheet6{12,2} = tbl_kw_bi{2,5};
    sheet6{13,1} = 'Degrees of freedom';
    sheet6{13,2} = tbl_kw_bi{2,3};
    sheet6{14,1} = 'p-value';
    sheet6{14,2} = p_kw_bi;
    sheet6{14,3} = format_significance(p_kw_bi);
    sheet6{15,1} = '';

    row_idx = 16;
    if p_kw_bi < 0.05
        sheet6{row_idx,1} = 'Post-hoc pairwise comparisons (Dunn-Sidak)';
        row_idx = row_idx + 1;
        sheet6{row_idx,1} = 'Comparison';
        sheet6{row_idx,2} = 'p-value';
        sheet6{row_idx,3} = 'Significance';
        row_idx = row_idx + 1;

        c_bi = multcompare(stats_kw_bi, 'Display', 'off', 'CType', 'dunn-sidak');
        region_names = {'LA_PN', 'LA_IN', 'BA_PN', 'BA_IN', 'AStria', 'CeA'};
        for i = 1:size(c_bi, 1)
            sheet6{row_idx,1} = [region_names{c_bi(i,1)} ' vs ' region_names{c_bi(i,2)}];
            sheet6{row_idx,2} = c_bi(i,6);
            sheet6{row_idx,3} = format_significance(c_bi(i,6));
            row_idx = row_idx + 1;
        end
    else
        sheet6{row_idx,1} = 'No post-hoc tests (K-W p >= 0.05)';
    end

    writecell(sheet6, output_filename, 'Sheet', 'CrossRegion_BurstIndex');
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
