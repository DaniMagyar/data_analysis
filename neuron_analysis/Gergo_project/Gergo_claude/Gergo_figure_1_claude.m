%% Gergo Figure 1 - Individual Panel Export
% Creates and saves each panel as a separate PNG file
% All fonts are Arial 10pt for publication
% No panel labels (A, B, C, etc.)

g.mainFolder = 'C:\Users\dmagyar\Desktop\Gergo';
g.subfolders.PFC = 'C:\Users\dmagyar\Desktop\Gergo\Opto 1 Gergo mPFCbe vetito\Opto 1';
g.subfolders.PFC_tagging  = 'C:\Users\dmagyar\Desktop\Gergo\NagyGergo_Opto1_PrL';
folderList_PFC = dir(g.subfolders.PFC);
folderList_PFC = folderList_PFC([folderList_PFC.isdir] & ~ismember({folderList_PFC.name}, {'.', '..'}));
g.fullPaths_PFC = fullfile(g.subfolders.PFC, {folderList_PFC.name});

g.subfolders.DMS = 'C:\Users\dmagyar\Desktop\Gergo\Opto 2 Gergo DMSbe vetito\Opto 2';
g.subfolders.DMS_tagging = 'C:\Users\dmagyar\Desktop\Gergo\NagyGergo_Opto2_DMS';
folderList_DMS = dir(g.subfolders.DMS);
folderList_DMS = folderList_DMS([folderList_DMS.isdir] & ~ismember({folderList_DMS.name}, {'.', '..'}));
g.fullPaths_DMS = fullfile(g.subfolders.DMS, {folderList_DMS.name});

g.pre_time = 1;
g.post_time = 1;
g.bin_time = 0.001;
g.smoothvalue = 201;
g.timeaxis = -g.pre_time:g.bin_time:g.post_time;
g.colors = BAfc_colors;
g.test_time = 0.6;
g.exctestvalue = 2;
g.inhtestvalue = 2;
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
% Font sizes set to 10pt for 4x3.5cm panels
g.fontSize2 = 10;  % Axis labels and ticks
g.fontSize1 = 10;  % Title font
g.xlinewidth = 1.5;  % Stimulus line width
g.axisLineWidth = 1;  % Axis line width
g.markerSize = 3;  % Raster plot marker size
g.fontName = 'Arial';  % Font type
g.optopre = 0.02;
g.optopost = 0.02;
g.optobin = 0.001;
g.optotimeaxis = -g.optopre:g.optobin:g.optopost;

% Example neuron cell IDs
g.example_cellID_PFC = 4;
g.example_cellID_DMS = 27;
g.example_cellID_PFC_opto = 4;
g.example_cellID_DMS_opto = 16;

% Optotagging TTL limit
g.opto_ttl_limit = [500 1000];

% Output folder
g.outputFolder = 'C:\Users\dmagyar\Documents\data_analysis\neuron_analysis\Gergo_project\Gergo_claude';
if ~exist(g.outputFolder, 'dir')
    mkdir(g.outputFolder);
end

% Panel dimensions - sized for 4 panels to fit side-by-side on a Word page
% Standard page width ~17cm, so each panel ~4.25cm width to fit 4 across
panel_width = 5;   % cm
panel_height = 4.4;  % cm (maintaining reasonable aspect ratio)

% Standard axes position for all panels (normalized units: [left bottom width height])
% This ensures all plot areas are the same size regardless of labels/colorbars
% Need substantial space on left and bottom for labels to be fully visible
axes_pos_standard = [0.30 0.25 0.60 0.65];  % Standard panels
axes_pos_heatmap = [0.30 0.25 0.45 0.65];   % All heatmaps (need room for colorbar)
axes_pos_heatmap2 = [0.30 0.25 0.60 0.65];   % All heatmaps (need room for colorbar)

%% Create cell_metrics
if ~isfield(g, 'cell_metrics')
    g.cell_metrics.spikes.times = {};
    g.cell_metrics.projection = {};
    g.cell_metrics.general.shockTTL = {};
    g.cell_metrics.general.optoTTL = {};
    g.cell_metrics.cellID = {};
    % Load PFC projection neurons
    for ii = 1:size(g.fullPaths_PFC,2)
        cd(g.fullPaths_PFC{ii})
        cellIDs = dir;
        cellIDs = cellIDs(~[cellIDs.isdir]);
        cellIDs = cellIDs(contains({cellIDs.name},'GR'));
        for jj = 1:size(cellIDs,1)
            spikes = load(cellIDs(jj).name, 'TS');
            ttl = load('shocKTTL.mat', 'shockTTL');
            optofilename = dir('*events*');
            [~, timestamps, ~] = load_open_ephys_data(optofilename.name);
            timestamps = timestamps(1:2:end);
            opto_ttl = setdiff(timestamps, ttl.shockTTL);
            g.cell_metrics.general.optoTTL{end+1} = opto_ttl;
            g.cell_metrics.spikes.times{end+1} = spikes.TS/10000;
            g.cell_metrics.projection{end+1} = 'PFC';
            g.cell_metrics.general.shockTTL{end+1} = ttl.shockTTL;
            g.cell_metrics.cellID{end+1} = size(g.cell_metrics.cellID,2)+1;
        end
    end
    % Load DMS projection neurons
    for ii = 1:size(g.fullPaths_DMS,2)
        cd(g.fullPaths_DMS{ii})
        cellIDs = dir;
        cellIDs = cellIDs(~[cellIDs.isdir]);
        cellIDs = cellIDs(contains({cellIDs.name},'GR'));
        for jj = 1:size(cellIDs,1)
            spikes = load(cellIDs(jj).name, 'TS');
            ttl = load('shocKTTL.mat', 'shockTTL');
            optofilename = dir('*events*');
            [~, timestamps, ~] = load_open_ephys_data(optofilename.name);
            timestamps = timestamps(1:2:end);
            opto_ttl = setdiff(timestamps, ttl.shockTTL);
            g.cell_metrics.general.optoTTL{end+1} = opto_ttl;
            g.cell_metrics.spikes.times{end+1} = spikes.TS/10000;
            g.cell_metrics.projection{end+1} = 'DMS';
            g.cell_metrics.general.shockTTL{end+1} = ttl.shockTTL;
            g.cell_metrics.cellID{end+1} = size(g.cell_metrics.cellID,2)+1;
        end
    end
end
clearvars -except g panel_width panel_height axes_pos_standard axes_pos_heatmap axes_pos_heatmap2

%% PSTH
psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'shockTTL', 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

% Z-score based on baseline period only
baseline_bins = 1:round(g.pre_time/g.bin_time);
baseline_mean = mean(psth_spx(:, baseline_bins), 2);
baseline_std = std(psth_spx(:, baseline_bins), 0, 2);
psth_spx_zscore = (psth_spx - baseline_mean) ./ baseline_std;
psth_spx_zscore = smoothdata(psth_spx_zscore, 2, 'sgolay', g.smoothvalue);
idx_PFC = strcmp(g.cell_metrics.projection, 'PFC')';
idx_DMS = strcmp(g.cell_metrics.projection, 'DMS')';
psth_spx_zscore_PFC = psth_spx_zscore(idx_PFC,:);
psth_spx_zscore_DMS = psth_spx_zscore(idx_DMS,:);
psth_spx_PFC = psth_spx(idx_PFC,:);
psth_spx_DMS = psth_spx(idx_DMS,:);

% Calculate response onset latency for each neuron
onset_latency = nan(size(psth_spx_zscore,1),1);
for ii = 1:size(psth_spx_zscore,1)
    datasegment = psth_spx_zscore(ii,g.roi);
    if max(datasegment) >= g.exctestvalue
        onset_latency(ii) = find(datasegment >= g.exctestvalue, 1, 'first') * g.bin_time * 1000;
    elseif min(datasegment) <= -g.inhtestvalue
        onset_latency(ii) = find(datasegment <= -g.inhtestvalue, 1, 'first') * g.bin_time * 1000;
    else
        onset_latency(ii) = 999999;
    end
end

% Sort by onset latency
[~, order_PFC] = sort(onset_latency(idx_PFC), 'ascend');
[~, order_DMS] = sort(onset_latency(idx_DMS), 'ascend');

psth_spx_zscore_PFC = psth_spx_zscore_PFC(order_PFC,:);
psth_spx_zscore_DMS = psth_spx_zscore_DMS(order_DMS,:);
psth_spx_PFC = psth_spx_PFC(order_PFC,:);
psth_spx_DMS = psth_spx_DMS(order_DMS,:);

% Set color limits
g.clim = [prctile(psth_spx_zscore(:), 0.5) prctile(psth_spx_zscore(:), 99.5)];

% Calculate responsive neurons
n_responsive_PFC = sum(onset_latency(idx_PFC) < 999999);
n_total_PFC = sum(idx_PFC);
pct_responsive_PFC = (n_responsive_PFC / n_total_PFC) * 100;
pct_nonresponsive_PFC = 100 - pct_responsive_PFC;

n_responsive_DMS = sum(onset_latency(idx_DMS) < 999999);
n_total_DMS = sum(idx_DMS);
pct_responsive_DMS = (n_responsive_DMS / n_total_DMS) * 100;
pct_nonresponsive_DMS = 100 - pct_responsive_DMS;

fprintf('\n========================================\n');
fprintf('Creating and saving individual panels...\n');
fprintf('========================================\n\n');

%% Panel 1: PFC example raster
fprintf('Panel 1: PFC example raster...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
spike_times_PFC = g.cell_metrics.spikes.times{g.example_cellID_PFC};
stimulus_times_PFC = g.cell_metrics.general.shockTTL{g.example_cellID_PFC};
stimulus_times_PFC = stimulus_times_PFC(1:min(80, length(stimulus_times_PFC)));
hold on;
for trial = 1:length(stimulus_times_PFC)
    aligned_spikes = spike_times_PFC - stimulus_times_PFC(trial);
    valid_spikes = aligned_spikes(aligned_spikes >= -g.pre_time & aligned_spikes <= g.post_time);
    scatter(valid_spikes, trial * ones(size(valid_spikes)), g.markerSize, 'k', 'filled');
end
xlim([-g.pre_time g.post_time]);
xticks([-g.pre_time 0 g.post_time]);
% ylabel removed for cleaner panels
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Set standard axes position
box off;
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_1_PFC_example_raster.png'), 'Resolution', 300);
close(fig);

%% Panel 2: PFC heatmap
fprintf('Panel 2: PFC heatmap...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
imagesc(g.timeaxis,1:size(psth_spx_zscore_PFC,1),psth_spx_zscore_PFC);
caxis(g.clim);
colormap(gca, g.colors.Heatmap);
hold on;
xticks([-g.pre_time 0 g.pre_time]);
xticklabels({num2str(-g.pre_time), '0', num2str(g.pre_time)});
yticks([1 size(psth_spx_zscore_PFC,1)]);
% ylabel removed for cleaner panels
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Same position as all other panels
box off;
plot([0 0], [0.5 size(psth_spx_zscore_PFC,1)+0.5], '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_2_PFC_heatmap.png'), 'Resolution', 300);
close(fig);

%% Panel 3: PFC firing rate
fprintf('Panel 3: PFC firing rate...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
time = g.timeaxis(2:end);
meanData_PFC = mean(psth_spx_PFC, 1)/g.bin_time/size(psth_spx_PFC,1);
meanData_PFC = smoothdata(meanData_PFC, 2, 'sgolay', g.smoothvalue);
semData_PFC = std(psth_spx_PFC, 0, 1)/g.bin_time/size(psth_spx_PFC,1) / sqrt(size(psth_spx_PFC,1));
semData_PFC = smoothdata(semData_PFC, 2, 'sgolay', g.smoothvalue);
upper_pfc = meanData_PFC + semData_PFC;
lower_pfc = meanData_PFC - semData_PFC;
hold on;
fill([time, fliplr(time)], [upper_pfc, fliplr(lower_pfc)], [0.8 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(time, meanData_PFC, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 3);
xlim([-g.pre_time g.post_time]);
ylim([-5 60]);
xticks([-g.pre_time 0 g.pre_time]);
xticklabels({num2str(-g.pre_time), '0', num2str(g.pre_time)});
% ylabel removed for cleaner panels
yticks([0 30 60]);
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Set standard axes position
box off;
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_3_PFC_firing_rate.png'), 'Resolution', 300);
close(fig);

%% Panel 4: PFC responsive percentage
fprintf('Panel 4: PFC responsiveness...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
b = bar([1], [pct_nonresponsive_PFC; pct_responsive_PFC]', 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.5);
b(1).CData = [0.85 0.85 0.85];
b(2).CData = [0.8 0.2 0.2];
ylim([0 100]);
% ylabel removed for cleaner panels
yticks([0 50 100]);
set(gca, 'XTickLabel', {''}, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Set standard axes position
box off;
if round(pct_nonresponsive_PFC) > 0
    text(1, pct_nonresponsive_PFC/2, sprintf('%.0f%%', pct_nonresponsive_PFC), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k', 'FontName', g.fontName);
end
if round(pct_responsive_PFC) > 0
    text(1, pct_nonresponsive_PFC + pct_responsive_PFC/2, sprintf('%.0f%%', pct_responsive_PFC), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'w', 'FontName', g.fontName);
end
text(0.80, 0.45, sprintf('(n=%d)', n_total_PFC), 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', g.fontSize2, 'FontWeight', 'normal', 'FontName', g.fontName);
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_4_PFC_responsiveness.png'), 'Resolution', 300);
close(fig);

%% Panel 5: DMS example raster
fprintf('Panel 5: DMS example raster...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
spike_times_DMS = g.cell_metrics.spikes.times{g.example_cellID_DMS};
stimulus_times_DMS = g.cell_metrics.general.shockTTL{g.example_cellID_DMS};
stimulus_times_DMS = stimulus_times_DMS(1:min(80, length(stimulus_times_DMS)));
hold on;
for trial = 1:length(stimulus_times_DMS)
    aligned_spikes = spike_times_DMS - stimulus_times_DMS(trial);
    valid_spikes = aligned_spikes(aligned_spikes >= -g.pre_time & aligned_spikes <= g.post_time);
    scatter(valid_spikes, trial * ones(size(valid_spikes)), g.markerSize, 'k', 'filled');
end
xlim([-g.pre_time g.post_time]);
xticks([-g.pre_time 0 g.post_time]);
xticklabels({num2str(-g.pre_time), '0', num2str(g.pre_time)});
% xlabel and ylabel removed for cleaner panels
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Set standard axes position
box off;
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_5_DMS_example_raster.png'), 'Resolution', 300);
close(fig);

%% Panel 6: DMS heatmap
fprintf('Panel 6: DMS heatmap...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
imagesc(g.timeaxis,1:size(psth_spx_zscore_DMS,1),psth_spx_zscore_DMS);
caxis(g.clim);
colormap(gca, g.colors.Heatmap);
hold on;
% xlabel removed for cleaner panels
xticks([-g.pre_time 0 g.pre_time]);
xticklabels({num2str(-g.pre_time), '0', num2str(g.pre_time)});
yticks([1 size(psth_spx_zscore_DMS,1)]);
% ylabel removed for cleaner panels
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Same position as all other panels
box off;
plot([0 0], [0.5 size(psth_spx_zscore_DMS,1)+0.5], '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_6_DMS_heatmap.png'), 'Resolution', 300);
close(fig);

%% Panel 7: DMS firing rate
fprintf('Panel 7: DMS firing rate...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
meanData_DMS = mean(psth_spx_DMS, 1)/g.bin_time/size(psth_spx_DMS,1);
meanData_DMS = smoothdata(meanData_DMS, 2, 'sgolay', g.smoothvalue);
semData_DMS = std(psth_spx_DMS, 0, 1)/g.bin_time/size(psth_spx_DMS,1) / sqrt(size(psth_spx_DMS,1));
semData_DMS = smoothdata(semData_DMS, 2, 'sgolay', g.smoothvalue);
upper_dms = meanData_DMS + semData_DMS;
lower_dms = meanData_DMS - semData_DMS;
hold on;
fill([time, fliplr(time)], [upper_dms, fliplr(lower_dms)], [0.8 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(time, meanData_DMS, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 3);
xlim([-g.pre_time g.post_time]);
ylim([-5 60]);
xticks([-g.pre_time 0 g.pre_time]);
xticklabels({num2str(-g.pre_time), '0', num2str(g.pre_time)});
% xlabel and ylabel removed for cleaner panels
yticks([0 30 60]);
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Set standard axes position
box off;
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_7_DMS_firing_rate.png'), 'Resolution', 300);
close(fig);

%% Panel 8: DMS responsive percentage
fprintf('Panel 8: DMS responsiveness...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
b = bar([1], [pct_nonresponsive_DMS; pct_responsive_DMS]', 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.5);
b(1).CData = [0.85 0.85 0.85];
b(2).CData = [0.8 0.2 0.2];
ylim([0 100]);
% ylabel removed for cleaner panels
yticks([0 50 100]);
% xlabel removed for cleaner panels
set(gca, 'XTickLabel', {''}, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Set standard axes position
box off;
text(1, pct_nonresponsive_DMS/2, sprintf('%.0f%%', pct_nonresponsive_DMS), ...
    'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k', 'FontName', g.fontName);
text(1, pct_nonresponsive_DMS + pct_responsive_DMS/2, sprintf('%.0f%%', pct_responsive_DMS), ...
    'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'w', 'FontName', g.fontName);
text(0.80, 0.45, sprintf('(n=%d)', n_total_DMS), 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', g.fontSize2, 'FontWeight', 'normal', 'FontName', g.fontName);
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_8_DMS_responsiveness.png'), 'Resolution', 300);
close(fig);

fprintf('\n========================================\n');
fprintf('Main figure panels complete!\n');
fprintf('========================================\n\n');

%% SUPPLEMENTARY FIGURE - Optotagging panels

% Limit optoTTL if specified
if isnumeric(g.opto_ttl_limit)
    g.cell_metrics_opto = g.cell_metrics;
    if length(g.opto_ttl_limit) == 1
        for ii = 1:length(g.cell_metrics_opto.general.optoTTL)
            if length(g.cell_metrics_opto.general.optoTTL{ii}) > g.opto_ttl_limit
                g.cell_metrics_opto.general.optoTTL{ii} = g.cell_metrics_opto.general.optoTTL{ii}(1:g.opto_ttl_limit);
            end
        end
    elseif length(g.opto_ttl_limit) == 2
        for ii = 1:length(g.cell_metrics_opto.general.optoTTL)
            ttl_length = length(g.cell_metrics_opto.general.optoTTL{ii});
            start_idx = min(g.opto_ttl_limit(1), ttl_length);
            end_idx = min(g.opto_ttl_limit(2), ttl_length);
            if start_idx <= end_idx && start_idx >= 1
                g.cell_metrics_opto.general.optoTTL{ii} = g.cell_metrics_opto.general.optoTTL{ii}(start_idx:end_idx);
            else
                g.cell_metrics_opto.general.optoTTL{ii} = [];
            end
        end
    end
else
    g.cell_metrics_opto = g.cell_metrics;
end

% Compute optotagging
psth_spx_opto =  BAfc_psth_spx('cell_metrics', g.cell_metrics_opto, 'ttl', 'optoTTL', 'pre_time', g.optopre, 'post_time', g.optopost, 'bin_time', g.optobin);
psth_spx_zscore_opto = zscore(psth_spx_opto,0,2);
psth_spx_zscore_opto = smoothdata(psth_spx_zscore_opto,2,'sgolay', 3);
psth_spx_zscore_opto_PFC = psth_spx_zscore_opto(idx_PFC,:);
psth_spx_zscore_opto_DMS = psth_spx_zscore_opto(idx_DMS,:);
psth_spx_zscore_opto_PFC = psth_spx_zscore_opto_PFC(order_PFC,:);
psth_spx_zscore_opto_DMS = psth_spx_zscore_opto_DMS(order_DMS,:);
g.clim_opto = [min(psth_spx_zscore_opto(:)) max(psth_spx_zscore_opto(:))];

%% Panel S1: PFC optotagging raster
fprintf('Panel S1: PFC opto raster...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
spike_times_PFC = g.cell_metrics_opto.spikes.times{g.example_cellID_PFC_opto};
stimulus_times_PFC_opto = g.cell_metrics_opto.general.optoTTL{g.example_cellID_PFC_opto};
hold on;
fill([0 0.01 0.01 0], [0 0 length(stimulus_times_PFC_opto)+1 length(stimulus_times_PFC_opto)+1], [0.5 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
for trial = 1:length(stimulus_times_PFC_opto)
    aligned_spikes = spike_times_PFC - stimulus_times_PFC_opto(trial);
    valid_spikes = aligned_spikes(aligned_spikes >= -g.optopre & aligned_spikes <= g.optopost);
    scatter(valid_spikes, trial * ones(size(valid_spikes)), g.markerSize, 'k', 'filled');
end
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
xlim([-g.optopre g.optopost]);
ylim([0 length(stimulus_times_PFC_opto)]);
xticks([-g.optopre 0 g.optopost]);
xticklabels({num2str(-g.optopre*1000), '0', num2str(g.optopost*1000)});
% ylabel removed for cleaner panels
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Set standard axes position
box off;
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_S1_PFC_opto_raster.png'), 'Resolution', 300);
close(fig);

%% Panel S2: PFC optotagging heatmap
fprintf('Panel S2: PFC opto heatmap...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
imagesc(g.optotimeaxis,1:size(psth_spx_zscore_opto_PFC,1),psth_spx_zscore_opto_PFC);
caxis(g.clim_opto);
colormap(gca, g.colors.Heatmap);
xlim([-g.optopre g.optopre]);  % Set consistent x-axis limits
hold on;
plot([0 0], [0.5 size(psth_spx_zscore_opto_PFC,1)+0.5], '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
xticks([-g.optopre 0 g.optopre]);
xticklabels({num2str(-g.optopre*1000), '0', num2str(g.optopre*1000)});
yticks([1 size(psth_spx_zscore_opto_PFC,1)]);
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Same position as all other panels
box off;
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_S2_PFC_opto_heatmap.png'), 'Resolution', 300);
close(fig);

%% Panel S3: DMS optotagging raster
fprintf('Panel S3: DMS opto raster...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
spike_times_DMS = g.cell_metrics_opto.spikes.times{g.example_cellID_DMS_opto};
stimulus_times_DMS_opto = g.cell_metrics_opto.general.optoTTL{g.example_cellID_DMS_opto};
hold on;
fill([0 0.01 0.01 0], [0 0 length(stimulus_times_DMS_opto)+1 length(stimulus_times_DMS_opto)+1], [0.5 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
for trial = 1:length(stimulus_times_DMS_opto)
    aligned_spikes = spike_times_DMS - stimulus_times_DMS_opto(trial);
    valid_spikes = aligned_spikes(aligned_spikes >= -g.optopre & aligned_spikes <= g.optopost);
    scatter(valid_spikes, trial * ones(size(valid_spikes)), g.markerSize, 'k', 'filled');
end
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
xlim([-g.optopre g.optopost]);
ylim([0 length(stimulus_times_DMS_opto)]);
xticks([-g.optopre 0 g.optopost]);
xticklabels({num2str(-g.optopre*1000), '0', num2str(g.optopost*1000)});
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Set standard axes position
box off;
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_S3_DMS_opto_raster.png'), 'Resolution', 300);
close(fig);

%% Panel S4: DMS optotagging heatmap
fprintf('Panel S4: DMS opto heatmap...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, panel_width, panel_height], 'Color', 'w', 'Visible', 'off');
imagesc(g.optotimeaxis,1:size(psth_spx_zscore_opto_DMS,1),psth_spx_zscore_opto_DMS);
caxis(g.clim_opto);
colormap(gca, g.colors.Heatmap);
xlim([-g.optopre g.optopre]);  % Set consistent x-axis limits
hold on;
plot([0 0], [0.5 size(psth_spx_zscore_opto_DMS,1)+0.5], '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
xticks([-g.optopre 0 g.optopre]);
xticklabels({num2str(-g.optopre*1000), '0', num2str(g.optopre*1000)});
yticks([1 size(psth_spx_zscore_opto_DMS,1)]);
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out', 'FontName', g.fontName);
set(gca, 'Position', axes_pos_standard);  % Same position as all other panels
box off;
exportgraphics(gcf, fullfile(g.outputFolder, 'Panel_S4_DMS_opto_heatmap.png'), 'Resolution', 300);
close(fig);

fprintf('\n========================================\n');
fprintf('Supplementary panels complete!\n');
fprintf('========================================\n\n');

%% Colorbar 1: Main figure heatmaps (US response)
fprintf('Creating colorbar 1 (US response)...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, 3, panel_height], 'Color', 'w', 'Visible', 'off');
ax = axes('Visible', 'off');
colormap(ax, g.colors.Heatmap);
set(ax, 'CLim', g.clim);
cb = colorbar(ax, 'Location', 'west');
cb.LineWidth = g.axisLineWidth;
set(cb, 'FontName', g.fontName, 'FontSize', g.fontSize2);
ylabel(cb, 'Z-score', 'FontSize', g.fontSize2, 'FontName', g.fontName);
exportgraphics(gcf, fullfile(g.outputFolder, 'Colorbar_1_US_response.png'), 'Resolution', 300);
close(fig);

%% Colorbar 2: Supplementary heatmaps (Opto response)
fprintf('Creating colorbar 2 (Opto response)...\n');
fig = figure('Units', 'centimeters', 'Position', [0, 0, 3, panel_height], 'Color', 'w', 'Visible', 'off');
ax = axes('Visible', 'off');
colormap(ax, g.colors.Heatmap);
set(ax, 'CLim', g.clim_opto);
cb = colorbar(ax, 'Location', 'west');
cb.LineWidth = g.axisLineWidth;
set(cb, 'FontName', g.fontName, 'FontSize', g.fontSize2);
ylabel(cb, 'Z-score', 'FontSize', g.fontSize2, 'FontName', g.fontName);
exportgraphics(gcf, fullfile(g.outputFolder, 'Colorbar_2_Opto_response.png'), 'Resolution', 300);
close(fig);

fprintf('\n========================================\n');
fprintf('ALL PANELS SAVED TO:\n%s\n\n', g.outputFolder);
fprintf('14 PNG files created (8 main + 4 supplementary + 2 colorbars)\n');
fprintf('========================================\n\n');
