%% function Gergo_figure_1
% exportgraphics(gcf, 'Gergo.tiff', 'Resolution', 300);

% [data, timestamps, info] = load_open_ephys_data('KE210611_02_R_all_channels.events')

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
% Font sizes scaled up for 50x20cm figure that will be scaled down in Word
% When figure is scaled to ~18cm width in Word, fonts will appear at proper size
g.fontSize2 = 22;  % Axis labels and ticks (will scale to ~6-7pt in Word)
g.fontSize1 = 25;  % Title font (will scale to ~8pt in Word)
g.xlinewidth = 3;  % Stimulus line width (will scale to ~1pt in Word)
g.axisLineWidth = 2;  % Axis line width (will scale to ~0.7pt in Word)
g.markerSize = 8;  % Raster plot marker size (will scale to ~3pt in Word)
g.optopre = 0.02;
g.optopost = 0.02;
g.optobin = 0.001;
g.optotimeaxis = -g.optopre:g.optobin:g.optopost;

% Example neuron cell IDs for raster plots (main figure - shock response)
g.example_cellID_PFC = 4;
g.example_cellID_DMS = 27; % 26 27 28 29 good

% Example neuron cell IDs for supplementary figure (optotagging)
g.example_cellID_PFC_opto = 4;
g.example_cellID_DMS_opto = 16; % 26 27 28 29 good

% Optotagging TTL limit (use number, [start end], or 'all')
%g.opto_ttl_limit = 500;  % Use first 500 TTLs
g.opto_ttl_limit = [500 1000];  % Use TTLs from 500 to 1000
% g.opto_ttl_limit = 'all';  % Use all TTLs


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
        cellIDs = cellIDs(~[cellIDs.isdir]); % remove directories
        cellIDs = cellIDs(contains({cellIDs.name},'GR'));
        for jj = 1:size(cellIDs,1)
            spikes = load(cellIDs(jj).name, 'TS');
            ttl = load('shocKTTL.mat', 'shockTTL');
            optofilename = dir('*events*');
            [~, timestamps, ~] = load_open_ephys_data(optofilename.name);
            timestamps = timestamps(1:2:end); % keep onsets only;
            opto_ttl = setdiff(timestamps, ttl.shockTTL);
            g.cell_metrics.general.optoTTL{end+1} = opto_ttl;
            g.cell_metrics.spikes.times{end+1} = spikes.TS/10000; % MClust TS must be divided by 10k
            g.cell_metrics.projection{end+1} = 'PFC';
            g.cell_metrics.general.shockTTL{end+1} = ttl.shockTTL;
            g.cell_metrics.cellID{end+1} = size(g.cell_metrics.cellID,2)+1; 
        end
    end
    % Load DMS projection neurons
    for ii = 1:size(g.fullPaths_DMS,2)
        cd(g.fullPaths_DMS{ii})
        cellIDs = dir;
        cellIDs = cellIDs(~[cellIDs.isdir]); % remove directories
        cellIDs = cellIDs(contains({cellIDs.name},'GR'));
        for jj = 1:size(cellIDs,1)
            spikes = load(cellIDs(jj).name, 'TS');
            ttl = load('shocKTTL.mat', 'shockTTL');
            optofilename = dir('*events*');
            [~, timestamps, ~] = load_open_ephys_data(optofilename.name);
            timestamps = timestamps(1:2:end); % keep onsets only;
            opto_ttl = setdiff(timestamps, ttl.shockTTL);
            g.cell_metrics.general.optoTTL{end+1} = opto_ttl;
            g.cell_metrics.spikes.times{end+1} = spikes.TS/10000; % MClust TS must be divided by 10k
            g.cell_metrics.projection{end+1} = 'DMS';
            g.cell_metrics.general.shockTTL{end+1} = ttl.shockTTL;
            g.cell_metrics.cellID{end+1} = size(g.cell_metrics.cellID,2)+1;
        end
    end
end
clearvars -except g

%% PSTH

psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'shockTTL', 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

% Z-score based on baseline period only (for US response heatmaps/histograms)
baseline_bins = 1:round(g.pre_time/g.bin_time);
baseline_mean = mean(psth_spx(:, baseline_bins), 2);
baseline_std = std(psth_spx(:, baseline_bins), 0, 2);
psth_spx_zscore = (psth_spx - baseline_mean) ./ baseline_std;
%psth_spx_zscore = zscore(psth_spx,[],2);
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
    % Check for excitation
    if max(datasegment) >= g.exctestvalue
        onset_latency(ii) = find(datasegment >= g.exctestvalue, 1, 'first') * g.bin_time * 1000; % in ms
    % Check for inhibition
    elseif min(datasegment) <= -g.inhtestvalue
        onset_latency(ii) = find(datasegment <= -g.inhtestvalue, 1, 'first') * g.bin_time * 1000; % in ms
    else
        % Non-responsive neurons: assign large value to sort last
        onset_latency(ii) = 999999;
    end
end

% Sort by onset latency (earliest first)
[~, order_PFC] = sort(onset_latency(idx_PFC), 'ascend');
[~, order_DMS] = sort(onset_latency(idx_DMS), 'ascend');

psth_spx_zscore_PFC = psth_spx_zscore_PFC(order_PFC,:);
psth_spx_zscore_DMS = psth_spx_zscore_DMS(order_DMS,:);
psth_spx_PFC = psth_spx_PFC(order_PFC,:);
psth_spx_DMS = psth_spx_DMS(order_DMS,:);

% Set color limits based on percentiles to avoid extreme values dominating.
g.clim = [min(psth_spx_zscore(:)) max(psth_spx_zscore(:))];
g.clim = [prctile(psth_spx_zscore(:), 0.5) prctile(psth_spx_zscore(:), 99.5)];


% Main figure - Large size for high resolution export
% Figure is 50x20cm for export at 300 DPI with proper resolution
% When inserted into Word and scaled to ~18-20cm width, fonts will be readable
fig1 = figure('Units', 'centimeters', 'Position', [2, 2, 60, 25], 'Color', 'w');
fig1.PaperPositionMode = 'auto';
fig1.PaperUnits = 'centimeters';
fig1.PaperSize = [60, 25];
tiledlayout(fig1, 2,4,'TileSpacing', 'tight', 'Padding', 'compact')

% PFC example raster
ax1 = nexttile;
spike_times_PFC = g.cell_metrics.spikes.times{g.example_cellID_PFC};
stimulus_times_PFC = g.cell_metrics.general.shockTTL{g.example_cellID_PFC};
stimulus_times_PFC = stimulus_times_PFC(1:min(80, length(stimulus_times_PFC)));
hold on;
for trial = 1:length(stimulus_times_PFC)
    aligned_spikes = spike_times_PFC - stimulus_times_PFC(trial);
    valid_spikes = aligned_spikes(aligned_spikes >= -g.pre_time & aligned_spikes <= g.post_time);
    scatter(valid_spikes, trial * ones(size(valid_spikes)), g.markerSize, 'k', 'filled');
end
xlim([-g.pre_time g.post_time])
xticks([-g.pre_time 0 g.post_time]);
ylabel({['\bf\fontsize{' num2str(g.fontSize2+8) '}BA→PFC']; ['\rm\fontsize{' num2str(g.fontSize2) '}Trials']})
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
title('Example unit', 'FontSize', g.fontSize1, 'FontWeight', 'Normal');
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
box off;
text(-0.25, 1.1, 'A', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');
% Add lightning symbol on top of everything else
annotation('textbox',  [0.14, 0.925, 0.01, 0.01], ...
    'String', '\bf\fontsize{35}\color{red}⚡', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FitBoxToText', 'off');

% PFC heatmap
ax2 = nexttile;
imagesc(g.timeaxis,1:size(psth_spx_zscore_PFC,1),psth_spx_zscore_PFC);
clim(g.clim)
colormap(gca, g.colors.Heatmap);
hold on;
xticks([-g.pre_time 0 g.pre_time]);
xticklabels({num2str(-g.pre_time), '0', num2str(g.pre_time)});
yticks([1 size(psth_spx_zscore_PFC,1)])
title('US response', 'FontSize', g.fontSize1, 'FontWeight', 'Normal');
ylabel('Neuron #', 'FontSize', g.fontSize2)
cb = colorbar('eastoutside', 'FontSize', g.fontSize2-1);
cb.LineWidth = g.axisLineWidth;
ylabel(cb, 'Z-score', 'FontSize', g.fontSize2);
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
box off;
plot([0 0], [0.5 size(psth_spx_zscore_PFC,1)+0.5], '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
text(-0.25, 1.1, 'B', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');
% Add lightning symbol on top of everything else
annotation('textbox', [0.388, 0.925, 0.01, 0.01], ...
    'String', '\bf\fontsize{35}\color{red}⚡', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FitBoxToText', 'off');

% PFC firing rate histogram
ax3 = nexttile;
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
xlim([-g.pre_time g.post_time])
ylim([-5 60])
xticks([-g.pre_time 0 g.pre_time]);
xticklabels({num2str(-g.pre_time), '0', num2str(g.pre_time)});
ylabel('Firing rate (Hz)', 'FontSize', g.fontSize2);
yticks([0 30 60])
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
title('Population activity', 'FontSize', g.fontSize1, 'FontWeight', 'Normal');
box off;
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
text(-0.25, 1.1, 'C', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');
% Add lightning symbol on top of everything else
annotation('textbox', [0.636, 0.925, 0.01, 0.01], ...
    'String', '\bf\fontsize{35}\color{red}⚡', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FitBoxToText', 'off');

% PFC responsive neurons percentage
ax4 = nexttile;
n_responsive_PFC = sum(onset_latency(idx_PFC) < 999999);
n_total_PFC = sum(idx_PFC);
pct_responsive_PFC = (n_responsive_PFC / n_total_PFC) * 100;
pct_nonresponsive_PFC = 100 - pct_responsive_PFC;

% Flip: non-responsive at bottom, responsive on top
b = bar([1], [pct_nonresponsive_PFC; pct_responsive_PFC]', 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.5);
b(1).CData = [0.85 0.85 0.85];  % b(1) = non-responsive (gray)
b(2).CData = [0.8 0.2 0.2];     % b(2) = responsive (red)
ylim([0 100])
ylabel('Percentage (%)', 'FontSize', g.fontSize2)
yticks([0 50 100])
set(gca, 'XTickLabel', {''}, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
title('Responsiveness to US', 'FontSize', g.fontSize1, 'FontWeight', 'Normal');
box off;
text(-0.25, 1.1, 'D', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');
% Only show non-responsive percentage if it would display as non-zero
if round(pct_nonresponsive_PFC) > 0
    text(1, pct_nonresponsive_PFC/2, sprintf('%.0f%%', pct_nonresponsive_PFC), ...
        'HorizontalAlignment', 'center', 'FontSize', g.fontSize2-2, 'FontWeight', 'bold', 'Color', 'k');
end
% Only show responsive percentage if it would display as non-zero
if round(pct_responsive_PFC) > 0
    text(1, pct_nonresponsive_PFC + pct_responsive_PFC/2, sprintf('%.0f%%', pct_responsive_PFC), ...
        'HorizontalAlignment', 'center', 'FontSize', g.fontSize2-2, 'FontWeight', 'bold', 'Color', 'w');
end
lg = legend(b, {'Non-resp.', 'Resp.'}, 'Location', 'northeast', 'FontSize', g.fontSize2-2);
legend('boxoff');
lg.Position(1) = lg.Position(1) + 0.03;  % Move legend slightly right but keep it close to bar
lg.Direction = 'reverse';
% Add total n below legend, inside the panel
text(0.85, 0.55, sprintf('(n=%d)', n_total_PFC), 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', g.fontSize2-1, 'FontWeight', 'normal');

% Row 2: DMS
% DMS example raster
ax5 = nexttile;
spike_times_DMS = g.cell_metrics.spikes.times{g.example_cellID_DMS};
stimulus_times_DMS = g.cell_metrics.general.shockTTL{g.example_cellID_DMS};
stimulus_times_DMS = stimulus_times_DMS(1:min(80, length(stimulus_times_DMS)));
hold on;
for trial = 1:length(stimulus_times_DMS)
    aligned_spikes = spike_times_DMS - stimulus_times_DMS(trial);
    valid_spikes = aligned_spikes(aligned_spikes >= -g.pre_time & aligned_spikes <= g.post_time);
    scatter(valid_spikes, trial * ones(size(valid_spikes)), g.markerSize, 'k', 'filled');
end
xlim([-g.pre_time g.post_time])
xticks([-g.pre_time 0 g.post_time]);
xticklabels({num2str(-g.pre_time), '0', num2str(g.pre_time)});
xlabel('Time (s)', 'FontSize', g.fontSize2)
ylabel({['\bf\fontsize{' num2str(g.fontSize2+8) '}BA→DMS']; ['\rm\fontsize{' num2str(g.fontSize2) '}Trials']})
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
title('Example unit', 'FontSize', g.fontSize1, 'FontWeight', 'Normal');
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
box off;
text(-0.25, 1.1, 'E', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');
% Add lightning symbol on top of everything else
annotation('textbox', [0.14, 0.45, 0.01, 0.01], ...
    'String', '\bf\fontsize{35}\color{red}⚡', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FitBoxToText', 'off');

% DMS heatmap
ax6 = nexttile;
imagesc(g.timeaxis,1:size(psth_spx_zscore_DMS,1),psth_spx_zscore_DMS);
clim(g.clim)
colormap(gca, g.colors.Heatmap);
hold on;
xlabel('Time (s)', 'FontSize', g.fontSize2)
xticks([-g.pre_time 0 g.pre_time]);
xticklabels({num2str(-g.pre_time), '0', num2str(g.pre_time)});
yticks([1 size(psth_spx_zscore_DMS,1)])
title('US response', 'FontSize', g.fontSize1, 'FontWeight', 'Normal');
ylabel('Neuron #', 'FontSize', g.fontSize2)
cb = colorbar('eastoutside', 'FontSize', g.fontSize2-1);
cb.LineWidth = g.axisLineWidth;
ylabel(cb, 'Z-score', 'FontSize', g.fontSize2);
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
box off;
plot([0 0], [0.5 size(psth_spx_zscore_DMS,1)+0.5], '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
text(-0.25, 1.1, 'F', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');
% Add lightning symbol on top of everything else
annotation('textbox', [0.388, 0.45, 0.01, 0.01], ...
    'String', '\bf\fontsize{35}\color{red}⚡', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FitBoxToText', 'off');

% DMS firing rate histogram
ax7 = nexttile;
time = g.timeaxis(2:end);
meanData_DMS = mean(psth_spx_DMS, 1)/g.bin_time/size(psth_spx_DMS,1);
meanData_DMS = smoothdata(meanData_DMS, 2, 'sgolay', g.smoothvalue);
semData_DMS = std(psth_spx_DMS, 0, 1)/g.bin_time/size(psth_spx_DMS,1) / sqrt(size(psth_spx_DMS,1));
semData_DMS = smoothdata(semData_DMS, 2, 'sgolay', g.smoothvalue);
upper_dms = meanData_DMS + semData_DMS;
lower_dms = meanData_DMS - semData_DMS;
hold on;
fill([time, fliplr(time)], [upper_dms, fliplr(lower_dms)], [0.8 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(time, meanData_DMS, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 3);
xlim([-g.pre_time g.post_time])
ylim([-5 60])
xticks([-g.pre_time 0 g.pre_time]);
xticklabels({num2str(-g.pre_time), '0', num2str(g.pre_time)});
xlabel('Time (s)', 'FontSize', g.fontSize2);
ylabel('Firing rate (Hz)', 'FontSize', g.fontSize2);
yticks([0 30 60])
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
title('Population activity', 'FontSize', g.fontSize1, 'FontWeight', 'Normal');
box off;
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
text(-0.25, 1.1, 'G', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');
% Add lightning symbol on top of everything else
annotation('textbox', [0.636, 0.45, 0.01, 0.01], ...
    'String', '\bf\fontsize{35}\color{red}⚡', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FitBoxToText', 'off');

% DMS responsive neurons percentage
ax8 = nexttile;
n_responsive_DMS = sum(onset_latency(idx_DMS) < 999999);
n_total_DMS = sum(idx_DMS);
pct_responsive_DMS = (n_responsive_DMS / n_total_DMS) * 100;
pct_nonresponsive_DMS = 100 - pct_responsive_DMS;

% Flip: non-responsive at bottom, responsive on top
b = bar([1], [pct_nonresponsive_DMS; pct_responsive_DMS]', 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.5);
b(1).CData = [0.85 0.85 0.85];  % b(1) = non-responsive (gray)
b(2).CData = [0.8 0.2 0.2];     % b(2) = responsive (red)
ylim([0 100])
ylabel('Percentage (%)', 'FontSize', g.fontSize2)
yticks([0 50 100])
xlabel('Time (s)', 'FontSize', g.fontSize2);
set(gca, 'XTickLabel', {''}, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
title('Responsiveness to US', 'FontSize', g.fontSize1, 'FontWeight', 'Normal');
box off;
text(-0.25, 1.1, 'H', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');
text(1, pct_nonresponsive_DMS/2, sprintf('%.0f%%', pct_nonresponsive_DMS), ...
    'HorizontalAlignment', 'center', 'FontSize', g.fontSize2-2, 'FontWeight', 'bold', 'Color', 'k');
text(1, pct_nonresponsive_DMS + pct_responsive_DMS/2, sprintf('%.0f%%', pct_responsive_DMS), ...
    'HorizontalAlignment', 'center', 'FontSize', g.fontSize2-2, 'FontWeight', 'bold', 'Color', 'w');
lg = legend(b, {'Non-resp.', 'Resp.'}, 'Location', 'northeast', 'FontSize', g.fontSize2-2);
legend('boxoff');
lg.Position(1) = lg.Position(1) + 0.03;  % Move legend slightly right but keep it close to bar
lg.Direction = 'reverse';
% Add total n below legend, inside the panel
text(0.85, 0.55, sprintf('(n=%d)', n_total_DMS), 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', g.fontSize2-1, 'FontWeight', 'normal');

%% Supplementary figure: Optotagging and Z-score histograms
% Nature standard size: 180mm (double column) x ~100mm height for supplementary
% Set Units to cm for consistent sizing across screen and export
fig2 = figure('Units', 'centimeters', 'Position', [2, 2, 35, 25], 'Color', 'w');
fig2.PaperPositionMode = 'auto';
fig2.PaperUnits = 'centimeters';
fig2.PaperSize = [35, 25];
tiledlayout(fig2, 2,2,'TileSpacing', 'tight', 'Padding', 'compact')

% Limit optoTTL if specified
if isnumeric(g.opto_ttl_limit)
    g.cell_metrics_opto = g.cell_metrics;
    if length(g.opto_ttl_limit) == 1
        % Single number: use first N TTLs
        for ii = 1:length(g.cell_metrics_opto.general.optoTTL)
            if length(g.cell_metrics_opto.general.optoTTL{ii}) > g.opto_ttl_limit
                g.cell_metrics_opto.general.optoTTL{ii} = g.cell_metrics_opto.general.optoTTL{ii}(1:g.opto_ttl_limit);
            end
        end
    elseif length(g.opto_ttl_limit) == 2
        % Range [start end]: use TTLs from start to end
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
idx_PFC = strcmp(g.cell_metrics_opto.projection, 'PFC')';
idx_DMS = strcmp(g.cell_metrics_opto.projection, 'DMS')';
psth_spx_zscore_opto_PFC = psth_spx_zscore_opto(idx_PFC,:);
psth_spx_zscore_opto_DMS = psth_spx_zscore_opto(idx_DMS,:);
psth_spx_zscore_opto_PFC = psth_spx_zscore_opto_PFC(order_PFC,:);
psth_spx_zscore_opto_DMS = psth_spx_zscore_opto_DMS(order_DMS,:);
g.clim_opto = [min(psth_spx_zscore_opto(:)) max(psth_spx_zscore_opto(:))];

% Row 1: PFC
% PFC example raster (optotagging)
ax9 = nexttile;
spike_times_PFC = g.cell_metrics_opto.spikes.times{g.example_cellID_PFC_opto};
stimulus_times_PFC_opto = g.cell_metrics_opto.general.optoTTL{g.example_cellID_PFC_opto};
hold on;
% Add shaded area for light pulse (0-10ms)
fill([0 0.01 0.01 0], [0 0 length(stimulus_times_PFC_opto)+1 length(stimulus_times_PFC_opto)+1], [0.5 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
for trial = 1:length(stimulus_times_PFC_opto)
    aligned_spikes = spike_times_PFC - stimulus_times_PFC_opto(trial);
    valid_spikes = aligned_spikes(aligned_spikes >= -g.optopre & aligned_spikes <= g.optopost);
    scatter(valid_spikes, trial * ones(size(valid_spikes)), g.markerSize, 'k', 'filled');
end
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
xlim([-g.optopre g.optopost])
ylim([0 length(stimulus_times_PFC_opto)])  % Set y-axis to actual number of trials
xticks([-g.optopre 0 g.optopost]);
xticklabels({num2str(-g.optopre*1000), '0', num2str(g.optopost*1000)});
ylabel({['\bf\fontsize{' num2str(g.fontSize2+8) '}BA→PFC']; ['\rm\fontsize{' num2str(g.fontSize2) '}Trials']})
title('Optotagging', 'FontSize', g.fontSize1, 'FontWeight', 'Normal')
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
box off;
text(-0.25, 1.1, 'A', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');

% PFC optotagging heatmap
ax10 = nexttile;
imagesc(g.optotimeaxis,1:size(psth_spx_zscore_opto_PFC,1),psth_spx_zscore_opto_PFC);
clim(g.clim_opto)
colormap(gca, g.colors.Heatmap);
hold on;
plot([0 0], [0.5 size(psth_spx_zscore_opto_PFC,1)+0.5], '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
xticks([-g.optopre 0 g.optopre]);
xticklabels({num2str(-g.optopre*1000), '0', num2str(g.optopre*1000)});
yticks([1 size(psth_spx_zscore_opto_PFC,1)])
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
title('Light response', 'FontSize', g.fontSize1, 'FontWeight', 'Normal')
ylabel('Neuron #', 'FontSize', g.fontSize2)
cb = colorbar('eastoutside', 'FontSize', g.fontSize2-1);
cb.LineWidth = g.axisLineWidth;
ylabel(cb, 'Z-score', 'FontSize', g.fontSize2);
text(-0.25, 1.1, 'B', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');

% Row 2: DMS
% DMS example raster (optotagging)
ax12 = nexttile;
spike_times_DMS = g.cell_metrics_opto.spikes.times{g.example_cellID_DMS_opto};
stimulus_times_DMS_opto = g.cell_metrics_opto.general.optoTTL{g.example_cellID_DMS_opto};
hold on;
% Add shaded area for light pulse (0-10ms)
fill([0 0.01 0.01 0], [0 0 length(stimulus_times_DMS_opto)+1 length(stimulus_times_DMS_opto)+1], [0.5 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
for trial = 1:length(stimulus_times_DMS_opto)
    aligned_spikes = spike_times_DMS - stimulus_times_DMS_opto(trial);
    valid_spikes = aligned_spikes(aligned_spikes >= -g.optopre & aligned_spikes <= g.optopost);
    scatter(valid_spikes, trial * ones(size(valid_spikes)), g.markerSize, 'k', 'filled');
end
xline(0, '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
xlim([-g.optopre g.optopost])
ylim([0 length(stimulus_times_DMS_opto)])  % Set y-axis to actual number of trials
xticks([-g.optopre 0 g.optopost]);
xticklabels({num2str(-g.optopre*1000), '0', num2str(g.optopost*1000)});
xlabel('Time (ms)', 'FontSize', g.fontSize2)
ylabel({['\bf\fontsize{' num2str(g.fontSize2+8) '}BA→DMS']; ['\rm\fontsize{' num2str(g.fontSize2) '}Trials']})
title('Optotagging', 'FontSize', g.fontSize1, 'FontWeight', 'Normal')
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
box off;
text(-0.25, 1.1, 'D', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');

% DMS optotagging heatmap
ax13 = nexttile;
imagesc(g.optotimeaxis,1:size(psth_spx_zscore_opto_DMS,1),psth_spx_zscore_opto_DMS);
clim(g.clim_opto)
colormap(gca, g.colors.Heatmap);
hold on;
plot([0 0], [0.5 size(psth_spx_zscore_opto_DMS,1)+0.5], '--', 'Color', 'k', 'LineWidth', g.xlinewidth);
xlabel('Time (ms)', 'FontSize', g.fontSize2)
xticks([-g.optopre 0 g.optopre]);
xticklabels({num2str(-g.optopre*1000), '0', num2str(g.optopre*1000)});
yticks([1 size(psth_spx_zscore_opto_DMS,1)])
set(gca, 'FontSize', g.fontSize2, 'LineWidth', g.axisLineWidth, 'TickDir', 'out');
title('Light response', 'FontSize', g.fontSize1, 'FontWeight', 'Normal')
ylabel('Neuron #', 'FontSize', g.fontSize2)
cb = colorbar('eastoutside', 'FontSize', g.fontSize2-1);
cb.LineWidth = g.axisLineWidth;
ylabel(cb, 'Z-score', 'FontSize', g.fontSize2);
text(-0.25, 1.1, 'D', 'Units', 'normalized', 'FontSize', g.fontSize1+2, 'FontWeight', 'bold');
