% function BAfc_figure_1 
% Plotting experimental design, waveform features and cell type  separation

recordings = {...
    'MD243_kilosort',... 
    'MD250_kilosort',... 
    'MD251_kilosort',... 
    'MD252_kilosort',...
    'MD253_kilosort',...
    'MD254_kilosort',...
    'MD266_kilosort',...
    'MD267_kilosort',...
    'MD268_kilosort',...
    'MD269_kilosort',...
    'MD275_kilosort_waveform800',...
    'MD276_kilosort_waveform800',...
    'MD277_kilosort_waveform800',...
    'MD278_kilosort_waveform800',...
    'MD288_kilosort',...
    'MD289_kilosort',...
    'MD290_kilosort_waveform800',...
    'MD291_kilosort',...
    'MD292_002_kilosort',...
    'MD293_kilosort',...
    'MD294_kilosort',...
    'MD295_kilosort',...
    'MD296_kilosort',...
    'MD297_kilosort'};
g.cell_metrics = BAfc_load_neurons('recordings', recordings);

clearvars -except g 
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 15;
g.fontSize2 = 12;

%% Initiate figure
fig = figure('Position', [400, 100, 1800, 1200]);
t = tiledlayout(fig,16,12,'TileSpacing', 'tight', 'Padding', 'none');

%% (1,1) - Schematic figure
ax1 = nexttile(t,1,[4 4]);
[img, cmap] = imread([g.mainFolder '\drawed_mouse.png']);
if ~isempty(cmap)
    img = ind2rgb(img, cmap);  % Convert to RGB
end
imshow(img, 'Parent', ax1);
title(ax1, 'Experimental setup', 'FontSize', g.fontSize1)
clearvars -except t g % clear variables

%% (1,2) - Example raw traces and PSTH
t2 = tiledlayout(t,2,4, 'TileSpacing', 'tigh');
t2.Layout.Tile = 5;
t2.Layout.TileSpan = [4 4];
ax2_1 = nexttile(t2, [1 4]);
rawMD292_002 = load([g.mainFolder '\MD292_002_kilosort\kilosort25preprocess\temp_wh_ch32.mat']); % cludID = 341 from MD292_002
rawMD292_002 = rawMD292_002.temp_wh_Ch32.values*0.195; % convert to mV
cellID = intersect(find(strcmp(g.cell_metrics.animal, 'MD292')), find(g.cell_metrics.cluID == 341));
ttlMD292_002 = load([g.mainFolder '\MD292_002_kilosort\kilosort25preprocess\TTLsKS.mat']);
ttlMD292_002 = ttlMD292_002.shocks([5 7 36 41 43 47]); % 34, 36, 41 43, 47?
twin = [-0.01 0.04];
fs = 30000;
hold on
for ii = 1:size(ttlMD292_002,1)
    sn_curr = ttlMD292_002(ii,1)*fs;
    timeaxis = linspace(twin(1),twin(2),(twin(2)-twin(1))*fs)*1000;%in ms
    plot(ax2_1,timeaxis, rawMD292_002(round(sn_curr+twin(1)*fs):round(sn_curr+twin(2)*fs-1)), 'Color', g.colors.raw_primary)
end
ylabel('Voltage (mV)', 'FontSize', g.fontSize2)
ylim([-0.5 0.3]);
ax2_1.XColor = 'none';
% Define shaded region (0 to 10 ms)
yLimits = get(ax2_1, 'YLim'); % get y-axis limits
xShade = [0 10 10 0];       % in ms
yShade = [yLimits(1) yLimits(1) yLimits(2) yLimits(2)];
fill(ax2_1, xShade, yShade, [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none') % red with transparency
hold off
title(ax2_1, 'Example data traces', 'FontSize', g.fontSize1)
% PSTH
ax2_2 = nexttile(t2, [1 4]);
% Parameters
binSize = 0.001; % Bin size in seconds (e.g., 10 ms)
window = twin; % Time window around stimulus (before and after)
% Preallocate
spikeTimes = g.cell_metrics.spikes.times{cellID};
stimTimes = load('C:\Users\dmagyar\Desktop\BA_fear_cond\MD292_002_kilosort\kilosort25preprocess\TTLsKS.mat', 'triptest_shocks_only');
stimTimes = stimTimes.triptest_shocks_only;
edges = window(1):binSize:window(2);
counts = zeros(length(stimTimes), length(edges)-1);
% Loop through stimuli
for i = 1:length(stimTimes)
    alignedSpikes = spikeTimes - stimTimes(i); % Align spikes to stim
    counts(i,:) = histcounts(alignedSpikes, edges);
end
% Average across trials
meanCounts = mean(counts, 1);
% Plot PSTH
b = bar(edges(1:end-1)*1000+binSize*1000/2, meanCounts/binSize, 'hist'); % spikes/sec
b.FaceColor =  g.colors.raw_primary;
b.EdgeColor = 'none';  
hold on
yLimits = get(ax2_2, 'YLim'); % get y-axis limits
xShade = [0 10 10 0];       % in ms
yShade = [yLimits(1) yLimits(1) yLimits(2) yLimits(2)];
fill(ax2_2, xShade, yShade, [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none') % red with transparency
hold off
xlabel('Time (ms)', 'FontSize', g.fontSize1);
ylabel('Firing Rate (Hz)', 'FontSize', g.fontSize2);
ylim([0 500]);
ax2_1.Layer = 'top'; % to fix graphical glithc
ax2_2.Layer = 'top'; 
ax2_1.Box = 'off';
ax2_2.Box = 'off';
clearvars -except t g % clear variables

%% (1,3) - Trough to peak distribution
ax3 = nexttile(t, 9, [4 4]);
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);
idx_LA = find(contains(g.cell_metrics.brainRegion, 'LA'));
idx_BA = find(contains(g.cell_metrics.brainRegion, 'BA'));
data = g.cell_metrics.troughToPeak([idx_LA idx_BA])';
gm = fitgmdist(data, 2); % Fit a 2-component Gaussian mixture
% Generate x values for PDF
x = linspace(min(data), max(data), 1000);
y = pdf(gm, x');
% Plot histogram
h = histogram(ax3,data, 'Normalization', 'pdf', 'BinWidth', 0.01);
hold on;
% Overlay the GMM fit
plot(x, y, 'r', 'LineWidth', 2);
xline(0.475, '--r', '475 ms', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle', 'FontSize', g.fontSize1);
hold off;
title('Trough to peak distribution', 'FontSize', g.fontSize1);
xlabel('Trough to peak time (ms)', 'FontSize', g.fontSize1);
ylabel('Probability Density', 'FontSize', g.fontSize1);
ax3.Box = 'off';
clearvars -except t g % clear variables

%% (2,1) - Example PN and IN ISI with ACG and waveform
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);
cellID_1 = intersect(find(strcmp(g.cell_metrics.animal, 'MD289')), find(g.cell_metrics.cluID == 363));
cellID_2 = intersect(find(strcmp(g.cell_metrics.animal, 'MD250')), find(g.cell_metrics.cluID == 321));
spike_timestamps_n = g.cell_metrics.spikes.times{cellID_1}; % example neuron
ISI_n = diff(spike_timestamps_n); 
spike_timestamps_w = g.cell_metrics.spikes.times{cellID_2}; % example neuron
ISI_w = diff(spike_timestamps_w);
num_bins = 50; % Increase the number of bins
bin_edges = logspace(log10(min([ISI_n;ISI_w])), log10(max([ISI_n;ISI_w])), num_bins);   
% PN
ax4_1 = nexttile(t, 49, [4 3]);
histogram(ISI_w, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', g.colors.PN_primary);
xlim([0.001 80])
set(gca, 'XScale', 'log'); % Set x-axis to log scale
xlabel('Time (ms)', 'FontSize', g.fontSize1);
ylim([0 0.1])
ylabel('Probability', 'FontSize', g.fontSize1);
title('Example PN ISI', 'FontSize', g.fontSize1)
p1 = ax4_1.Position;
acg_w = axes('Position', [p1(1)+0.15 p1(2)+0.15 p1(3)*0.25 p1(4)*0.25]);
bar(acg_w,g.cell_metrics.acg.narrow(75:125,cellID_1), 'FaceColor', g.colors.PN_primary, 'EdgeColor', 'none')
waveform_w = axes('Position', [p1(1)+0.15 p1(2)+0.08 p1(3)*0.25 p1(4)*0.25]);
plot(waveform_w,g.cell_metrics.waveforms.filt{cellID_1}, 'Color', g.colors.PN_primary, 'LineWidth', 2);  
waveform_w.XColor = 'none';
waveform_w.YColor = 'none';
acg_w.XColor = 'none';
acg_w.YColor = 'none';
% IN
ax4_2 = nexttile(t, 52, [4 3]);
histogram(ISI_n, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', g.colors.IN_primary);  
xlim([0.001 80])
set(gca, 'XScale', 'log'); % Set x-axis to log scale
xlabel('Time (ms)', 'FontSize', g.fontSize1);
ylim([0 0.1])    
title('Example IN ISI', 'FontSize', g.fontSize1)
ax4_2.YTickLabel = [];
p2 = ax4_2.Position;
acg_n = axes('Position', [p2(1)+0.15 p2(2)+0.15 p2(3)*0.25 p2(4)*0.25]);
bar(acg_n,g.cell_metrics.acg.narrow(75:125,cellID_2), 'FaceColor', g.colors.IN_primary, 'EdgeColor', 'none')
waveform_n = axes('Position', [p2(1)+0.15 p2(2)+0.08 p2(3)*0.25 p2(4)*0.25]);
plot(waveform_n,g.cell_metrics.waveforms.filt{cellID_2}, 'Color', g.colors.IN_primary, 'LineWidth', 2);  
waveform_n.XColor = 'none';
waveform_n.YColor = 'none';
acg_n.XColor = 'none';
acg_n.YColor = 'none';
ax4_1.Box = 'off';
ax4_2.Box = 'off';
clearvars -except t g % clear variables

%% (2,2) - Spike feature accross unit type
ax5 = nexttile(t, 55, [4 3]);
idx_PN = strcmp(g.cell_metrics.putativeCellType,'PN');
idx_IN = strcmp(g.cell_metrics.putativeCellType,'IN');
idx_unknown = strcmp(g.cell_metrics.putativeCellType,'unknown');
hold on
plot(g.cell_metrics.troughToPeak(idx_IN),g.cell_metrics.firingRate(idx_IN), 'o', 'Color',g.colors.IN_primary)
plot(g.cell_metrics.troughToPeak(idx_PN),g.cell_metrics.firingRate(idx_PN), 'o', 'Color',g.colors.PN_primary)
plot(g.cell_metrics.troughToPeak(idx_unknown),g.cell_metrics.firingRate(idx_unknown), 'o', 'Color',[.7 .7 .7])
set(gca, 'YScale', 'log')  % Makes the Y-axis logarithmic
yline(6, '--k', 'LineWidth', 1);
xline(0.475, '--k', 'LineWidth', 1);   
ylabel('Firing rate (Hz)', 'FontSize', g.fontSize1);
xlabel('Trough to peak time (ms)', 'FontSize', g.fontSize1);
title('Spike features accross unit type', 'FontSize', g.fontSize1);
%% (2,3) - Normalized waveforms

ax6 = nexttile(t, 58, [4 3]);
timeaxis = g.cell_metrics.waveforms.time{1};
waveforms = zeros(size(g.cell_metrics.waveforms.filt,2),48);
for ii = 1:size(g.cell_metrics.waveforms.filt,2)
    if size(g.cell_metrics.waveforms.time{ii},2) == 48
        waveforms(ii,:) = zscore(g.cell_metrics.waveforms.filt{ii});
    elseif size(g.cell_metrics.waveforms.time{ii},2) == 96
        waveforms(ii,:) = zscore(g.cell_metrics.waveforms.filt{ii}(30:77)); % 30:77 is the best, others shift towrads left or right
    end
end
% flipping reverse spikes
for ii = 1:size(waveforms,1)
    if median(waveforms(ii,23:27)) > 1 && strcmp(g.cell_metrics.putativeCellType{ii}, 'Narrow Interneuron')
        waveforms(ii,:) = -waveforms(ii,:);
    end
end
% finding outliers based on Euclidean distance
wide_idx = find(idx_PN == 1);
mean_wide = mean(waveforms(wide_idx,:),1);
narrow_idx = find(idx_IN == 1);
mean_narrow = mean(waveforms(narrow_idx,:),1);

distances = zeros(size(waveforms,1), 1);
for ii = wide_idx
    distances(ii) = sqrt(sum((waveforms(ii, :) - mean_wide).^2));
end
for ii = narrow_idx
    distances(ii) = sqrt(sum((waveforms(ii, :) - mean_narrow).^2));
end    
waveforms(distances>7,:) = NaN;    
% Compute mean and standard deviation
mean_IN = mean(waveforms(idx_IN, 13:48), 1, 'omitnan');
std_IN = std(waveforms(idx_IN, 13:48), 0, 1, 'omitnan');
mean_PN = mean(waveforms(idx_PN, 13:48), 1, 'omitnan');
std_PN = std(waveforms(idx_PN, 13:48), 0, 1, 'omitnan');
% Plot mean waveform with shaded standard deviation for IN
hold on;
fill([timeaxis(13:48), fliplr(timeaxis(13:48))], ...
     [mean_IN + std_IN, fliplr(mean_IN - std_IN)], ...
     g.colors.IN_secondary, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shaded region for SD
plot(timeaxis(13:48), mean_IN, 'Color', g.colors.IN_primary, 'LineWidth', 1.5); % Mean waveform
hold off;
% Plot mean waveform with shaded standard deviation for IN
%figure;
hold on;
fill([timeaxis(13:48), fliplr(timeaxis(13:48))], ...
     [mean_PN + std_PN, fliplr(mean_PN - std_PN)], ...
     g.colors.PN_secondary, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shaded region for SD
plot(timeaxis(13:48), mean_PN, 'Color', g.colors.PN_primary, 'LineWidth', 1.5); % Mean waveform
hold off;     
xlabel('Time (ms)', 'FontSize', g.fontSize1);
ylabel('Z-score', 'FontSize', g.fontSize1)
title('Normalized waveforms accross unit type', 'FontSize', g.fontSize1);
clearvars -except t g % clear variables

%% (3,1) - Example LA trace
ax7 = nexttile(t, 97, [8 4]);
[img, cmap] = imread([g.mainFolder '\example_LA_trace.png']);
if ~isempty(cmap)
    img = ind2rgb(img, cmap);  % Convert to RGB
end
imshow(img, 'Parent', ax7);
title(ax7, 'Example LA trace', 'FontSize', g.fontSize1)
clearvars -except t g % clear variables

%% (3,2) - Example BA trace
ax8 = nexttile(t, 101, [8 4]);
[img, cmap] = imread([g.mainFolder '\example_BA_trace.png']);
if ~isempty(cmap)
    img = ind2rgb(img, cmap);  % Convert to RGB
end
imshow(img, 'Parent', ax8);
title(ax8, 'Example BA trace', 'FontSize', g.fontSize1)
clearvars -except t g % clear variables

%% (3,3) - All traces
ax9 = nexttile(t, 105, [8 4]);
[img, cmap] = imread([g.mainFolder '\all_traces.png']);
if ~isempty(cmap)
    img = ind2rgb(img, cmap);  % Convert to RGB
end
imshow(img, 'Parent', ax9);
title(ax9, 'All recorded traces', 'FontSize', g.fontSize1)
clearvars -except t g % clear variables