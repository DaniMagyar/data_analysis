%function BAfc_figure_1

% recordings = {...
%     'MD243_kilosort',... 
%     'MD250_kilosort',... 
%     'MD251_kilosort',... 
%     'MD252_kilosort',...
%     'MD253_kilosort',...
%     'MD254_kilosort',...
%     'MD266_kilosort',...
%     'MD267_kilosort',...
%     'MD268_kilosort',...
%     'MD269_kilosort',...
%     'MD275_kilosort_waveform800',...
%     'MD276_kilosort_waveform800',...
%     'MD277_kilosort_waveform800',...
%     'MD278_kilosort_waveform800',...
%     'MD288_kilosort',...
%     'MD289_kilosort',...
%     'MD292_002_kilosort',...
%     'MD290_kilosort_waveform800',...
%     'MD291_kilosort',...
%     'MD293_kilosort',...
%     'MD294_kilosort',...
%     'MD295_kilosort',...
%     'MD296_kilosort',...
%     'MD297_kilosort'};
% cell_metrics = BAfc_load_neurons('recordings', recordings);

clearvars -except cell_metrics
mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
colors = BAfc_colors;

%% (1,1) - Schematic figure
fig = figure('Position', [400, 100, 2000, 1200]);
t = tiledlayout(fig,3,3,'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(t);
[img, cmap] = imread([mainFolder '\drawed_mouse.png']);
if ~isempty(cmap)
    img = ind2rgb(img, cmap);  % Convert to RGB
end
imshow(img, 'Parent', ax1);
title(ax1, 'Experimental setup')
clearvars -except cell_metrics t mainFolder colors % clear variables

%% (1,2) - Example raw traces and PSTH
t1 = tiledlayout(t,2,1,'TileSpacing', 'compact', 'Padding', 'compact');
t1.Layout.Tile = 2;
ax2_1 = nexttile(t1);
rawMD292_002 = load([mainFolder '\MD292_002_kilosort\kilosort25preprocess\temp_wh_ch32.mat']); % cludID = 341 from MD292_002
rawMD292_002 = rawMD292_002.temp_wh_Ch32.values*0.195; % convert to mV
cellID = intersect(find(strcmp(cell_metrics.animal, 'MD292')), find(cell_metrics.cluID == 341));
ttlMD292_002 = load([mainFolder '\MD292_002_kilosort\kilosort25preprocess\TTLsKS.mat']);
ttlMD292_002 = ttlMD292_002.shocks([5 7 36 41 43 47]); % 34, 36, 41 43, 47?
twin = [-0.01 0.04];
fs = 30000;
hold on
for ii = 1:size(ttlMD292_002,1)
    sn_curr = ttlMD292_002(ii,1)*fs;
    timeaxis = linspace(twin(1),twin(2),(twin(2)-twin(1))*fs)*1000;%in ms
    plot(ax2_1,timeaxis, rawMD292_002(round(sn_curr+twin(1)*fs):round(sn_curr+twin(2)*fs-1)), 'Color', colors.raw_primary)
end
ylabel('Voltage (mV)', 'FontSize', 12)
ylim([-0.5 0.3]);
ax2_1.XColor = 'none';
% Define shaded region (0 to 10 ms)
yLimits = get(ax2_1, 'YLim'); % get y-axis limits
xShade = [0 10 10 0];       % in ms
yShade = [yLimits(1) yLimits(1) yLimits(2) yLimits(2)];
fill(ax2_1, xShade, yShade, [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none') % red with transparency
hold off
title(ax2_1, 'Example data traces')
% PSTH
ax2_2 = nexttile(t1);
% Parameters
binSize = 0.001; % Bin size in seconds (e.g., 10 ms)
window = twin; % Time window around stimulus (before and after)
% Preallocate
spikeTimes = cell_metrics.spikes.times{cellID};
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
b.FaceColor =  colors.raw_primary;
b.EdgeColor = 'none';  
hold on
yLimits = get(ax2_2, 'YLim'); % get y-axis limits
xShade = [0 10 10 0];       % in ms
yShade = [yLimits(1) yLimits(1) yLimits(2) yLimits(2)];
fill(ax2_2, xShade, yShade, [1 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none') % red with transparency
hold off
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Firing Rate (Hz)', 'FontSize', 12);
ylim([0 500]);
ax2_1.Layer = 'top'; % to fix graphical glithc
ax2_2.Layer = 'top'; 
ax2_1.Box = 'off';
ax2_2.Box = 'off';
clearvars -except cell_metrics t mainFolder colors % clear variables

%% (1,3) - Trough to peak distribution
ax3 = nexttile(t);
cell_metrics = BAfc_putative_cellTypes('cell_metrics', cell_metrics);
idx_LA = find(contains(cell_metrics.brainRegion, 'LA'));
idx_BA = find(contains(cell_metrics.brainRegion, 'BA'));
data = cell_metrics.troughToPeak([idx_LA idx_BA])';
gm = fitgmdist(data, 2); % Fit a 2-component Gaussian mixture
% Generate x values for PDF
x = linspace(min(data), max(data), 1000);
y = pdf(gm, x');
% Plot histogram
h = histogram(ax3,data, 'Normalization', 'pdf', 'BinWidth', 0.01);
hold on;
% Overlay the GMM fit
plot(x, y, 'r', 'LineWidth', 2);
xline(0.475, '--r', '475 ms', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle', 'FontSize', 12);
hold off;
title('Trough to peak distribution');
xlabel('Trough to peak time (ms)');
ylabel('Probability Density');
ax3.Box = 'off';
clearvars -except cell_metrics t mainFolder colors % clear variables

%% (2,1) - Example PN and IN ISI with ACG and waveform
cell_metrics = BAfc_putative_cellTypes('cell_metrics', cell_metrics);
t2 = tiledlayout(t,1,2,'TileSpacing', 'compact', 'Padding', 'compact');
t2.Layout.Tile = 4;
cellID_1 = intersect(find(strcmp(cell_metrics.animal, 'MD250')), find(cell_metrics.cluID == 321));
cellID_2 = intersect(find(strcmp(cell_metrics.animal, 'MD289')), find(cell_metrics.cluID == 363));
spike_timestamps_n = cell_metrics.spikes.times{cellID_1}; % example neuron
ISI_n = diff(spike_timestamps_n); 
spike_timestamps_w = cell_metrics.spikes.times{cellID_2}; % example neuron
ISI_w = diff(spike_timestamps_w);
num_bins = 50; % Increase the number of bins
bin_edges = logspace(log10(min([ISI_n;ISI_w])), log10(max([ISI_n;ISI_w])), num_bins);   
% PN
ax4_1 = nexttile(t2);
histogram(ISI_w, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', colors.PN_primary);
xlim([0.001 80])
set(gca, 'XScale', 'log'); % Set x-axis to log scale
xlabel('Time (ms)', 'FontSize', 12);
ylim([0 0.1])
ylabel('Probability', 'FontSize', 12);
title('Example PN ISI', 'FontSize', 12)
p1 = t2.Position;
acg_w = axes('Position', [p1(1)+0.08 p1(2)+0.2 p1(3)*0.15 p1(4)*0.2]);
bar(acg_w,cell_metrics.acg.narrow(50:150,cellID_1), 'FaceColor', colors.PN_primary, 'EdgeColor', 'none')
acg_w.XColor = 'none';
acg_w.YColor = 'none';
% IN
ax4_2 = nexttile(t2);
histogram(ISI_n, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', colors.IN_primary);  
xlim([0.001 80])
set(gca, 'XScale', 'log'); % Set x-axis to log scale
xlabel('Time (ms)', 'FontSize', 12);
ylim([0 0.1])    
title('Example IN ISI', 'FontSize', 12)
ax4_2.YTickLabel = [];
p2 = t2.Position;
acg_n = axes('Position', [p2(1)+0.23 p2(2)+0.2 p2(3)*0.15 p2(4)*0.2]);
bar(acg_n,cell_metrics.acg.narrow(50:150,cellID_2), 'FaceColor', colors.IN_primary, 'EdgeColor', 'none')
acg_n.XColor = 'none';
acg_n.YColor = 'none';

ax4_1.Box = 'off';
ax4_2.Box = 'off';