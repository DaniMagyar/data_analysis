%% function Gergo_figure_2

% Anesthetized data

clear all

g.mainFolder = 'C:\Users\dmagyar\Desktop\Gergo';
g.subfolders.LA = 'C:\Users\dmagyar\Desktop\Gergo\LA_urethane Gergo\LA_urethane';
folderList_LA = dir(g.subfolders.LA);
folderList_LA = folderList_LA([folderList_LA.isdir] & ~ismember({folderList_LA.name}, {'.', '..'}));
g.fullPaths_LA = fullfile(g.subfolders.LA, {folderList_LA.name});

g.subfolders.BA = 'C:\Users\dmagyar\Desktop\Gergo\BA_urethane Gergo\BA_urethane';
folderList_BA = dir(g.subfolders.BA);
folderList_BA = folderList_BA([folderList_BA.isdir] & ~ismember({folderList_BA.name}, {'.', '..'}));
g.fullPaths_BA = fullfile(g.subfolders.BA, {folderList_BA.name});

g.pre_time = 1;
g.post_time = 1;
g.bin_time = 0.005;
g.smoothvalue = 20;
g.timeaxis = -g.pre_time:g.bin_time:g.post_time;
g.colors = BAfc_colors;
g.test_time = 0.6;
g.exctestvalue = 1; % if the mean change is negative, but there is a peak where zscore is above this, than counted as excited
g.inhtestvalue = 1;
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
g.fontSize2 = 12;
g.fontSize1 = 15;
g.xlinewidth = 2;

%% Create cell_metrics
g.cell_metrics.spikes.times = {};
g.cell_metrics.brainRegion = {};
g.cell_metrics.general.shockTTL = {};
g.cell_metrics.cellID = {};
% Load LA neurons
for ii = 1:size(g.fullPaths_LA,2)
    cd(g.fullPaths_LA{ii})
    cellIDs = dir;
    cellIDs = cellIDs(~[cellIDs.isdir]); % remove directories
    cellIDs = cellIDs(contains({cellIDs.name},'GR'));
    for jj = 1:size(cellIDs,1)
        spikes = load(cellIDs(jj).name, 'TS');
        ttl = load('shocKTTL.mat', 'shockTTL');
        g.cell_metrics.spikes.times{end+1} = spikes.TS/10000; % MClust TS must be divided by 10k
        g.cell_metrics.brainRegion{end+1} = 'LA';
        g.cell_metrics.general.shockTTL{end+1} = ttl.shockTTL;
        g.cell_metrics.cellID{end+1} = size(g.cell_metrics.cellID,2)+1; 
    end
end
% Load BA neurons
for ii = 1:size(g.fullPaths_BA,2)
    cd(g.fullPaths_BA{ii})
    cellIDs = dir;
    cellIDs = cellIDs(~[cellIDs.isdir]); % remove directories
    cellIDs = cellIDs(contains({cellIDs.name},'GR'));
    for jj = 1:size(cellIDs,1)
        spikes = load(cellIDs(jj).name, 'TS');
        ttl = load('shocKTTL.mat', 'shockTTL');
        g.cell_metrics.spikes.times{end+1} = spikes.TS/10000; % MClust TS must be divided by 10k
        g.cell_metrics.brainRegion{end+1} = 'BA';
        g.cell_metrics.general.shockTTL{end+1} = ttl.shockTTL;
        g.cell_metrics.cellID{end+1} = size(g.cell_metrics.cellID,2)+1;
    end
end

clearvars -except g

%% PSTH

psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'shockTTL', 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
psth_spx_zscore = zscore(psth_spx,0,2);
psth_spx_zscore = smoothdata(psth_spx_zscore,2,'sgolay', g.smoothvalue);
idx_LA = strcmp(g.cell_metrics.brainRegion, 'LA')';
idx_BA = strcmp(g.cell_metrics.brainRegion, 'BA')';
psth_spx_zscore_LA = psth_spx_zscore(idx_LA,:);
psth_spx_zscore_BA = psth_spx_zscore(idx_BA,:);
psth_spx_LA = psth_spx(idx_LA,:);
psth_spx_BA = psth_spx(idx_BA,:);

[results, clusters] = Gergo_psth_sorter(g,psth_spx_zscore);

[~, order_LA] = sortrows([results.respDir(idx_LA) clusters(idx_LA) results.onsetIdx(idx_LA)*g.bin_time*1000 results.offsetIdx(idx_LA)*g.bin_time*1000], [1 2 3 4], 'ascend');
[~, order_BA] = sortrows([results.respDir(idx_BA) clusters(idx_BA) results.onsetIdx(idx_BA)*g.bin_time*1000 results.offsetIdx(idx_BA)*g.bin_time*1000], [1 2 3 4], 'ascend');

psth_spx_zscore_LA = psth_spx_zscore_LA(order_LA,:);
psth_spx_zscore_BA = psth_spx_zscore_BA(order_BA,:);

g.clim = [min(psth_spx_zscore(:)) max(psth_spx_zscore(:))];


fig = figure('Position', [400, 100, 700, 800]);
t = tiledlayout(fig, 3,2,'TileSpacing', 'compact', 'Padding', 'none');

% Schematic figure

ax1 = nexttile(t,1,[1 2]);
[img, cmap] = imread([g.mainFolder '\mouse_drawing_anesthetized.png']);
if ~isempty(cmap)
    img = ind2rgb(img, cmap);  % Convert to RGB
end
imshow(img, 'Parent', ax1);
title(ax1, 'Experimental setup - anesthetized', 'FontSize', g.fontSize1)


% LA heatmap
nexttile 
imagesc(g.timeaxis,1:size(psth_spx_zscore_LA,1),psth_spx_zscore_LA);
clim(g.clim)
colormap(g.colors.Heatmap);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
%xlabel('Time (s)')
xticks([-1 0 1]);
yticks([5 10 size(psth_spx_zscore_LA,1)])
set(gca, 'FontSize', g.fontSize2);
title('LA neurons', 'FontSize', g.fontSize1)

ylabel('Neuron #')

% BA heatmapt
nexttile 
imagesc(g.timeaxis,1:size(psth_spx_zscore_BA,1),psth_spx_zscore_BA);
clim(g.clim)
colormap(g.colors.Heatmap);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
%xlabel('Time (s)')
xticks([-1 0 1]);
yticks([5 10 size(psth_spx_zscore_BA,1)])
set(gca, 'FontSize', g.fontSize2);
title('BA neurons', 'FontSize', g.fontSize1)

cb = colorbar('eastoutside', 'FontSize', g.fontSize2);

% LA firing rate histogram
nexttile
time = g.timeaxis(2:end);
meanData_LA = mean(psth_spx_LA, 1)/g.bin_time/size(psth_spx_LA,1);         % Mean across trials
stdData_PFC = std(psth_spx_LA, 0, 1)/g.bin_time/size(psth_spx_LA,1);        % Standard deviation across trials
upper_la = meanData_LA + stdData_PFC;
lower_la = meanData_LA - stdData_PFC;
fill([time, fliplr(time)], [upper_la, fliplr(lower_la)], [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);    % gray shading
hold on;
plot(time, meanData_LA, 'r', 'LineWidth', 2);   % black line
ylim([-10 20])
xticks([-1 0 1]);
xlabel('Time (s)');
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
ylabel('Firing rate (Hz)');
yticks([0 10 20])
set(gca, 'FontSize', g.fontSize2);
%title('PFC projecting group', 'FontSize', g.fontSize1);
box off;

% BA firing rate histogram
nexttile
time = g.timeaxis(2:end);
meanData_BA = mean(psth_spx_BA, 1)/g.bin_time/size(psth_spx_BA,1);         % Mean across trials
stdData_BA = std(psth_spx_BA, 0, 1)/g.bin_time/size(psth_spx_BA,1);        % Standard deviation across trials
upper_ba = meanData_BA + stdData_BA;
lower_ba = meanData_BA - stdData_BA;
fill([time, fliplr(time)], [upper_ba, fliplr(lower_ba)], [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);    % gray shading
hold on;
plot(time, meanData_BA, 'r', 'LineWidth', 2);   % black line
ylim([-10 20])
xticks([-1 0 1]);
xlabel('Time (s)');
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
yticks([0 10 20])
set(gca, 'FontSize', g.fontSize2);
%title('DMS projecting group', 'FontSize', g.fontSize1);
box off;
