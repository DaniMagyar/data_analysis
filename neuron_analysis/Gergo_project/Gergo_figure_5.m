%% function Gergo_figure_5
% Aneshetized recordings, only responsive neurons

clear all
recordings = {'MD313_kilosort'};
g.mainFolder = 'C:\Users\dmagyar\Desktop\Gergo';
g.cell_metrics = BAfc_load_neurons('mainFolder', g.mainFolder, 'recordings', recordings, 'ttl', {'triptest_shocks_only'});
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);

clearvars -except g 

g.colors = BAfc_colors;
g.fontSize1 = 13;
g.fontSize2 = 12;
g.pre_time = 0.2;
g.post_time = 0.2;
g.bin_time = 0.01;

%% Initiate figure
fig = figure('Position', [400, 100, 800, 1200]);
t = tiledlayout(fig,12,8,'TileSpacing', 'tight', 'Padding', 'none');

% %% (1,1) - Schematic figure
% ax1 = nexttile(t,1,[4 4]);
% [img, cmap] = imread([g.mainFolder '\probe_shank_neurons.png']);
% if ~isempty(cmap)
%     img = ind2rgb(img, cmap);  % Convert to RGB
% end
% imshow(img, 'Parent', ax1);
% title(ax1, 'Experimental setup', 'FontSize', g.fontSize1)
% clearvars -except t g % clear variables
% 
% cellID_1 = intersect(find(strcmp(g.cell_metrics.animal, 'MD313')), find(g.cell_metrics.cluID == 530)); 
% cellID_2 = intersect(find(strcmp(g.cell_metrics.animal, 'MD313')), find(g.cell_metrics.cluID == 468)); 
% 
% spike_times{1} = g.cell_metrics.spikes.times{cellID_1}; 
% spike_times{2} = g.cell_metrics.spikes.times{cellID_2}; 
% 
% stimulus_times{1} = g.cell_metrics.general.triptest_shocks_only{cellID_1};
% stimulus_times{2} = g.cell_metrics.general.triptest_shocks_only{cellID_2};
% 
% % Plot raster
% time_bins = -g.pre_time:g.bin_time:g.post_time;    
% spike_counts = zeros(1, length(time_bins) - 1);   
% rasterloc = [5 21];
% rastertitles = {'Example LA neurons (Cell 1)','Example BA neurons (Cell 2)'};
% for ii = 1:2
%     ax = nexttile(t,rasterloc(ii),[2 4]);
%     hold on;    
%     for trial = 1:length(stimulus_times{ii})
%         aligned_spikes = spike_times{ii} - stimulus_times{ii}(trial);
%         valid_spikes = aligned_spikes(aligned_spikes >= -g.pre_time & aligned_spikes <= g.post_time);
%         scatter(valid_spikes, trial * ones(size(valid_spikes)), 10, 'k', 'filled');
%         spike_counts = spike_counts + histcounts(valid_spikes, time_bins);
%     end
%     xline(0, 'r--', 'LineWidth', 2)
%     xlim([-g.pre_time g.post_time])
%     % xlabel('Time (s)')
%     xticks([-g.pre_time 0 g.post_time]);
%     title(rastertitles{ii})
%     ylabel('Trials')
%     set(gca, 'FontSize', g.fontSize2)
% end


%% function Gergo_figure_2

% Anesthetized data

clearvars -except t

g.mainFolder = 'C:\Users\dmagyar\Desktop\Gergo';
g.subfolders.LA = 'C:\Users\dmagyar\Desktop\Gergo\LA_urethane Gergo\LA_urethane';
folderList_LA = dir(g.subfolders.LA);
folderList_LA = folderList_LA([folderList_LA.isdir] & ~ismember({folderList_LA.name}, {'.', '..'}));
g.fullPaths_LA = fullfile(g.subfolders.LA, {folderList_LA.name});

g.subfolders.BA = 'C:\Users\dmagyar\Desktop\Gergo\BA_urethane Gergo\BA_urethane';
folderList_BA = dir(g.subfolders.BA);
folderList_BA = folderList_BA([folderList_BA.isdir] & ~ismember({folderList_BA.name}, {'.', '..'}));
g.fullPaths_BA = fullfile(g.subfolders.BA, {folderList_BA.name});


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

clearvars -except g t


g.pre_time = 0.2;
g.post_time = 0.2;
g.bin_time = 0.005;
g.smoothvalue = 7;
g.timeaxis = -g.pre_time:g.bin_time:g.post_time;
g.colors = BAfc_colors;
g.test_time = 0.2;
g.exctestvalue = 1; % if the mean change is negative, but there is a peak where zscore is above this, than counted as excited
g.inhtestvalue = 1;
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
g.fontSize2 = 12;
g.fontSize1 = 15;
g.xlinewidth = 2;

% PSTH

psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'shockTTL', 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
psth_spx(:, g.pre_time/g.bin_time+1:g.pre_time/g.bin_time+2) = 0;
disp('binsize must be 5ms, manually removed spikes from psth_spx')
psth_spx_zscore = zscore(psth_spx,0,2);
psth_spx_zscore = smoothdata(psth_spx_zscore,2,'sgolay', g.smoothvalue);

n_stim = cellfun(@numel,g.cell_metrics.general.shockTTL)';
n_post_spx = sum(psth_spx(:, (g.pre_time/g.bin_time)+1:end),2);

%g.cell_metrics.brainRegion((n_post_spx./n_stim)<0.25) = {'low'};
%g.cell_metrics.brainRegion(n_post_spx<5) = {'low'};

idx_LA = strcmp(g.cell_metrics.brainRegion, 'LA')';
idx_BA = strcmp(g.cell_metrics.brainRegion, 'BA')';
psth_spx_zscore_LA = psth_spx_zscore(idx_LA,:);
psth_spx_zscore_BA = psth_spx_zscore(idx_BA,:);

[results, clusters] = Gergo_psth_sorter(g,psth_spx_zscore);

% [~, order_LA] = sortrows([results.respDir(idx_LA) clusters(idx_LA) results.onsetIdx(idx_LA)*g.bin_time*1000 results.offsetIdx(idx_LA)*g.bin_time*1000], [1 2 3 4], 'ascend');
% [~, order_BA] = sortrows([results.respDir(idx_BA) clusters(idx_BA) results.onsetIdx(idx_BA)*g.bin_time*1000 results.offsetIdx(idx_BA)*g.bin_time*1000], [1 2 3 4], 'ascend');

[~, order_LA] = sortrows([results.respDir(idx_LA) results.onsetIdx(idx_LA)*g.bin_time*1000], [1 2], 'ascend');
[~, order_BA] = sortrows([results.respDir(idx_BA) results.onsetIdx(idx_BA)*g.bin_time*1000], [1 2], 'ascend');

n_LA_excited = sum(clusters(idx_LA) == 1) + sum(clusters(idx_LA) == 2);
n_BA_excited = sum(clusters(idx_BA) == 1) + sum(clusters(idx_BA) == 2);

psth_spx_zscore_LA = psth_spx_zscore_LA(order_LA(1:n_LA_excited),:);
psth_spx_zscore_BA = psth_spx_zscore_BA(order_BA(1:n_BA_excited),:);

g.clim = [min(psth_spx_zscore(:)) max(psth_spx_zscore(:))*2];

% LA heatmap
nexttile(t,33,[4 4]);
imagesc(g.timeaxis,1:size(psth_spx_zscore_LA,1),psth_spx_zscore_LA);
clim(g.clim)
colormap(g.colors.Heatmap);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
xlabel('Time (s)')
xticks([-g.pre_time 0 g.post_time]);
yticks([1 size(psth_spx_zscore_LA,1)])
set(gca, 'FontSize', g.fontSize2);
title('LA neurons', 'FontSize', g.fontSize1)

ylabel('Neuron #')

% BA heatmapt
nexttile(t,37,[4 4]);
imagesc(g.timeaxis,1:size(psth_spx_zscore_BA,1),psth_spx_zscore_BA);
clim(g.clim)
colormap(g.colors.Heatmap);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
xlabel('Time (s)')
xticks([-g.pre_time 0 g.post_time]);
yticks([1 size(psth_spx_zscore_BA,1)])
set(gca, 'FontSize', g.fontSize2);
title('BA neurons', 'FontSize', g.fontSize1)

cb = colorbar('eastoutside', 'FontSize', g.fontSize2);




barbinsize = 0.002;
bartimeaxis = -g.pre_time:barbinsize:g.post_time;
psth_spx_bar =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'shockTTL', 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', barbinsize);
psth_spx_bar(:, g.pre_time/barbinsize+1:g.pre_time/barbinsize+5) = 0;
% psth_spx_bar = smoothdata(psth_spx_bar,2,'movmean', 5);
psth_spx_LA = psth_spx_bar(idx_LA,:);
psth_spx_BA = psth_spx_bar(idx_BA,:);


% LA firing rate histogram
nexttile(t,65,[4 4]);
psth_spx_excited_LA = psth_spx_LA(order_LA(1:n_LA_excited),:);
time = bartimeaxis(2:end);
allspikes_LA = sum(psth_spx_excited_LA,1);       % Mean across trials
bar(time, allspikes_LA)
xlabel('Time (s)')
xticks([-g.pre_time 0 g.post_time]);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
set(gca, 'FontSize', g.fontSize2);
%title('PFC projecting group', 'FontSize', g.fontSize1);
box off;

% BA firing rate histogram
nexttile(t,69,[4 4]);
psth_spx_excited_BA = psth_spx_BA(order_BA(1:n_BA_excited),:);
time = bartimeaxis(2:end);
allspikes_BA = sum(psth_spx_excited_BA,1);       % Mean across trials
bar(time, allspikes_BA)
xlabel('Time (s)')
xticks([-g.pre_time 0 g.post_time]);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
set(gca, 'FontSize', g.fontSize2);
%title('PFC projecting group', 'FontSize', g.fontSize1);
box off;



LA_neurons_sorted_idx = find(idx_LA);
LA_neurons_sorted_idx = LA_neurons_sorted_idx(order_LA(1:n_LA_excited));

BA_neurons_sorted_idx = find(idx_BA);
BA_neurons_sorted_idx = BA_neurons_sorted_idx(order_BA(1:n_BA_excited));