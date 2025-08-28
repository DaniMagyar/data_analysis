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
g.optopre = 0.02;
g.optopost = 0.02;
g.optobin = 0.001;
g.optotimeaxis = -g.optopre:g.optobin:g.optopost;


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
psth_spx_zscore = zscore(psth_spx,0,2);
psth_spx_zscore = smoothdata(psth_spx_zscore,2,'sgolay', g.smoothvalue);
idx_PFC = strcmp(g.cell_metrics.projection, 'PFC')';
idx_DMS = strcmp(g.cell_metrics.projection, 'DMS')';
psth_spx_zscore_PFC = psth_spx_zscore(idx_PFC,:);
psth_spx_zscore_DMS = psth_spx_zscore(idx_DMS,:);
psth_spx_PFC = psth_spx(idx_PFC,:);
psth_spx_DMS = psth_spx(idx_DMS,:);

[results, clusters] = Gergo_psth_sorter(g,psth_spx_zscore);

[~, order_PFC] = sortrows([results.respDir(idx_PFC) clusters(idx_PFC) results.onsetIdx(idx_PFC)*g.bin_time*1000 results.offsetIdx(idx_PFC)*g.bin_time*1000], [1 2 3 4], 'ascend');
[~, order_DMS] = sortrows([results.respDir(idx_DMS) clusters(idx_DMS) results.onsetIdx(idx_DMS)*g.bin_time*1000 results.offsetIdx(idx_DMS)*g.bin_time*1000], [1 2 3 4], 'ascend');

psth_spx_zscore_PFC = psth_spx_zscore_PFC(order_PFC,:);
psth_spx_zscore_DMS = psth_spx_zscore_DMS(order_DMS,:);


g.clim = [min(psth_spx_zscore(:)) max(psth_spx_zscore(:))];


fig = figure('Position', [400, 100, 900, 1200]);
tiledlayout(fig, 4,2,'TileSpacing', 'compact', 'Padding', 'none')

% PFC heatmap
nexttile 
imagesc(g.timeaxis,1:size(psth_spx_zscore_PFC,1),psth_spx_zscore_PFC);
clim(g.clim)
colormap(g.colors.Heatmap);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
%xlabel('Time (s)')
xticks([-g.pre_time 0 g.post_time]);
yticks([5 10 size(psth_spx_zscore_PFC,1)])
set(gca, 'FontSize', g.fontSize2);
title('PFC projecting neurons - US reponse', 'FontSize', g.fontSize1)

ylabel('Neuron #')

% DMS heatmapt
nexttile 
imagesc(g.timeaxis,1:size(psth_spx_zscore_DMS,1),psth_spx_zscore_DMS);
clim(g.clim)
colormap(g.colors.Heatmap);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
%xlabel('Time (s)')
xticks([-g.pre_time 0 g.post_time]);
yticks([5 10 size(psth_spx_zscore_DMS,1)])
set(gca, 'FontSize', g.fontSize2);
title('DMS projecting neurons - US reponse', 'FontSize', g.fontSize1)

cb = colorbar('eastoutside', 'FontSize', g.fontSize2);

% PFC zscore histogram
nexttile
time = g.timeaxis(2:end);
meanData_PFC_z = mean(psth_spx_zscore_PFC, 1);         % Mean across trials
stdData_PFC_z = std(psth_spx_zscore_PFC, 0, 1);        % Standard deviation across trials
upper_pfc_z = meanData_PFC_z + stdData_PFC_z;
lower_pfc_z = meanData_PFC_z - stdData_PFC_z;
fill([time, fliplr(time)], [upper_pfc_z, fliplr(lower_pfc_z)], [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);    % gray shading
hold on;
plot(time, meanData_PFC_z, 'r', 'LineWidth', 2);   % black line
ylim([-1 2.5])
xticks([-g.pre_time 0 g.post_time]);
%xlabel('Time (s)');
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
ylabel('Z-scored firing rate');
yticks([-1 0 1 2])
set(gca, 'FontSize', g.fontSize2);
%title('PFC projecting group', 'FontSize', g.fontSize1);
box off;

% DMS zscore histogram
nexttile
time = g.timeaxis(2:end);
meanData_DMS_z = mean(psth_spx_zscore_DMS, 1);         % Mean across trials
stdData_DMS_z = std(psth_spx_zscore_DMS, 0, 1);        % Standard deviation across trials
upper_dms_z = meanData_DMS_z + stdData_DMS_z;
lower_dms_z = meanData_DMS_z - stdData_DMS_z;
fill([time, fliplr(time)], [upper_dms_z, fliplr(lower_dms_z)], [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);    % gray shading
hold on;
plot(time, meanData_DMS_z, 'r', 'LineWidth', 2);   % black line
ylim([-1 2.5])
xticks([-g.pre_time 0 g.post_time]);
%xlabel('Time (s)');
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
yticks([-1 0 1 2])
set(gca, 'FontSize', g.fontSize2);
%title('DMS projecting group', 'FontSize', g.fontSize1);
box off;

% PFC firing rate histogram
nexttile
time = g.timeaxis(2:end);
meanData_PFC = mean(psth_spx_PFC, 1)/g.bin_time/size(psth_spx_PFC,1);         % Mean across trials
stdData_PFC = std(psth_spx_PFC, 0, 1)/g.bin_time/size(psth_spx_PFC,1);        % Standard deviation across trials
upper_pfc = meanData_PFC + stdData_PFC;
lower_pfc = meanData_PFC - stdData_PFC;
fill([time, fliplr(time)], [upper_pfc, fliplr(lower_pfc)], [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);    % gray shading
hold on;
plot(time, meanData_PFC, 'r', 'LineWidth', 2);   % black line
ylim([-30 120])
xticks([-g.pre_time 0 g.post_time]);
%xlabel('Time (s)');
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
ylabel('Firing rate (Hz)');
yticks([ 0 50 100])
set(gca, 'FontSize', g.fontSize2);
%title('PFC projecting group', 'FontSize', g.fontSize1);
box off;

% DMS firing rate histogram
nexttile
time = g.timeaxis(2:end);
meanData_DMS = mean(psth_spx_DMS, 1)/g.bin_time/size(psth_spx_DMS,1);         % Mean across trials
stdData_DMS = std(psth_spx_DMS, 0, 1)/g.bin_time/size(psth_spx_DMS,1);        % Standard deviation across trials
upper_dms = meanData_DMS + stdData_DMS;
lower_dms = meanData_DMS - stdData_DMS;
fill([time, fliplr(time)], [upper_dms, fliplr(lower_dms)], [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);    % gray shading
hold on;
plot(time, meanData_DMS, 'r', 'LineWidth', 2);   % black line
ylim([-30 120])
xticks([-g.pre_time 0 g.post_time]);
%xlabel('Time (s)');
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
yticks([ 0 50 100])
set(gca, 'FontSize', g.fontSize2);
%title('DMS projecting group', 'FontSize', g.fontSize1);
box off;

clearvars -except g order_DMS order_PFC

%% Opto-tagging
%% PSTH
psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'optoTTL', 'pre_time', g.optopre, 'post_time', g.optopost, 'bin_time', g.optobin);
psth_spx_zscore = zscore(psth_spx,0,2);
%psth_spx_zscore = smoothdata(psth_spx_zscore,2,'sgolay', 3);
idx_PFC = strcmp(g.cell_metrics.projection, 'PFC')';
idx_DMS = strcmp(g.cell_metrics.projection, 'DMS')';
psth_spx_zscore_PFC = psth_spx_zscore(idx_PFC,:);
psth_spx_zscore_DMS = psth_spx_zscore(idx_DMS,:);
psth_spx_PFC = psth_spx(idx_PFC,:);
psth_spx_DMS = psth_spx(idx_DMS,:);
% 
% [results, clusters] = Gergo_psth_sorter(g,psth_spx_zscore);
% 
% [~, order_PFC] = sortrows([results.respDir(idx_PFC) clusters(idx_PFC) results.onsetIdx(idx_PFC)*g.bin_time*1000 results.offsetIdx(idx_PFC)*g.bin_time*1000], [1 2 3 4], 'ascend');
% [~, order_DMS] = sortrows([results.respDir(idx_DMS) clusters(idx_DMS) results.onsetIdx(idx_DMS)*g.bin_time*1000 results.offsetIdx(idx_DMS)*g.bin_time*1000], [1 2 3 4], 'ascend');
% 
psth_spx_zscore_PFC = psth_spx_zscore_PFC(order_PFC,:);
psth_spx_zscore_DMS = psth_spx_zscore_DMS(order_DMS,:);


g.clim = [min(psth_spx_zscore(:)) max(psth_spx_zscore(:))];


% PFC heatmap
nexttile 
imagesc(g.optotimeaxis,1:size(psth_spx_zscore_PFC,1),psth_spx_zscore_PFC);
clim(g.clim)
colormap(g.colors.Heatmap);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
xlabel('Time (s)')
xticks([-g.optopre 0 g.optopre]);
yticks([5 10 size(psth_spx_zscore_PFC,1)])
set(gca, 'FontSize', g.fontSize2);
title('Optotagging', 'FontSize', g.fontSize1)

ylabel('Neuron #')

% DMS heatmapt
nexttile 
imagesc(g.optotimeaxis,1:size(psth_spx_zscore_DMS,1),psth_spx_zscore_DMS);
clim(g.clim)
colormap(g.colors.Heatmap);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
xlabel('Time (s)')
xticks([-g.optopre 0 g.optopre]);
yticks([5 10 size(psth_spx_zscore_DMS,1)])
set(gca, 'FontSize', g.fontSize2);
title('Optotagging', 'FontSize', g.fontSize1)

cb = colorbar('eastoutside', 'FontSize', g.fontSize2);
