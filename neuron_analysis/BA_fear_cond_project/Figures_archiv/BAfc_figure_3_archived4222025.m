% function BAfc_figure_3
% % Plotting region specific responses on heatmaps
% % 
% recordings = {...
%     'MD292_002_kilosort',...
%     'MD293_kilosort',...
%     'MD294_kilosort',...
%     'MD295_kilosort',...
%     'MD296_kilosort',...
%     'MD297_kilosort'};
% g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'});

clearvars -except g
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 12;
g.fontSize2 = 12;
g.fontsize3 = 10;
g.pre_time = 5;
g.post_time = 1;
g.bin_time = 0.001;
g.test_win = [0.012 0.1];
g.test_time = 0.5;
g.smoothvalue = 5;
g.testvalue = 1.5;


%% Initiate figure
fig = figure('Position', [400, 100, 1200, 1200]);
t = tiledlayout(fig,6,3,'TileSpacing', 'compact', 'Padding', 'none');
%% (1:4,1:4) - Example neurons
cluIDs = [273 106 315 465];
animals = {'MD297', 'MD295', 'MD294', 'MD292'}; % LA PN, LA IN, BA PN, BA IN
ttl = {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'};
for ii = 1:4
    cellID = intersect(find(strcmp(g.cell_metrics.animal, animals{ii})), find(g.cell_metrics.cluID == cluIDs(ii)));
    for jj = 1:3
        ax = nexttile;
        plotPSTHLocal(g.cell_metrics.spikes.times{cellID}, g.cell_metrics.general.(ttl{jj}){cellID}, [-0.5 1], 0.01)
    end
end

clearvars -except g t

%% All responsive neurons comparison
ttl = {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'};
for ii = 1:3
    [psth_spx_all.(ttl{ii}), results_all.(ttl{ii}), zResp_median.(ttl{ii})] =  BAfc_find_response('cell_metrics', g.cell_metrics, 'ttl', ttl{ii}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time, 'test_win', g.test_win, 'smoothvalue', g.smoothvalue, 'test_value',g.testvalue);
    latencies.(ttl{ii}) =  calcLatencyLocal(psth_spx_all.(ttl{ii}),g);
end

idx_LA_PN = find((strcmp(results_all.(ttl{1}), 'exc') | strcmp(results_all.(ttl{2}), 'exc') | strcmp(results_all.(ttl{3}), 'exc')) & ...
    (strcmp(g.cell_metrics.brainRegion, 'LA')' & strcmp(g.cell_metrics.putativeCellType, 'PN')'));
idx_LA_IN = find((strcmp(results_all.(ttl{1}), 'exc') | strcmp(results_all.(ttl{2}), 'exc') | strcmp(results_all.(ttl{3}), 'exc')) & ...
    (strcmp(g.cell_metrics.brainRegion, 'LA')' & strcmp(g.cell_metrics.putativeCellType, 'IN')'));
idx_BA_PN = find((strcmp(results_all.(ttl{1}), 'exc') | strcmp(results_all.(ttl{2}), 'exc') | strcmp(results_all.(ttl{3}), 'exc')) & ...
    (strcmp(g.cell_metrics.brainRegion, 'BA')' & strcmp(g.cell_metrics.putativeCellType, 'PN')'));
idx_BA_IN = find((strcmp(results_all.(ttl{1}), 'exc') | strcmp(results_all.(ttl{2}), 'exc') | strcmp(results_all.(ttl{3}), 'exc')) & ...
    (strcmp(g.cell_metrics.brainRegion, 'BA')' & strcmp(g.cell_metrics.putativeCellType, 'IN')'));


% Plot LA PN and IN 
ax = nexttile([2 1]);
hold on;
for ii = idx_LA_PN
    plot([1, 2, 3], [zResp_median.triptest_sound_only(ii), zResp_median.triptest_shocks_only(ii), zResp_median.triptest_both(ii)], 'Color', g.colors.PN_primary); % Line connecting paired points
    plot(1, zResp_median.triptest_sound_only(ii), 'x', 'Color', g.colors.PN_primary); % Data point from vector1
    plot(2, zResp_median.triptest_shocks_only(ii), 'x', 'Color', g.colors.PN_primary); % Data point from vector2
    plot(3, zResp_median.triptest_both(ii), 'x', 'Color', g.colors.PN_primary); % Data point from vector2
end
for ii = idx_LA_IN
    plot([1, 2, 3], [zResp_median.triptest_sound_only(ii), zResp_median.triptest_shocks_only(ii), zResp_median.triptest_both(ii)], 'Color', g.colors.IN_primary); % Line connecting paired points
    plot(1, zResp_median.triptest_sound_only(ii), 'x', 'Color', g.colors.IN_primary); % Data point from vector1
    plot(2, zResp_median.triptest_shocks_only(ii), 'x', 'Color', g.colors.IN_primary); % Data point from vector2
    plot(3, zResp_median.triptest_both(ii), 'x', 'Color', g.colors.IN_primary); % Data point from vector2
end
% Plot BA PN and IN 
ax = nexttile([2 1]);
hold on;
for ii = idx_BA_PN
    plot([1, 2, 3], [zResp_median.triptest_sound_only(ii), zResp_median.triptest_shocks_only(ii), zResp_median.triptest_both(ii)], 'Color', g.colors.PN_primary); % Line connecting paired points
    plot(1, zResp_median.triptest_sound_only(ii), 'x', 'Color', g.colors.PN_primary); % Data point from vector1
    plot(2, zResp_median.triptest_shocks_only(ii), 'x', 'Color', g.colors.PN_primary); % Data point from vector2
    plot(3, zResp_median.triptest_both(ii), 'x', 'Color', g.colors.PN_primary); % Data point from vector2
end
for ii = idx_BA_IN
    plot([1, 2, 3], [zResp_median.triptest_sound_only(ii), zResp_median.triptest_shocks_only(ii), zResp_median.triptest_both(ii)], 'Color', g.colors.IN_primary); % Line connecting paired points
    plot(1, zResp_median.triptest_sound_only(ii), 'x', 'Color', g.colors.IN_primary); % Data point from vector1
    plot(2, zResp_median.triptest_shocks_only(ii), 'x', 'Color', g.colors.IN_primary); % Data point from vector2
    plot(3, zResp_median.triptest_both(ii), 'x', 'Color', g.colors.IN_primary); % Data point from vector2
end
ax = nexttile([2 1]);

data{1} = latencies.triptest_sound_only(find(strcmp(results_all.triptest_sound_only, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'LA')' & strcmp(g.cell_metrics.putativeCellType, 'PN')'))); 
data{2} = latencies.triptest_sound_only(find(strcmp(results_all.triptest_sound_only, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'LA')' & strcmp(g.cell_metrics.putativeCellType, 'IN')')));
data{3} = latencies.triptest_sound_only(find(strcmp(results_all.triptest_sound_only, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'BA')' & strcmp(g.cell_metrics.putativeCellType, 'PN')')));
data{4} = latencies.triptest_sound_only(find(strcmp(results_all.triptest_sound_only, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'BA')' & strcmp(g.cell_metrics.putativeCellType, 'IN')')));

data{5} = latencies.triptest_shocks_only(find(strcmp(results_all.triptest_shocks_only, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'LA')' & strcmp(g.cell_metrics.putativeCellType, 'PN')'))); 
data{6} = latencies.triptest_shocks_only(find(strcmp(results_all.triptest_shocks_only, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'LA')' & strcmp(g.cell_metrics.putativeCellType, 'IN')'))); 
data{7} = latencies.triptest_shocks_only(find(strcmp(results_all.triptest_shocks_only, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'BA')' & strcmp(g.cell_metrics.putativeCellType, 'PN')'))); 
data{8} = latencies.triptest_shocks_only(find(strcmp(results_all.triptest_shocks_only, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'BA')' & strcmp(g.cell_metrics.putativeCellType, 'IN')'))); 

data{9} = latencies.triptest_both(find(strcmp(results_all.triptest_both, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'LA')' & strcmp(g.cell_metrics.putativeCellType, 'PN')'))); 
data{10} = latencies.triptest_both(find(strcmp(results_all.triptest_both, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'LA')' & strcmp(g.cell_metrics.putativeCellType, 'IN')'))); 
data{11} = latencies.triptest_both(find(strcmp(results_all.triptest_both, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'BA')' & strcmp(g.cell_metrics.putativeCellType, 'PN')'))); 
data{12} = latencies.triptest_both(find(strcmp(results_all.triptest_both, 'exc') & (strcmp(g.cell_metrics.brainRegion, 'BA')' & strcmp(g.cell_metrics.putativeCellType, 'IN')'))); 

allData = [];
group = [];
for i = 1:length(data)
    currentData = data{i};
    allData = [allData; currentData(:)];
    group = [group; i * ones(length(currentData), 1)];
end

boxplot(allData, group)
xlabel('Group')
ylabel('Values')
title('Boxplot of Data in First Row of Cell Array')



function plotPSTHLocal(spikeTimes, stimTimes, window, binSize)
    edges = window(1):binSize:window(2);
    counts = zeros(1, length(edges)-1);
    for trial = 1:length(stimTimes)
        alignedSpikes = spikeTimes - stimTimes(trial);
        trialSpikes = alignedSpikes(alignedSpikes >= window(1) & alignedSpikes <= window(2));
        counts = counts + histcounts(trialSpikes, edges);
    end
    % Convert counts to firing rate (spikes per second)
    firingRate = counts / length(stimTimes) / binSize;
    binCenters = edges(1:end-1) + binSize/2;
    bar(binCenters, firingRate, 'k', 'BarWidth', 1);
    xlabel('Time (s)');
    ylabel('Firing rate (Hz)');
    title('Peri-Stimulus Time Histogram');
end


function [latency] = calcLatencyLocal(psth_spx,g)
    roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
    psth_zscore = zeros(size(psth_spx));
    latency = zeros(size(psth_spx,1),1);
    for ii = 1:size(psth_spx,1)
        psth_zscore(ii,:) = (psth_spx(ii,:) - mean(psth_spx(ii,1:g.pre_time/g.bin_time)))/std(psth_spx(ii,1:g.pre_time/g.bin_time));
        psth_zscore(ii,:) = smoothdata(psth_zscore(ii,:), 'movmean', g.smoothvalue);
        datasegment = (psth_zscore(ii,roi));
        datasegment = round(datasegment*10000)/10000;
        if all(datasegment == datasegment(1)) || all(datasegment<0)
            latency(ii,1) = NaN;
        else
            resp_peak = max(datasegment);
            response_start = find(datasegment>=resp_peak/2,1,'first'); % begining of the response
            latency(ii,1) = round(response_start*g.bin_time*1000);
        end
    end
end