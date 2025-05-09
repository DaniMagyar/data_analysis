% function BAfc_figure_3
% % Plotting region specific responses on heatmaps
clear all
recordings = {...
    'MD292_002_kilosort',...
    'MD293_kilosort',...
    'MD294_kilosort',...
    'MD295_kilosort',...
    'MD296_kilosort',...
    'MD297_kilosort'};
g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'});

clearvars -except g
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 12;
g.fontSize2 = 12;
g.fontsize3 = 10;
g.pre_time = 5;
g.post_time = 5;
g.bin_time = 0.001;
g.test_time = 0.05;
g.smoothvalue = 11;
g.testvalue = 3;
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;


%% Initiate figure
fig = figure('Position', [400, 100, 1400, 1200]);
t = tiledlayout(fig,3,3,'TileSpacing', 'compact', 'Padding', 'none');
%% (1:4,1:4) - Example neurons
cluIDs = [259 329];% LA PN, LA PN
animals = {'MD295', 'MD294'}; % LA PN, LA PN

ttl = {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'};
for ii = 1:2
    cellID = intersect(find(strcmp(g.cell_metrics.animal, animals{ii})), find(g.cell_metrics.cluID == cluIDs(ii)));
    
    titles = {'CS', 'US', 'CS + US'};
    for jj = 1:3
        ax = nexttile;
        plotPSTHLocal(g.cell_metrics.spikes.times{cellID}, g.cell_metrics.general.(ttl{jj}){cellID}, [-0.5 1], 0.01);
        xlabel('Time (s)');
        if ii == 1 && jj == 1
            ylabel({'Example neuron 1', 'Firing Rate (Hz)'})
        elseif ii == 2 && jj == 1
            ylabel({'Example neuron 2', 'Firing Rate (Hz)'})
        end
        yticks([0 50 100])
        xticks([-0.5 0 0.5 1])
        set(gca, 'FontSize', g.fontSize1);
        box off
        if ii == 1 
            title(titles{jj}, 'FontSize', 20)
        end
    end
end

clearvars -except g t ttl

%% All responsive neurons comparison
for hmp = 1:3
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    psth_spx(:,g.pre_time/g.bin_time+1:(g.pre_time+0.012)/g.bin_time) = NaN;
    psth_spx = smoothdata(psth_spx,2,'sgolay', g.smoothvalue);
    psth_spx(:,g.pre_time/g.bin_time+1:(g.pre_time+0.012)/g.bin_time) = NaN;
    %psth_spx_zscore = zscore(psth_spx,0,2);       
    mu = mean(psth_spx, 2, 'omitnan');
    sigma = std(psth_spx, 0, 2, 'omitnan');
    psth_spx_zscore = (psth_spx - mu)./sigma;
    
    idx_PN = strcmp(g.cell_metrics.brainRegion,'LA') & strcmp(g.cell_metrics.putativeCellType,'PN');

    for ii = 1:size(psth_spx_zscore,1)
        datasegment = psth_spx_zscore(ii,g.roi);
        if all(datasegment == datasegment(1)) ||  max(datasegment)<g.testvalue
            results.onsetIdx(ii,1) = NaN;
            results.offsetIdx(ii,1) = NaN; 
        else
            results.onsetIdx(ii,1) = find(datasegment>=g.testvalue,1,'first'); % begining of the response
            results.offsetIdx(ii,1) = find(datasegment>=g.testvalue,1,'last'); % end of response
            responseMean = mean(datasegment(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
            k = 0;
            while k == 0
                if responseMean < g.testvalue/2
                    datasegment(results.offsetIdx(ii,1):end) = NaN;
                    results.offsetIdx(ii,1) = find(datasegment>=g.testvalue,1,'last'); % end of response
                    responseMean = mean(datasegment(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
                else 
                    k = 1;
                end
            end
        end
    end

    data.(ttl{hmp}).results = results;
    data.(ttl{hmp}).num_spx50 = sum(psth_spx(:,(g.pre_time+0.012)/g.bin_time:(g.pre_time+0.05)/g.bin_time),2);
    data.(ttl{hmp}).zscore = mean(psth_spx_zscore(:,(g.pre_time+0.012)/g.bin_time:(g.pre_time+0.05)/g.bin_time),2);
end
latencylimit = 40;
nspikelimit = 10;
excited = ...
    (data.triptest_sound_only.results.onsetIdx <= latencylimit & data.triptest_sound_only.num_spx50 >=nspikelimit) | ...
    (data.triptest_shocks_only.results.onsetIdx <= latencylimit & data.triptest_shocks_only.num_spx50 >=nspikelimit) | ...
    (data.triptest_both.results.onsetIdx <= latencylimit & data.triptest_both.num_spx50 >=nspikelimit);


idx_LA_PN = find(excited & strcmp(g.cell_metrics.brainRegion,'LA')' & strcmp(g.cell_metrics.putativeCellType,'PN')');
idx_LA_IN = find(excited & strcmp(g.cell_metrics.brainRegion,'LA')' & strcmp(g.cell_metrics.putativeCellType,'IN')');

% Plot LA PN and IN 

data1{1} = data.triptest_sound_only.zscore(idx_LA_PN);
data1{2} = data.triptest_shocks_only.zscore(idx_LA_PN);
data1{3} = data.triptest_both.zscore(idx_LA_PN);

data1{4} = data.triptest_sound_only.zscore(idx_LA_IN);
data1{5} = data.triptest_shocks_only.zscore(idx_LA_IN);
data1{6} = data.triptest_both.zscore(idx_LA_IN);

% Create boxplot
allData1 = [];
group1 = [];
for i = 1:length(data1)
    currentData = data1{i};
    allData1 = [allData1; currentData(:)];
    group1 = [group1; i * ones(length(currentData), 1)];
end

nexttile(t)
h1 = boxplot(allData1, group1);
ylim([-2 15])
for i = 1:6
    if i <= 3
        col =  g.colors.PN_primary; % blue
    else
        col =  g.colors.IN_primary; % red
    end
    % Set color for all elements of the i-th box
    for j = 1:7
        set(h1(j,i), 'Color', col);
    end
end
hold on

% Plot individual points and connect them
jitterAmount = 0.2;

for ii = 1:2  % 1: LA_PN (groups 1–3), 2: LA_IN (groups 4–6)
    baseIdx = (ii-1)*3;
    nPoints = min(cellfun(@length, data1(baseIdx + (1:3)))); % ensure equal length

    for i = 1:nPoints
        x = baseIdx + (1:3);
        y = [data1{baseIdx+1}(i), data1{baseIdx+2}(i), data1{baseIdx+3}(i)];
        xJitter = x + (rand(1,3)-0.5) * jitterAmount;
        if ii ==1 
            plot(xJitter, y, '-', 'Color', [g.colors.PN_primary, 0.5], 'LineWidth', 1)
            scatter(xJitter, y, 25, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', g.colors.PN_primary)
        else
            plot(xJitter, y, '-', 'Color', [g.colors.IN_primary, 0.5], 'LineWidth', 1)
            scatter(xJitter, y, 25, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', g.colors.IN_primary)
        end
    end
end
ylim([-2 15])
ylabel('Z-score')
title('Response amplitude of individual neurons')
xticks([1 2 3 4 5 6])                          % set the tick locations
xticklabels({'CS', 'US', 'CS+US', 'CS', 'US', 'CS+US'})  % custom labels
c1 = scatter(NaN, NaN, 100, g.colors.PN_primary, 'filled', 'MarkerFaceAlpha', 0.5);
set(c1, 'MarkerFaceAlpha', 0.5);
hold on;
c2 = scatter(NaN, NaN, 100, g.colors.IN_primary, 'filled', 'MarkerFaceAlpha', 0.5);
set(c2, 'MarkerFaceAlpha', 0.5);
% Add legend
legend([c1 c2], {'PN', 'IN'}, 'Location', 'northeast');
set(gca, 'FontSize', g.fontSize1);



data2{1} = data.triptest_sound_only.results.onsetIdx(idx_LA_PN);
data2{2} = data.triptest_shocks_only.results.onsetIdx(idx_LA_PN);
data2{3} = data.triptest_both.results.onsetIdx(idx_LA_PN);

data2{4} = data.triptest_sound_only.results.onsetIdx(idx_LA_IN);
data2{5} = data.triptest_shocks_only.results.onsetIdx(idx_LA_IN);
data2{6} = data.triptest_both.results.onsetIdx(idx_LA_IN);
allData2 = [];
group2 = [];
for i = 1:length(data2)
    currentData = data2{i};
    allData2 = [allData2; currentData(:)];
    group2 = [group2; i * ones(length(currentData), 1)];
end

nexttile
h2 = boxplot(allData2, group2, 'Symbol', '');
for i = 1:6
    if i <= 3
        col =  g.colors.PN_primary; % blue
    else
        col =  g.colors.IN_primary; % red
    end
    % Set color for all elements of the i-th box
    for j = 1:7
        set(h2(j,i), 'Color', col)
    end
end


hold on

% Plot individual data points with jitter
for i = 1:length(data2)
    x = i + 0.2*(rand(size(data2{i})) - 0.5); % jitter around group index
    if i == 1 || i ==2 || i == 3
        scatter(x, data2{i}, 25, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', g.colors.PN_primary)
    else
        scatter(x, data2{i}, 25, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', g.colors.IN_primary)
    end
end
ylim([0 50])
ylabel('Latency (ms)')
title('Response latency of individual neurons')
xticks([1 2 3 4 5 6])                          % set the tick locations
xticklabels({'CS', 'US', 'CS+US', 'CS', 'US', 'CS+US'})  % custom labels
c1 = scatter(NaN, NaN, 100, g.colors.PN_primary, 'filled', 'MarkerFaceAlpha', 0.5);
set(c1, 'MarkerFaceAlpha', 0.5);
hold on;
c2 = scatter(NaN, NaN, 100, g.colors.IN_primary, 'filled', 'MarkerFaceAlpha', 0.5);
set(c2, 'MarkerFaceAlpha', 0.5);
% Add legend
legend([c1 c2], {'PN', 'IN'}, 'Location', 'northeast');
set(gca, 'FontSize', g.fontSize1);



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
    ylim([0 100])
end