% function BAfc_figure_2_v3
% % Plotting region specific responses on heatmaps, proper clustering
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
g.test_time = 0.2; % 1 sec is good, beacuse captures better than 0.5
g.testvalue = 3; % if the mean change is negative, but there is a peak where zscore is above this, than counted as excited
g.bin_time = 0.005; % 0.01 volt sokat 5 os smoothal
g.smoothvalue = 5;
g.plotwin = [0.5 1];
g.spxwin = [0.02 0.1];
g.spxbin = 0.002;
g.timeaxis = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.clim = [-5.5 13];
g.clustnum = 4;
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
g.roi_pca = round((g.pre_time-0.5)/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);

%% Initiate figure
fig = figure('Position', [400, 100, 1400, 1200]);
t = tiledlayout(fig,6,4,'TileSpacing', 'compact', 'Padding', 'none');

%% (1,1:4) - Heatmaps

ttl = {'triptest_sound_only', 'triptest_shocks_only','triptest_both'};
hmptitles = {'CS', 'US', 'CS + US'};
PSTHall = [];
for hmp = 1:size(ttl,2)
    psth_spx_og =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    psth_spx = zscore(psth_spx_og,0,2);   
    psth_spx = smoothdata(psth_spx,2,'sgolay', g.smoothvalue);
    idx_PN = strcmp(g.cell_metrics.brainRegion,'LA') & strcmp(g.cell_metrics.putativeCellType,'PN');
    psth_spx_og_PN = psth_spx_og(idx_PN,:);
    psth_spx_PN = psth_spx(idx_PN,:);
    hmpdata.(['psth_' num2str(hmp)]) = psth_spx_PN(:,g.roi_pca);
    PSTHall = [PSTHall hmpdata.(['psth_' num2str(hmp)])]; 
end

[~,PCA1]  = pca(PSTHall); % running pricipal component analysis
PCA2      = PCA1(:,1:2); % extracting the first PCA_num pricinpal components 
inNAN = nanmedian(PCA2,1); % to avoid NaN values, set them to the median
for jj = 1:size(PSTHall,1)
if isnan(PCA2(jj,1)) == 1
    PCA2(jj,1:3) = inNAN;
end
end
Dend      = linkage(PCA2,'complete','mahalanobis'); % calculating the dendrogram
Clusters  = cluster(Dend,'maxclust',g.clustnum); % clustering based on the tree
D         = pdist(PCA2); %euclidean distrance between point in the artifical space
leafOrder = optimalleaforder(Dend,D)'; % optimal order for plotting    
PCA        = PCA2;
Dendrogram = Dend;
leafOrderfl = leafOrder;
%leafOrderfl([54 57]) = leafOrderfl([57 54]);
mploc = [1 2 3];
for hmp = 1:size(ttl,2)
    ax = nexttile(t,mploc(hmp),[2 1]);
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    psth_spx = zscore(psth_spx,0,2);   
    psth_spx = smoothdata(psth_spx,2,'sgolay', g.smoothvalue);
    idx_PN = strcmp(g.cell_metrics.brainRegion,'LA') & strcmp(g.cell_metrics.putativeCellType,'PN');
    psth_spx_PN = psth_spx(idx_PN,:);
    psth_spx = [psth_spx_PN(leafOrderfl,:)];
    matrix = psth_spx(:,(g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    imagesc(g.timeaxis,1:size(matrix,1),matrix);
    title(hmptitles{hmp}, 'FontSize', 20)
    clim(g.clim)
    disp([-max(max(psth_spx,[],1))/2.5 max(max(psth_spx,[],1))])   
    colormap(g.colors.Heatmap); 
    yline(sum(idx_PN)+0.5, 'k', 'LineWidth', 1);
    xlabel('Time (s)', 'FontSize', g.fontSize2)
    ylabel('Cell number', 'FontSize', g.fontSize2)
    if hmp == 1
        cb = colorbar('westoutside', 'FontSize',15);
    end
    hold on
    k1 = 0;
    for ii = 1:g.clustnum
        k2 = sum(Clusters == ii);
        yline(k1+k2+0.5, 'k', 'LineWidth', 1);
        k1 = k1 + k2;
    end
    hold off
end
nexttile(t,4,[2 1])
cutoff = Dend(end-g.clustnum+2,3); % cutting the tree to get 'clustnum' clusters
h = dendrogram(Dend,0,'Orientation','right','reorder',flipud(leafOrder),...
    'ColorThreshold',cutoff); %plotting dendrogram
axis off


%% (1:4,1:4) - Example neurons
%  cellID 193: CS+, US-, CSUS+ CluID: 298 Animal: MD296
%  cellID 196: CS+, US-, CSUS- CluID: 330 Animal: MD296
%  cellID 241: CS-, US+, CSUS+ CluID: 338 Animal: MD297
%  cellID 218: CS-, US-, CSUS- Cluid: 238 Anumal: MD297

cluIDs = [298 330 338 238];% LA PN, LA PN
animals = {'MD296', 'MD296', 'MD297', 'MD297'}; % LA PN, LA PN

ttl = {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'};
for ii = 1:4
    cellID = intersect(find(strcmp(g.cell_metrics.animal, animals{ii})), find(g.cell_metrics.cluID == cluIDs(ii)));
    
    
    for jj = 1:3
        ax = nexttile;
        plotPSTHLocal(g.cell_metrics.spikes.times{cellID}, g.cell_metrics.general.(ttl{jj}){cellID}, [-0.5 1], 0.01);
        xlabel('Time (s)');
        if ii == 1
            if jj == 1
                ylabel({'Cluster 1 neuron', 'Firing Rate (Hz)'})
            end
            yticks([0 50 100])
            xticks([-0.5 0 0.5 1])
            ylim([0 100])
        elseif ii == 2
            if jj == 1
                ylabel({'Cluster 2 neuron', 'Firing Rate (Hz)'})
            end
            yticks([0 50 100])
            xticks([-0.5 0 0.5 1])
            ylim([0 100])
        elseif ii == 3
            if jj == 1
                ylabel({'Cluster 4 neuron', 'Firing Rate (Hz)'})
            end
            yticks([0 50 100])
            xticks([-0.5 0 0.5 1])
            ylim([0 100])
        elseif ii == 4
            if jj == 1
                ylabel({'Cluster 4 neuron', 'Firing Rate (Hz)'})
            end
            yticks([0 10 20])
            xticks([-0.5 0 0.5 1])
            ylim([0 20])
        end

        set(gca, 'FontSize', g.fontSize1);
        box off
    end
    ax = nexttile;
    ax.Visible = 'off';
end

mx = nexttile(t,12,[2,1]);
similarity_matrix = corr(PSTHall(leafOrder,:)');
imagesc(similarity_matrix); 
colormap(mx,'default')
colorbar

clearvars -except g t ttl


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
end






