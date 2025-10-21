% function BAfc_figure_4
% % Plotting region specific responses on heatmaps, proper clustering
clear all
% recordings = {...
%     'MD292_002_kilosort',...
%     'MD293_kilosort',...
%     'MD294_kilosort',...
%     'MD295_kilosort',...
%     'MD296_kilosort',...
%     'MD297_kilosort',...
%     'MD298_kilosort',...
%     'MD299_kilosort',...
%     'MD300_kilosort',...
%     'MD304_kilosort'};

recordings = {...
    'MD292_002_kilosort',...
    'MD293_001_kilosort',...
    'MD294_001_kilosort',...
    'MD295_001_kilosort',...
    'MD296_001_kilosort',...
    'MD297_001_kilosort',...
    'MD298_001_kilosort',...
    'MD299_001_kilosort',...
    'MD300_001_kilosort',...
    'MD304_001_kilosort',...
    'MD307_001_kilosort',...
    'MD309_001_kilosort',...
    'MD311_002_kilosort',...
    'MD312_001_kilosort',...
    'MD313_001_kilosort',...
    'MD314_001_kilosort',...
    'MD315_001_kilosort',...
    'MD316_002_kilosort',...
    'MD317_001_kilosort',...
    'MD318_001_kilosort',...
    'MD318_002_kilosort',...
    'MD319_003_kilosort'};

g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'});

clearvars -except g
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 13;
g.fontSize2 = 12;
g.xlinewidth = 2;
g.pre_time = 5;
g.post_time = 5;
g.test_time = 0.5; % 1 sec is good, beacuse captures better than 0.5
g.testvalue = 3; % if the mean change is negative, but there is a peak where zscore is above this, than counted as excited
g.bin_time = 0.01; % 0.01 volt sokat 5 os smoothal
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
t = tiledlayout(fig,(2+g.clustnum),4,'TileSpacing', 'compact', 'Padding', 'none');

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
PCA2      = PCA1(:,1:3); % extracting the first PCA_num pricinpal components 
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
    psth_spx_sorted = [psth_spx_PN(leafOrderfl,:)];
    matrix = psth_spx_sorted(:,(g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    imagesc(g.timeaxis,1:size(matrix,1),matrix);
    clim(g.clim)
    disp([-max(max(psth_spx_sorted,[],1))/2.5 max(max(psth_spx_sorted,[],1))])   
    colormap(g.colors.Heatmap); 
    yline(sum(idx_PN)+0.5, 'k', 'LineWidth', 1);
    xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
    ylabel('Cell number')
    set(gca, 'FontSize', g.fontSize2);
    title(hmptitles{hmp}, 'FontSize', g.fontSize1)
    if hmp == 1
        cb = colorbar('westoutside', 'FontSize', g.fontSize2);
    end
    hold on
    n_clu = find(diff(Clusters(leafOrderfl))~=0);
    yline(n_clu+0.5, 'Color', 'k', 'LineWidth', 1);
    hold off
    
    strt = 1;
    for clus = 1:g.clustnum
        barloc = (8:4:(8+g.clustnum*(g.clustnum-1)))+hmp;
        allspikes = [];
        if clus == 1
            allspikes = sum(psth_spx_sorted(1:n_clu(clus),:),1);
        elseif clus == g.clustnum
            allspikes = sum(psth_spx_sorted(n_clu(clus-1)+1:end,:),1);
        else
            allspikes = sum(psth_spx_sorted(n_clu(clus-1)+1:n_clu(clus),:),1);
        end
        nexttile(t, barloc(clus), [1 1])
        b = bar(g.timeaxis(2:end),allspikes(:,(g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time), 'k'); % 'k' for black bars
        b.FaceColor =  g.colors.barcolor;
        b.EdgeColor = g.colors.barcolor; 
        ylim([-30 100])
        yticks([-25 0 50 100])
        if hmp == 1
            ylabel(['Cluster ' num2str(clus) ' average'])
        end
        if clus == 4
            xlabel('Time(s)')
        end
        set(gca, 'FontSize', g.fontSize2);
    end



end
nexttile(t,4,[2 1])
cutoff = Dend(end-g.clustnum+2,3); % cutting the tree to get 'clustnum' clusters
h = dendrogram(Dend,0,'Orientation','right','reorder',flipud(leafOrder),...
    'ColorThreshold',cutoff); %plotting dendrogram
axis off


