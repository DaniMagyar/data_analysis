function [PSTHall, data] = BAfc_slopes_intercept_v3(varargin)
% v3 different methond, AUC not used anymore

% Default params
prs =  inputParser;
addRequired(prs,'ttl',@ischar)
addRequired(prs,'trial_total',@isnumeric) 
addRequired(prs,'trial_per_block',@isnumeric) 
addParameter(prs,'pre_time',0, @isnumeric)
addParameter(prs,'post_time',1, @isnumeric)
addParameter(prs,'test_time',0.2, @isnumeric)
addParameter(prs,'bin_time',0.05, @isnumeric)
addParameter(prs,'smoothvalue',5, @isnumeric)
addParameter(prs,'idxS',[], @isnumeric) % previously selected indices
parse(prs,varargin{:})
g = prs.Results;

parfor ii = 1:g.trial_total/g.trial_per_block
    psth_spx{ii} =  BAfc_psth_spx_v2(g.ttl, g.pre_time, g.post_time, g.bin_time, 'TTLinclude', ((ii-1)*g.trial_per_block+1:ii*g.trial_per_block));
    psth_spx_base{ii} = BAfc_psth_spx_v2(g.ttl, 10, 10, g.bin_time, 'TTLinclude', ((ii-1)*g.trial_per_block+1:ii*g.trial_per_block), 'TTLshift', -10);
end
 % calculating zscore per block, each block has it's own baseline block
for jj = 1:size(psth_spx,2)
    for kk = 1:size(psth_spx{jj},1)
        psth_spx_zscore{jj}(kk,:) = (psth_spx{jj}(kk,:) -  mean(psth_spx_base{jj}(kk,:)))/std(psth_spx_base{jj}(kk,:));
    end
    spxNumTestWin(:,jj) = sum(psth_spx{jj}(:,g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time),2);
    peakZscoreTestWin(:,jj) = max(psth_spx_zscore{jj}(:,g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time),[],2);
end

spxNumZ = zscore(spxNumTestWin,[],2);
% slope of zscores (we need zscore because if the beaseline activity
% increases, zscore reamins more stable if neuron non-responsive)
for i = 1:size(spxNumZ, 1) 
    neuronZ = spxNumZ(i, :);
    [p, S] = polyfit(1:length(neuronZ), neuronZ, 1); 
    slopes(i,1) = p(1); 
    intercepts(i,1) = p(2);  
end

plotdata = zscore(horzcat(psth_spx{:}),[],2);
% for plotting
for iii = 1:size(plotdata,1)
    PSTHall(iii,:) = smoothdata(plotdata(iii,:), 'gaussian', g.smoothvalue);
end

zslopes = zscore(slopes);
zintercepst = zscore(intercepts);



if isempty(g.idxS)
    g.idxS = 1:size(PSTHall,1);
    disp('input indices not found')
end

[~,PCA1]  = pca(zscore(spxNumZ(g.idxS,:),[],2)); % running pricipal component analysis
PCA2      = PCA1(:,1:2); % extracting the first PCA_num pricinpal components 

%data = [zintercepst(g.idxS), zslopes(g.idxS), zscore(spxNumZ(g.idxS,:),[],2), spxNumZ(g.idxS,:)]; % az AUC_matrix kiszedi a lowFR sejteket amiknek csak 1 trail miatt van slopejuk
data = [PCA2];
data = [zintercepst(g.idxS) zslopes(g.idxS)];

k = 3; % Number of clusters
[idx, C] = kmeans(data, k); % 'idx' contains cluster labels
PSTHall_idxS = PSTHall(g.idxS,:);

figure()
for ii=1:k
    subplot(k,1,ii)
    imagesc( 1:size(PSTHall_idxS,2)/20+1, 1:size(PSTHall_idxS(find(idx==ii),:),1), PSTHall_idxS(find(idx==ii),:));
    clim([-max(max(PSTHall_idxS,[],1))/2.5 max(max(PSTHall_idxS,[],1))]);
    mycolormap = getColormap;
    colormap(mycolormap); 
end
figure()
scatter(data(:,1), data(:,2), 50, idx, 'filled'); % Color points by cluster
xlabel('Property 1');
ylabel('Property 2');

% 
% cell_metrics = BAfc_load_neurons;
% cell_metrics.labels(1:end) = {'0'};
% 
% 
% cell_metrics.labels(g.idxS) = cellstr(num2str(idx)); % kmeans
% %cell_metrics.labels(intersect(find(zintercepts>1.3), find(zslopes<-1))) = cellstr(num2str(-1.3)); % without kmeans
% cell_metrics.SlopeShock = slopes';
% cell_metrics.InterceptShock = intercepts';
% cell_metrics.general.psth_spx_zscore = psth_spx_zscore_cell;
% cell_metrics = CellExplorer('metrics',cell_metrics);
% 
% 
% 





clustnum = 3;
% PCA analysis ------------------------------------------------------------
[~,PCA1]  = pca(PSTHall_idxS); % running pricipal component analysis
PCA2      = PCA1(:,1); % extracting the first PCA_num pricinpal components 
inNAN = nanmedian(PCA2,1); % to avoid NaN values, set them to the median
for jj = 1:size(PSTHall_idxS,1)
    if isnan(PCA2(jj,1)) == 1
        PCA2(jj,1:3) = inNAN;
    end
end

features = PCA2;
timewin = 1:11;
Dend      = linkage(features,'complete','mahalanobis'); % calculating the dendrogram
Clusters  = cluster(Dend,'maxclust',clustnum); % clustering based on the tree
D         = pdist(features); %euclidean distrance between point in the artifical space
leafOrder = optimalleaforder(Dend,D); % optimal order for plotting    
PCA        = features;
Dendrogram = Dend;
figure; 
subplot(1,4,1:2)
imagesc(timewin, 1:size(PSTHall_idxS,1), PSTHall_idxS(leafOrder,:)); %plotting sorted psth matrix 
clim([-max(max(PSTHall_idxS,[],1))/2.5 max(max(PSTHall_idxS,[],1))]);
colormap(mycolormap); 
hold on;
%plot([-0.00125 -0.00125],[1 size(PSTHall,1)],'r');
hold off;
ylabel('# Cell');
xlabel('Time (s)') 
subplot(1,4,3:4)
cutoff = Dend(end-clustnum+2,3); % cutting the tree to get 'clustnum' clusters
h = dendrogram(Dend,0,'Orientation','right','reorder',leafOrder,...
    'ColorThreshold',cutoff); %plotting dendrogram
set(h,'LineWidth',1)
set(gca,'Ydir','reverse');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'visible','off')



%[~,leafOrder] = sort(intercepts, 'descend'); hierarchical plotting