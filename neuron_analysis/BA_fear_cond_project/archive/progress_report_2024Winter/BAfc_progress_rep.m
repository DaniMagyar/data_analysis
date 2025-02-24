% Clustering HABITUATION overall responses


pre_time = 5;
post_time = 10;
bin_time = 0.1;


[psth_spx] =  BAfc_psth_spx_v2('TTL_tone_habit_first', pre_time, post_time, bin_time);
[~, idx_exc, idx_inh] =  BAfc_find_sound_v04('TTL_tone_habit_first', 1:10);
PSTHall1 = smoothdata(zscore(psth_spx,[],2));
[~,PCA1] = pca(PSTHall1);
features1 = PCA1(:,1);
clustnum = 3;
[Kidx, Kcentroids] = kmeans(features1, clustnum, ...
    'Distance', 'sqeuclidean', 'Replicates', 10);
figure()
for ii=1:clustnum
    subplot(clustnum,1,ii)
    imagesc( 1:size(PSTHall1,2), 1:size(PSTHall1(find(Kidx==ii),:),1), PSTHall1(find(Kidx==ii),:));
    clim([-max(max(PSTHall1,[],1))/2.5 max(max(PSTHall1,[],1))]);
    mycolormap = getColormap;
    colormap(mycolormap); 
end


timewin = 1:100;
Dend      = linkage(features1,'complete','mahalanobis'); % calculating the dendrogram
Clusters  = cluster(Dend,'maxclust',clustnum); % clustering based on the tree
D         = pdist(features1); %euclidean distrance between point in the artifical space
leafOrder = optimalleaforder(Dend,D); % optimal order for plotting    
PCA        = features1;
Dendrogram = Dend;
figure; 
subplot(1,4,1:2)
imagesc(timewin, 1:size(PSTHall1,1), PSTHall1(leafOrder,:)); %plotting sorted psth matrix 
clim([-max(max(PSTHall1,[],1))/2.5 max(max(PSTHall1,[],1))]);
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


[psth_spx] =  BAfc_psth_spx_v2('TTL_tone_recall_first', pre_time, post_time, bin_time);
[~, idx_exc, idx_inh] =  BAfc_find_sound_v04('TTL_tone_recall_first', 1:10);
PSTHall2 = smoothdata(zscore(psth_spx,[],2));
[~,PCA1] = pca(PSTHall2);
features2 = PCA1(:,1);
clustnum = 3;
[Kidx, Kcentroids] = kmeans(features2, clustnum, ...
    'Distance', 'sqeuclidean', 'Replicates', 10);
figure()
for ii=1:clustnum
    subplot(clustnum,1,ii)
    imagesc( 1:size(PSTHall2,2), 1:size(PSTHall2(find(Kidx==ii),:),1), PSTHall2(find(Kidx==ii),:));
    clim([-max(max(PSTHall2,[],1))/2.5 max(max(PSTHall2,[],1))]);
    mycolormap = getColormap;
    colormap(mycolormap); 
end


timewin = 1:100;
Dend      = linkage(features2,'complete','mahalanobis'); % calculating the dendrogram
Clusters  = cluster(Dend,'maxclust',clustnum); % clustering based on the tree
D         = pdist(features2); %euclidean distrance between point in the artifical space
leafOrder = optimalleaforder(Dend,D); % optimal order for plotting    
PCA        = features2;
Dendrogram = Dend;
figure; 
subplot(1,4,1:2)
imagesc(timewin, 1:size(PSTHall2,1), PSTHall2(leafOrder,:)); %plotting sorted psth matrix 
clim([-max(max(PSTHall2,[],1))/2.5 max(max(PSTHall2,[],1))]);
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




[psth_spx] =  BAfc_psth_spx_v2('TTL_tone_cond_first', pre_time, post_time, bin_time);
[~, idx_exc, idx_inh] =  BAfc_find_sound_v04('TTL_tone_cond_first', 1:10);
PSTHall3 = smoothdata(zscore(psth_spx,[],2));
[~,PCA1] = pca(PSTHall3);
features3 = PCA1(:,1);
clustnum = 3;
[Kidx, Kcentroids] = kmeans(features3, clustnum, ...
    'Distance', 'sqeuclidean', 'Replicates', 10);
figure()
for ii=1:clustnum
    subplot(clustnum,1,ii)
    imagesc( 1:size(PSTHall3,2), 1:size(PSTHall3(find(Kidx==ii),:),1), PSTHall3(find(Kidx==ii),:));
    clim([-max(max(PSTHall3,[],1))/2.5 max(max(PSTHall3,[],1))]);
    mycolormap = getColormap;
    colormap(mycolormap); 
end


timewin = 1:100;
Dend      = linkage(features3,'complete','mahalanobis'); % calculating the dendrogram
Clusters  = cluster(Dend,'maxclust',clustnum); % clustering based on the tree
D         = pdist(features3); %euclidean distrance between point in the artifical space
leafOrder = optimalleaforder(Dend,D); % optimal order for plotting    
PCA        = features3;
Dendrogram = Dend;
figure; 
subplot(1,4,1:2)
imagesc(timewin, 1:size(PSTHall3,1), PSTHall3(leafOrder,:)); %plotting sorted psth matrix 
clim([-max(max(PSTHall3,[],1))/2.5 max(max(PSTHall3,[],1))]);
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
















































features3 = [features1 features3];
PSTHall3 = [PSTHall1 PSTHall3];

timewin = 1:100;
Dend      = linkage(features3,'complete','mahalanobis'); % calculating the dendrogram
Clusters  = cluster(Dend,'maxclust',clustnum); % clustering based on the tree
D         = pdist(features3); %euclidean distrance between point in the artifical space
%leafOrder = optimalleaforder(Dend,D); % optimal order for plotting    
PCA        = features3;
Dendrogram = Dend;
figure; 
subplot(1,4,1:2)
imagesc(timewin, 1:size(PSTHall3,1), PSTHall3(leafOrder,:)); %plotting sorted psth matrix 
clim([-max(max(PSTHall3,[],1))/2.5 max(max(PSTHall3,[],1))]);
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
