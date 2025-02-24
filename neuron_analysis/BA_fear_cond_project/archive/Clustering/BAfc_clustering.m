pre_time = 14.5; % in sec
post_time = 5;
baseline_time = 5;
bin_time = 0.001; % in seconds
clustnum = 4;
timewin = 1:pre_time/bin_time+post_time/bin_time;
respwin = pre_time/bin_time+1:pre_time/bin_time+post_time/bin_time;
baselinewin = (pre_time-baseline_time)/bin_time+1:pre_time/bin_time;
bin_factor = 20; % Adjust this factor as needed
% Colormap blue:white:red -------------------------------------------------
c1 = 1/255*[0,115,185];
c2 = 1/255*[239,239,239];
c3 = 1/255*[217,84,26];
mycolormap(1:20,1)  = linspace(c1(1),c2(1),20);
mycolormap(21:64,1) = linspace(c2(1),c3(1),44);
mycolormap(1:20,2)  = linspace(c1(2),c2(2),20);
mycolormap(21:64,2) = linspace(c2(2),c3(2),44);
mycolormap(1:20,3)  = linspace(c1(3),c2(3),20);
mycolormap(21:64,3) = linspace(c2(3),c3(3),44);  


for ii = 1:2
    if ii == 1
        [psth_spx_rec, idx_exc_r, idx_inh_r] =  BAfc_find_sound_v04('TTL_tone_recall_first', (1:10));
        [psth_spx_hab, idx_exc_h, idx_inh_h] =  BAfc_find_sound_v04('TTL_tone_habit_first', (1:10));
        %idx = unique([idx_exc_r; idx_exc_h; idx_inh_r; idx_inh_h]);
        idx = 1:size(psth_spx_rec,1)';
        PSTHall = psth_spx_rec(idx, timewin); 

    elseif ii == 2      
        PSTHall = psth_spx_hab(idx, timewin); 
    end               
    new_bin_count = size(PSTHall, 2) / bin_factor;
    % Check if downsampling is feasible
    if mod(size(PSTHall, 2), bin_factor) ~= 0
        error('Number of bins is not evenly divisible by bin_factor.');
    end   
    % Downsample by averaging every bin_factor bins
    if bin_factor > 1
        psth_matrix_downsampled = reshape(mean(reshape(PSTHall', bin_factor, [])), new_bin_count, [])';
    else
        psth_matrix_downsampled = PSTHall;
    end    
    % Normalize each row (optional)
    PSTHall = zscore(psth_matrix_downsampled, 0, 2);   
    % Smoothing (optional)    
    for i = 1:size(PSTHall,1)
        PSTHall(i,:) = smoothdata(PSTHall(i,:), 'movmean', 20);
    end

    if ii == 1 
        PSTHall_rec = PSTHall;
    elseif ii == 2 
        PSTHall_hab = PSTHall;
    end
end

dFR_baseline = zscore((sum(psth_spx_rec(idx, baselinewin),2) - sum(psth_spx_hab(idx, baselinewin),2))/(numel(baselinewin)*bin_time)); % firing rate in HZ
dFR_response = zscore((sum(psth_spx_rec(idx, respwin),2) - sum(psth_spx_hab(idx, respwin),2))/(numel(respwin)*bin_time));
PSTHall = [PSTHall_hab PSTHall_rec];

% PCA analysis ------------------------------------------------------------
[~,PCA1]  = pca(PSTHall); % running pricipal component analysis
PCA2      = PCA1(:,1:3); % extracting the first PCA_num pricinpal components 
inNAN = nanmedian(PCA2,1); % to avoid NaN values, set them to the median
for jj = 1:size(PSTHall,1)
    if isnan(PCA2(jj,1)) == 1
        PCA2(jj,1:3) = inNAN;
    end
end

features = PCA2;
features(:,3) = dFR_baseline;
features(:,4) = dFR_response;



Dend      = linkage(features,'complete','mahalanobis'); % calculating the dendrogram
Clusters  = cluster(Dend,'maxclust',clustnum); % clustering based on the tree
D         = pdist(features); %euclidean distrance between point in the artifical space
leafOrder = optimalleaforder(Dend,D); % optimal order for plotting    
PCA        = features;
Dendrogram = Dend;
figure; 
subplot(1,4,1:2)
imagesc(timewin, 1:size(PSTHall,1), PSTHall(leafOrder,:)); %plotting sorted psth matrix 
clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
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



clustnum = 5;
[Kidx, Kcentroids] = kmeans(features, clustnum, ...
    'Distance', 'sqeuclidean', 'Replicates', 500);
figure()
for ii=1:clustnum
    subplot(clustnum,1,ii)
    imagesc(timewin, 1:size(PSTHall(find(Kidx==ii),:),1), PSTHall(find(Kidx==ii),:));
    clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
    colormap(mycolormap); 
end

cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'0'};
cell_metrics.labels(idx) = cellstr(num2str(Clusters));

cell_metrics.labels(idx) = cellstr(num2str(Kidx));


cell_metrics = CellExplorer('metrics',cell_metrics);