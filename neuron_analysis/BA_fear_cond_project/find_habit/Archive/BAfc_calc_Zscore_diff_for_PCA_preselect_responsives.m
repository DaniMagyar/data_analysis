function [cell_metrics] = BAfc_calc_Zscore_diff_for_PCA_preselect_responsives(varargin)

prs =  inputParser;
addParameter(prs,'mainFolder','C:\Users\dmagyar\Desktop\BA_fear_cond',@ischar) % Path of main folder. (e.g. 'C:\Users\dmagyar\Desktop\M2_shock_response')
addParameter(prs,'Stims',@iscell) % TTLs to use from TTLsKS.mat. (e.g. 'shocks')
addParameter(prs,'TTLselect',[0],@isnumeric) % Selected TTLs (e.g. [1:3])
addParameter(prs,'selVariable','none',@ischar) % Variable used for neuron selection from cell_metrics. (e.g. 'labels')
addParameter(prs,'selValue','none',@ischar) % Value used for neuron selection from 'selectVariable'. (e.g. '2')
addParameter(prs,'neuronType','none',@ischar) % Options: 'Narrow Interneuron', 'Wide Interneuron', 'Pyramidal Cell'. Default 'none' loads all types. 
addParameter(prs,'fs',30000,@isnumeric) % Sampling rate, default 30000.
addParameter(prs,'psth_bin',3000,@isnumeric) % Binning window calculated from fs. (e.g. 1ms=30, 20ms=600(20x30), 50ms=1500(50x30))
addParameter(prs,'int',[-4.5 4.5],@isnumeric) % Plotting window [-pre post] in seconds.  
addParameter(prs,'Wcx_win',[-5 5],@isnumeric) % Time window [-pre post] for Wilcoxon in seconds. Two sides should be equal.
addParameter(prs,'Wcx_alpha',0.05,@isnumeric) % alpha value for Wilcoxon, default 0.05
addParameter(prs,'norm',1,@isnumeric) % Normalise data (Z-score). Options: 1 or 0.
addParameter(prs,'offset',1,@isnumeric) % Offset to pre stimulus firing rate ( int(1)).
addParameter(prs,'Sorted',0,@isnumeric) % If 1, Load 'SortIDX' from 'BAparams.mat'. Default 0.
addParameter(prs,'SignificantZscore',0,@isnumeric) % Defines significance based on provided Z-score value. Default 0.
addParameter(prs,'PCA_num',2,@isnumeric) % Number of principal components used
addParameter(prs,'clustnum',8,@isnumeric) % Number of clusters to create
addParameter(prs,'LKmethod','complete',@ischar)
addParameter(prs,'LKmetric','euclidean',@ischar)
parse(prs,varargin{:})
g = prs.Results;

cell_metrics = BAfc_load_neurons;

% Select neurons from cell_metrics based on structure variable
if strcmp(g.selVariable, 'none')
    cellIdx_selVariable = 1:numel(cell_metrics.cellID); % all neurons
else
    cellIdx_selVariable = find(strcmp(cell_metrics.(g.selVariable), g.selValue)); 
end
% Select neuron type
if strcmp(g.neuronType, 'none')
    cellIdx_type = 1:numel(cell_metrics.cellID); % all neurons
else
    cellIdx_type = find(strcmp(cell_metrics.putativeCellType, g.neuronType)); 
end
% Remove 'junk' 
cellIdx_keep =  1:numel(cell_metrics.cellID); % all neurons
cellIdx_keep(find(strcmp(cell_metrics.putativeCellType, 'junk'))) = []; % neurons labeled as 'junk'
% Find neurons that match selected variable and celltype and not junk
cellIdx_selected = intersect(intersect(cellIdx_selVariable, cellIdx_type), cellIdx_keep);

% Load PSTH matrix
prevBatch = 0;
loadBar = waitbar(0,'Loading neurons...');
Significance(1:numel(cellIdx_selected)) = 0;
for hh = cellIdx_selected
    TTL = [];
    AP = cell_metrics.spikes.times{hh};
    for ii = 1:length(g.Stims)
        ttl = cell_metrics.general.(g.Stims{ii}){hh};
        TTL = vertcat(TTL,ttl);
    end
    TTL = sort(TTL);
    if g.TTLselect ~= 0 
        TTL = TTL(g.TTLselect);
    end
    waitbar(hh/numel(cellIdx_selected),loadBar);
    % PSTH matrix
    bin_time = g.psth_bin/g.fs;     
    pre_time = abs(g.int(1));      
    post_time = g.int(2);      
    for jj = 1:numel(TTL) %Each TTL is a a column. 
         preAP{:,jj} = AP(AP>=(TTL(jj)-pre_time) & AP<TTL(jj)); % spikes before each TTL separately
         postAP{:,jj} = AP(AP>TTL(jj) & AP<(TTL(jj)+post_time)); % spikes after each TTL separately
    end
    for ll = 1:numel(TTL) %Each TTL is a a column. 
         preAP_norm{ll} = preAP{ll}-TTL(ll); % spikes relative to their own TTL
         postAP_norm{ll} = postAP{ll}-TTL(ll);
    end
    for tt = 1:numel(TTL) %Each TTL is a a column. 
        for nn = 1:(pre_time/bin_time) % number of timestamps in each bin.
             preAP_bin(nn,tt) = sum(preAP_norm{tt}>=(-pre_time+(nn-1)*bin_time) & preAP_norm{tt}<(-pre_time+nn*bin_time));
        end
        for oo = 1:(post_time/bin_time)
             postAP_bin(oo,tt) = sum(postAP_norm{tt}>=((oo-1)*bin_time) & postAP_norm{tt}<(oo*bin_time));
        end
    end  
    
    preFR = cellfun(@numel, preAP_norm)/pre_time;
    postFR = cellfun(@numel, postAP_norm)/post_time;
    diffFR = postFR - preFR;
    diffFRZscore = (diffFR-mean(preFR))/std(preFR);


    clear preAP postAP preAP_norm postAP_norm  % reset matrices, because different TTL num  
    temp = [(preFR-mean(preFR))/std(preFR); (diffFR-mean(preFR))/std(preFR)];
    temp2 = temp(:)';
    temp3 = [preFR-mean(preFR); postFR-mean(preFR)];
    temp4 = temp3(:)';
    PSTHall(find(cellIdx_selected==hh),:) = temp4;
    PSTHallIDs(find(cellIdx_selected==hh),1) = {['Batch ' num2str(cell_metrics.batchIDs(hh)) ' Cell ' num2str(cell_metrics.cellID(hh))]};
    PSTHallIDs(find(cellIdx_selected==hh),2) = cell_metrics.brainRegion(hh);

    for pp = 1:size(preAP_bin,2)
        [pRanksum(pp), resRanksum(pp)] = ranksum(preAP_bin(:,pp), postAP_bin(:,pp), 'alpha', 0.05);
    end

    if any(resRanksum == 1)
        resRanksumAll(hh) = 1;
    else
        resRanksumAll(hh) = 0;
    end
end

idx_ranksum = find(resRanksumAll==1);
PSTHall = PSTHall(idx_ranksum,:);



close(loadBar)
time = (g.int(1)*g.fs:g.psth_bin:g.int(2)*g.fs)';
time = time(1:end-1);
timeline = time/g.fs+g.psth_bin/g.fs;
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


% PCA analysis ------------------------------------------------------------
[~,PCA1]  = pca(PSTHall); % running pricipal component analysis
PCA2      = PCA1(:,1:g.PCA_num); % extracting the first PCA_num pricinpal components 
inNAN = nanmedian(PCA2,1); % to avoid NaN values, set them to the median
for jj = 1:size(PSTHall,1)
    if isnan(PCA2(jj,1)) == 1
        PCA2(jj,1:3) = inNAN;
    end
end
Dend      = linkage(PCA2,g.LKmethod,g.LKmetric); % calculating the dendrogram
Clusters  = cluster(Dend,'maxclust',g.clustnum); % clustering based on the tree
D         = pdist(PCA2); %euclidean distrance between point in the artifical space
leafOrder = optimalleaforder(Dend,D); % optimal order for plotting    
PCA        = PCA2;
Dendrogram = Dend;
figure; 
subplot(1,4,1:2)
imagesc(time/g.fs, 1:size(PSTHall,1), PSTHall(leafOrder,:)); %plotting sorted psth matrix 
clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
colormap(mycolormap); 
hold on;
%plot([-0.00125 -0.00125],[1 size(PSTHall,1)],'r');
hold off;
ylabel('# Cell');
xlabel('Time (s)') 
subplot(1,4,3:4)
cutoff = Dend(end-g.clustnum+2,3); % cutting the tree to get 'clustnum' clusters
h = dendrogram(Dend,0,'Orientation','right','reorder',leafOrder,...
    'ColorThreshold',cutoff); %plotting dendrogram
set(h,'LineWidth',1)
set(gca,'Ydir','reverse');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'visible','off')

cell_metrics.labels(1:end) = {'0'};
cell_metrics.labels(idx_ranksum) = cellstr(num2str(Clusters));