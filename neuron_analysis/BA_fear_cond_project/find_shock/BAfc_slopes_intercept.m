function [PSTHall, data] = BAfc_slopes_intercept(ttl, trial_total, trial_per_block)

pre_time = 0.4;
post_time = 1;
test_time = 0.4;
bin_time = 0.05;
% trial_total = 50;
% trial_per_block = 5; % 20 must be dividible with this number
smoothvalue = 8;

for ii = 1:trial_total/trial_per_block
    psth_spx{ii} =  BAfc_psth_spx(ttl, ((ii-1)*trial_per_block+1:ii*trial_per_block), pre_time, post_time, bin_time);
end


psth_spx_all = horzcat(psth_spx{:});
psth_spx_zscore = zscore(psth_spx_all,0,2);

numSplits = trial_total/trial_per_block;
psth_spx_zscore_all = cell(1, numSplits);

% Split the matrix into 5000-column chunks and store in cells
for i = 1:numSplits
    colStart = (i-1) * (pre_time+post_time)/bin_time + 1;       % Start column index
    colEnd = i * (pre_time+post_time)/bin_time;                % End column index
    psth_spx_zscore_cell{i} = psth_spx_zscore(:, colStart:colEnd); % Extract and store
    psth_spx_all_cell{i} = psth_spx_all(:, colStart:colEnd);
end


AUC_all_cells = cell(1, 5);

% Loop through each cell
for i = 1:length(psth_spx_zscore_cell)
    % Extract the test columns of the current cell matrix
    current_data = psth_spx_zscore_cell{i}(:, pre_time/bin_time+1:(pre_time+test_time)/bin_time);

    % Initialize a vector to store AUC values for the current cell
    AUC_current_cell = zeros(size(current_data, 1), 1);
    
    % Loop through each neuron (row) in the current cell
    for j = 1:size(current_data, 1)
        % Calculate the AUC for the neuron (row) using trapezoidal rule
        AUC_current_cell(j) = trapz(current_data(j, :));  % Calculate AUC for each neuron
    end
    
    % Store the AUC values for the current cell
    AUC_all_cells{i} = AUC_current_cell;
end

AUC_matrix = horzcat(AUC_all_cells{:});

% Initialize arrays to store slopes and intercepts
slopes = zeros(size(AUC_matrix, 1), 1);
intercepts = zeros(size(AUC_matrix, 1), 1);

% Loop through each row (neuron) in AUC_matrix
for i = 1:size(AUC_matrix, 1)
    % Get the AUC values for the current neuron (row)
    neuron_AUC = AUC_matrix(i, :);
    
    % Perform linear regression (polyfit) on the AUC values
    % x = 1:length(neuron_AUC) represents the time steps (1 to number of cells)
    [p, S] = polyfit(1:length(neuron_AUC), neuron_AUC, 1); % Linear fit (degree = 1)
    
    % Store the slope and intercept
    slopes(i) = p(1);  % Slope of the linear regression
    intercepts(i) = p(2);  % Intercept of the linear regression
end

for iii = 1:size(psth_spx_zscore,1)
    PSTHall(iii,:) = smoothdata(psth_spx_zscore(iii,:), 'gaussian', smoothvalue);
end

% count number of spikes in response window
for kk = 1:size(psth_spx_all_cell,2)
    num_spikes_per_block(:,kk) = sum(psth_spx_all_cell{kk}(:, pre_time/bin_time+1:(pre_time+test_time)/bin_time),2);
end


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


zslopes = zscore(slopes);
zintercepst = zscore(intercepts);
data = [zintercepst, zslopes, zscore(AUC_matrix), zscore(num_spikes_per_block,[],2)]; % az AUC_matrix kiszedi a lowFR sejteket amiknek csak 1 trail miatt van slopejuk

k = 5; % Number of clusters
[idx, C] = kmeans(data, k); % 'idx' contains cluster labels

scatter(data(:,1), data(:,2), 50, idx, 'filled'); % Color points by cluster
xlabel('Property 1');
ylabel('Property 2');

figure()
for ii=1:k
    subplot(k,1,ii)
    imagesc( 1:size(PSTHall,2), 1:size(PSTHall(find(idx==ii),:),1), PSTHall(find(idx==ii),:));
    clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
    colormap(mycolormap); 
end

% 
% cell_metrics = BAfc_load_neurons;
% cell_metrics.labels(1:end) = {'0'};
% 
% Clusters = idx;
% cell_metrics.labels(1:end) = cellstr(num2str(idx)); % kmeans
% %cell_metrics.labels(intersect(find(zintercepts>1.3), find(zslopes<-1))) = cellstr(num2str(-1.3)); % without kmeans
% cell_metrics.SlopeShock = slopes';
% cell_metrics.InterceptShock = intercepts';
% cell_metrics.general.psth_spx_zscore = psth_spx_zscore_cell;
% cell_metrics = CellExplorer('metrics',cell_metrics);
% 
% 
% 



% 
% 
% clustnum = 6;
% % PCA analysis ------------------------------------------------------------
% [~,PCA1]  = pca(PSTHall); % running pricipal component analysis
% PCA2      = PCA1(:,1:3); % extracting the first PCA_num pricinpal components 
% inNAN = nanmedian(PCA2,1); % to avoid NaN values, set them to the median
% for jj = 1:size(PSTHall,1)
%     if isnan(PCA2(jj,1)) == 1
%         PCA2(jj,1:3) = inNAN;
%     end
% end
% 
% features = PCA2;
% timewin = 1:pre_time/bin_time+post_time/bin_time;
% Dend      = linkage(features,'complete','mahalanobis'); % calculating the dendrogram
% Clusters  = cluster(Dend,'maxclust',clustnum); % clustering based on the tree
% D         = pdist(features); %euclidean distrance between point in the artifical space
% leafOrder = optimalleaforder(Dend,D); % optimal order for plotting    
% PCA        = features;
% Dendrogram = Dend;
% figure; 
% subplot(1,4,1:2)
% imagesc(timewin, 1:size(PSTHall,1), PSTHall(leafOrder,:)); %plotting sorted psth matrix 
% clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
% colormap(mycolormap); 
% hold on;
% %plot([-0.00125 -0.00125],[1 size(PSTHall,1)],'r');
% hold off;
% ylabel('# Cell');
% xlabel('Time (s)') 
% subplot(1,4,3:4)
% cutoff = Dend(end-clustnum+2,3); % cutting the tree to get 'clustnum' clusters
% h = dendrogram(Dend,0,'Orientation','right','reorder',leafOrder,...
%     'ColorThreshold',cutoff); %plotting dendrogram
% set(h,'LineWidth',1)
% set(gca,'Ydir','reverse');
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gca,'visible','off')
% 


%[~,leafOrder] = sort(intercepts, 'descend'); hierarchical plotting