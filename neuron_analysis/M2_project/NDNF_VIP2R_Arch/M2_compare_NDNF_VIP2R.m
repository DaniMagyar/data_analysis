% MD281,MD282, MD2823
basenames = {'temp_wh', 'temp_wh','temp_wh','temp_wh','temp_wh','temp_wh', 'temp_wh','temp_wh','temp_wh','temp_wh', 'temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD281_001_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD281_002_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD281_003_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD281_004_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD281_005_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD282_001_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD282_002_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD283_002_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD283_003_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

cell_metrics.general.TTL_shocks_only = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks_only;
    cell_metrics.general.TTL_shocks_only(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_shocks_light = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks_light;
    cell_metrics.general.TTL_shocks_light(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.labels(1:end) = {'neutral'};
pre_time = 5; % in sec
post_time = 5;    
bin_time = 0.001; % in seconds
% first TTL
% Preallocate arrays for results
num_cells = numel(cell_metrics.cellID);
psth_spx_1 = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly
smoothed_spike_train_zscore = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly
inputTTL = 'TTL_shocks_only';
parfor ii = 1:num_cells
    disp(ii)
    TTL = cell_metrics.general.(inputTTL){ii};
    AP = cell_metrics.spikes.times{ii};

    % Preallocate cell arrays for spikes
    preAP = cell(numel(TTL), 1);
    postAP = cell(numel(TTL), 1);
    preAP_norm = cell(numel(TTL), 1);
    postAP_norm = cell(numel(TTL), 1);

    % Binning variables
    preAP_bin = zeros(pre_time/bin_time, numel(TTL));
    postAP_bin = zeros(post_time/bin_time, numel(TTL));

    % Loop over TTL events
    for jj = 1:numel(TTL)
        % Spikes before and after each TTL
        preAP{jj} = AP(AP >= (TTL(jj) - pre_time) & AP < TTL(jj));
        postAP{jj} = AP(AP > TTL(jj) & AP < (TTL(jj) + post_time));

        % Normalize spike times to TTL
        preAP_norm{jj} = preAP{jj} - TTL(jj);
        postAP_norm{jj} = postAP{jj} - TTL(jj);
    end

    % Bin spike times
    for tt = 1:numel(TTL)
        preAP_bin(:, tt) = histcounts(preAP_norm{tt}, linspace(-pre_time, 0, pre_time/bin_time + 1));
        postAP_bin(:, tt) = histcounts(postAP_norm{tt}, linspace(0, post_time, post_time/bin_time + 1));
    end

    % Sum bins across TTLs
    psth_spx_1(ii,:) = sum([preAP_bin; postAP_bin], 2);

    % Smoothing and Z-scoring
    sigma = 20; % standard deviation of Gaussian kernel (ms)
    binWidth = 1; % bin width (ms)
    edges = -3*sigma:binWidth:3*sigma; % kernel edges
    gaussKernel = exp(-edges.^2 / (2*sigma^2)); % Gaussian kernel
    gaussKernel = gaussKernel / sum(gaussKernel); % Normalize
    smoothed_spike_train = conv(psth_spx_1(ii,:), gaussKernel, 'same'); % Convolve
    smoothed_spike_train_zscore(ii,:) = (smoothed_spike_train - mean(smoothed_spike_train(1:4000))) / std(smoothed_spike_train(1:4000));
end
% Find responses in 500ms time window
window200 = smoothed_spike_train_zscore(:,5001:5500);
[idx_exc,~] = find(window200>=5);
idx_exc = unique(idx_exc);
cell_metrics.labels(idx_exc) = {'exc'};
% second TTL
% Preallocate arrays for results
num_cells = numel(cell_metrics.cellID);
psth_spx_2 = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly
smoothed_spike_train_zscore = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly
inputTTL = 'TTL_shocks_light';
parfor ii = 1:num_cells
    disp(ii)
    TTL = cell_metrics.general.(inputTTL){ii};
    AP = cell_metrics.spikes.times{ii};

    % Preallocate cell arrays for spikes
    preAP = cell(numel(TTL), 1);
    postAP = cell(numel(TTL), 1);
    preAP_norm = cell(numel(TTL), 1);
    postAP_norm = cell(numel(TTL), 1);

    % Binning variables
    preAP_bin = zeros(pre_time/bin_time, numel(TTL));
    postAP_bin = zeros(post_time/bin_time, numel(TTL));

    % Loop over TTL events
    for jj = 1:numel(TTL)
        % Spikes before and after each TTL
        preAP{jj} = AP(AP >= (TTL(jj) - pre_time) & AP < TTL(jj));
        postAP{jj} = AP(AP > TTL(jj) & AP < (TTL(jj) + post_time));

        % Normalize spike times to TTL
        preAP_norm{jj} = preAP{jj} - TTL(jj);
        postAP_norm{jj} = postAP{jj} - TTL(jj);
    end

    % Bin spike times
    for tt = 1:numel(TTL)
        preAP_bin(:, tt) = histcounts(preAP_norm{tt}, linspace(-pre_time, 0, pre_time/bin_time + 1));
        postAP_bin(:, tt) = histcounts(postAP_norm{tt}, linspace(0, post_time, post_time/bin_time + 1));
    end

    % Sum bins across TTLs
    psth_spx_2(ii,:) = sum([preAP_bin; postAP_bin], 2);

    % Smoothing and Z-scoring
    sigma = 20; % standard deviation of Gaussian kernel (ms)
    binWidth = 1; % bin width (ms)
    edges = -3*sigma:binWidth:3*sigma; % kernel edges
    gaussKernel = exp(-edges.^2 / (2*sigma^2)); % Gaussian kernel
    gaussKernel = gaussKernel / sum(gaussKernel); % Normalize
    smoothed_spike_train = conv(psth_spx_2(ii,:), gaussKernel, 'same'); % Convolve
    smoothed_spike_train_zscore(ii,:) = (smoothed_spike_train - mean(smoothed_spike_train(1:4000))) / std(smoothed_spike_train(1:4000));
end
% Find responses in 500ms time window
window200 = smoothed_spike_train_zscore(:,5001:5500);
[idx_exc,~] = find(window200>=5);
idx_exc = unique(idx_exc);
cell_metrics.labels(idx_exc) = {'exc'};

prewin = 2501:5000;
poswin = 5001:7500;

numTTL = cell2mat(cellfun(@numel,cell_metrics.general.TTL_shocks_only,'UniformOutput',false)');
dFR1 = (sum(psth_spx_1(:, poswin),2) - sum(psth_spx_1(:, prewin),2))./numTTL;
dFR2 = (sum(psth_spx_2(:, poswin),2) - sum(psth_spx_2(:, prewin),2))./numTTL;
increase = find(dFR2>(dFR1+1));
x = dFR1(intersect(increase, idx_exc));
y = dFR2(intersect(increase, idx_exc));

% Assuming sumresp1 and sumresp2 are already calculated // myzscore
data1 = sum(psth_spx_1(intersect(increase, idx_exc), :),1);
data2 = sum(psth_spx_2(intersect(increase, idx_exc), :),1);
sumresp1 = smoothdata((data1 - mean(data1(4500:5000)))/std(data1(4500:5000)), 'movmean', 50);
sumresp2 = smoothdata((data2 - mean(data2(4500:5000)))/std(data2(4500:5000)), 'movmean', 50);

% Define the range you want to plot
range = 4500:6500;

% Create X-axis values for plotting
x_values = range;
x_values = linspace(-0.5, 1.5, length(range)); % From -0.5s to 1.5s over the range

% Plot sumresp1 as a line and fill the area under it with transparency
fill([x_values, fliplr(x_values)], [sumresp1(range), zeros(1, length(range))], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;

% Plot sumresp2 as a line and fill the area under it with transparency
fill([x_values, fliplr(x_values)], [sumresp2(range), zeros(1, length(range))], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the lines on top for clarity
plot(x_values, sumresp1(range), 'b', 'LineWidth', 2);
plot(x_values, sumresp2(range), 'r', 'LineWidth', 2);

% Add labels and title for clarity
xlabel('Time (s)');
ylabel('Z-scored response');
% title('Smoothed Z-scored Responses with Shaded Area');
set(gca, 'FontSize', 14); % Set the font size for the axis labels and ticks
legend({'Tail-shock', 'Tail-shock + Laser ON'}, 'FontSize', 14, 'Location', 'northeast');
text(1.2, 3, 'n=38', 'FontSize', 16);

% Customize the plot for better visibility
%grid on;
hold off;


cell_metrics.labels(intersect(increase, idx_exc)) = {'selected'};




%cell_metrics = CellExplorer('metrics',cell_metrics);
% cellID = 94;
% rasterPsthSingle(cell_metrics.spikes.times{cellID},cell_metrics.general.TTL_shocks_only{cellID}, 'window',[-2 2],'bin_time', 0.02)
% rasterPsthSingle(cell_metrics.spikes.times{cellID},cell_metrics.general.TTL_shocks_light{cellID}, 'window',[-2 2],'bin_time', 0.02)

% % Colormap blue:white:red -------------------------------------------------
% c1 = 1/255*[0,115,185];
% c2 = 1/255*[239,239,239];
% c3 = 1/255*[217,84,26];
% mycolormap(1:20,1)  = linspace(c1(1),c2(1),20);
% mycolormap(21:64,1) = linspace(c2(1),c3(1),44);
% mycolormap(1:20,2)  = linspace(c1(2),c2(2),20);
% mycolormap(21:64,2) = linspace(c2(2),c3(2),44);
% mycolormap(1:20,3)  = linspace(c1(3),c2(3),20);
% mycolormap(21:64,3) = linspace(c2(3),c3(3),44);  
% window = 4751:5250;
% PSTHall = smoothdata(zscore([psth_spx_1(idx_exc250, window)],[],2),2, 'movmean', 20);
% % PSTHall = smoothdata(zscore([psth_spx_2(idx_exc250, window)],[],2),2, 'movmean', 50)...
% %     - smoothdata(zscore([psth_spx_1(idx_exc250, window)],[],2),2, 'movmean', 50);
% PCA_num = 2;
% clustnum = 3;
% % PCA analysis ------------------------------------------------------------
% [~,PCA1]  = pca(PSTHall); % running pricipal component analysis
% PCA2      = PCA1(:,1:PCA_num); % extracting the first PCA_num pricinpal components 
% inNAN = nanmedian(PCA2,1); % to avoid NaN values, set them to the median
% for jj = 1:size(PSTHall,1)
%     if isnan(PCA2(jj,1)) == 1
%         PCA2(jj,1:3) = inNAN;
%     end
% end
% Dend      = linkage(PCA2,'complete','euclidean'); % calculating the dendrogram
% Clusters  = cluster(Dend,'maxclust',clustnum); % clustering based on the tree
% D         = pdist(PCA2); %euclidean distrance between point in the artifical space
% leafOrder = optimalleaforder(Dend,D); % optimal order for plotting    
% PCA        = PCA2;
% Dendrogram = Dend;
% figure; 
% subplot(1,4,1:2)
% imagesc(1:500, 1:size(PSTHall,1), PSTHall(leafOrder,:)); %plotting sorted psth matrix 
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
% 
% 
% 
% cell_metrics.labels(idx_exc250(find(Clusters == 2))) = {'exc3'};
% sumresp1 = sum(psth_spx_1(idx_exc250(find(Clusters == 3)), :),1);
% sumresp2 = sum(psth_spx_2(idx_exc250(find(Clusters == 3)), :),1);
% bar(sumresp1(4951:5250), 'b','FaceAlpha', 0.5)
% hold on
% bar(sumresp2(4951:5250), 'r','FaceAlpha', 0.5)
% hold off
% 

% 
% sumresp1 = smoothdata(zscore(sum(psth_spx_1(idx_exc250(find(Clusters == 1)), :),1),[],2), 'movmean', 10);
% sumresp2 = smoothdata(zscore(sum(psth_spx_2(idx_exc250(find(Clusters == 1)), :),1),[],2), 'movmean', 10);
% bar(sumresp1(4501:5500), 'b','FaceAlpha', 0.5)
% hold on
% bar(sumresp2(4501:5500), 'r','FaceAlpha', 0.5)
% hold off

