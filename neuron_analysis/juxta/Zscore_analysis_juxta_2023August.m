function Zscore_analysis_juxta_2023August(varargin)
% Default params
prs =  inputParser;
addParameter(prs,'mainFolder','Z:\HajosLab\Dani\Magyar_Daniel\experiments\juxta_BLA',@ischar) % Path of main folder. (e.g. 'C:\Users\dmagyar\Desktop\M2_shock_response')
addParameter(prs,'Stim','none',@ischar) % TTLs to use from TTLsKS.mat. (e.g. 'shocks')
addParameter(prs,'loadBatch','none',@ischar) % Use batch file instead of original data, default 'none'. (e.g. 'cell_metrics_batch.mat')
addParameter(prs,'Recordings',@iscell) % Use original files; list of included recordings in cell format. (e.g. {'MD127', 'MD128', 'MD129'}).
addParameter(prs,'selVariable','none',@ischar) % Variable used for neuron selection from cell_metrics. (e.g. 'labels')
addParameter(prs,'selValue','none',@ischar) % Value used for neuron selection from 'selectVariable'. (e.g. '2')
addParameter(prs,'neuronType','none',@ischar) % Options: 'Narrow Interneuron', 'Wide Interneuron', 'Pyramidal Cell'. Default 'none' loads all types. 
addParameter(prs,'plotType','Heatmap',@ischar) % Options: 'Heatmap', 'Lineplot'
addParameter(prs,'Sorted',0,@isnumeric) % If 1, Load 'SortIDX' from 'BAparams.mat'. Default 0.
addParameter(prs,'fs',20000,@isnumeric) % Sampling rate, default 30000.
addParameter(prs,'psth_bin',400,@isnumeric) % Binning window calculated from fs. (e.g. 1ms=30, 20ms=600(20x30), 50ms=1500(50x30))
addParameter(prs,'int',[-0.5 0.5],@isnumeric) % Plotting window [-pre post] in seconds.  
addParameter(prs,'Wcx_win',[-0.05 0.05],@isnumeric) % Time window [-pre post] for Wilcoxon in seconds. Two sides should be equal.
addParameter(prs,'Wcx_alpha',0.05,@isnumeric) % alpha value for Wilcoxon, default 0.05
addParameter(prs,'norm',1,@isnumeric) % Normalise data (Z-score). Options: 1 or 0.
addParameter(prs,'average',0,@isnumeric) % Plot excited/inhibited group averages. Options: 1 or 0. 
addParameter(prs,'SignificantZscore',0,@isnumeric) % Defines significance based on provided Z-score value. Default 0.
addParameter(prs,'PCA_num',2,@isnumeric) % Number of principal components used
addParameter(prs,'clustnum',4,@isnumeric) % Number of clusters to create
addParameter(prs,'LKmethod','complete',@ischar)
addParameter(prs,'LKmetric','mahalanobis',@ischar)
parse(prs,varargin{:})
g = prs.Results;
disp(g.mainFolder)
PSTHall=[]; 
Wilcoxon_results = []; 
% % % % % % Load cell_metrics structure
% % % % % if strcmp(g.loadBatch, 'none')
% % % % %     for ii = 1:length(g.Recordings)
% % % % %         basepaths(ii) = {[g.mainFolder '\' g.Recordings{ii} '_kilosort\kilosort25preprocess']};
% % % % %         basenames(ii) = {'temp_wh'};
% % % % %     end
% % % % %     cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
% % % % % else
% % % % %     cell_metrics = importdata(g.loadBatch);
% % % % % end
% % % % % % Select neurons from cell_metrics based on structure variable
% % % % % if strcmp(g.selVariable, 'none')
% % % % %     cellIdx_selVariable = 1:numel(cell_metrics.cellID); % all neurons
% % % % % else
% % % % %     cellIdx_selVariable = find(strcmp(cell_metrics.(g.selVariable), g.selValue)); 
% % % % % end
% % % % % % Select neuron type
% % % % % if strcmp(g.neuronType, 'none')
% % % % %     cellIdx_type = 1:numel(cell_metrics.cellID); % all neurons
% % % % % else
% % % % %     cellIdx_type = find(strcmp(cell_metrics.putativeCellType, g.neuronType)); 
% % % % % end
% % % % % % Find neurons that match selected variable and celltype
% % % % % cellIdx_selected = intersect(cellIdx_selVariable, cellIdx_type);
% % % % % % Load PSTH matrix
% % % % % prevBatch = 0;
loadBar = waitbar(0,'Loading neurons...');

Tab = dir('selected*'); 
Tab = {Tab.name}';
cellIdx_selected   = 1:size(Tab,1);

for hh = cellIdx_selected
    clear AP 
    load(Tab{hh})


    AP = Spikes;
    TTL = Shocks;
% % % %     if cell_metrics.batchIDs(hh) ~= prevBatch
% % % %         load([cell_metrics.general.basepaths{cell_metrics.batchIDs(hh)} '\TTLsKS.mat']);
% % % %         TTL = eval(g.Stim);
% % % %         prevBatch = cell_metrics.batchIDs(hh);
% % % %     end
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
     psth_spx_dani = vertcat(sum(preAP_bin,2), sum(postAP_bin,2));
     clear preAP_bin postAP_bin preAP postAP preAP_norm postAP_norm  % reset matrices, because different TTL num    
     % Each post point cointain all spikes during the whole illumination period. (Like one 2 secs long bin)
     for wcx = 1:numel(TTL) %Each TTL is a a column. 
         preAP_wcx{:,wcx} = AP(AP>=(TTL(wcx)-abs(g.Wcx_win(1))) & AP<TTL(wcx)); % spikes before each TTL separately
         postAP_wcx{:,wcx} = AP(AP>TTL(wcx) & AP<(TTL(wcx)+g.Wcx_win(2))); % spikes after each TTL separately
     end
     preAP_wcx_num = cellfun(@numel, preAP_wcx);
     postAP_wcx_num = cellfun(@numel, postAP_wcx);
     preAP_wcx_freq = preAP_wcx_num/abs(g.Wcx_win(1));
     postAP_wcx_freq = postAP_wcx_num/g.Wcx_win(2);
     [~,h] = signrank(preAP_wcx_freq, postAP_wcx_freq, 'alpha', g.Wcx_alpha);
     Wilcoxon_results{find(cellIdx_selected==hh),1} = hh; 
     Wilcoxon_results{find(cellIdx_selected==hh),2} = num2str(h);
     clear preAP_wcx postAP_wcx preAP_wcx_num postAP_wcx_num preAP_wcx_freq postAP_wcx_freq
            if g.norm == 1
                newpsth = zscore(psth_spx_dani);
            elseif g.norm == 0
                newpsth = psth_spx_dani;
            end
            PSTHall(find(cellIdx_selected==hh),:) = newpsth';
%            PSTHallIDs(find(cellIdx_selected==hh),1) = {['Batch ' num2str(cell_metrics.batchIDs(hh)) ' Cell ' num2str(cell_metrics.cellID(hh))]};
end
close(loadBar)
time = (g.int(1)*g.fs:g.psth_bin:g.int(2)*g.fs)';
time = time(1:end-1);
timeline = time/g.fs+g.psth_bin/g.fs;
switch g.plotType
    case 'Heatmap'
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
    % Calculating plotting order based on Z-score change
    testWindow_firstBin = abs(g.int(1))*g.fs/g.psth_bin+1;
    testBins = g.Wcx_win(2)*g.fs/g.psth_bin; % number of bins in Wcx test window. Z-score ordering, Wcx direction based on this
    testWindow_lastBin = testWindow_firstBin + (testBins-1);    
    testWindow = PSTHall(:, testWindow_firstBin:testWindow_lastBin);
    testMean = mean(testWindow,2);
    switch g.Sorted
        case 0
            [~,SortIDX] = sort(testMean, 'descend');
%            MyNewOrder = PSTHallIDs(SortIDX);
        case 1
            load ([g.mainFolder '\BAparams.mat'], 'SortIDX','MyNewOrder')
    end
    figure; 
    subplot(1,4,1:2)
    imagesc(time/g.fs, 1:size(PSTHall,1), PSTHall(SortIDX,:)); % plotting sorted psth matrix 
    clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
    colormap(mycolormap); 
    cb = colorbar('FontSize',15);
    cb.Ticks = linspace(-20,20,41);
    cb.Label.String = 'Z-score change';
    cb.Label.FontSize = 25;
    cb.Position = [0.55, 0.115, 0.01, 0.255];
    hold on;
    plot([-0.00125 -0.00125],[1-0.5 size(PSTHall,1)+0.5],'r', 'LineWidth', 2);
    hold off;
    ylabel('# Cell');
    xlabel('Time (s)');
    set(gca,'FontSize',35);   
    Wilcox_sorted = str2num(cell2mat(Wilcoxon_results(:,2)));
    Wilcox_sorted = Wilcox_sorted(SortIDX); 
    significant_idx = find(Wilcox_sorted == 1);
    non_significant_idx = find(Wilcox_sorted == 0);
    subplot(1,4,3:4)
    plot(zeros(1,length(non_significant_idx)), non_significant_idx, 'b.', 'MarkerSize', 30)
    hold on
    plot(zeros(1,length(significant_idx)), significant_idx, 'r.', 'MarkerSize', 30)
    hold off
    set(gca,'ydir','reverse')
    set(gca,'visible','off')
    set(subplot(1,4,3:4), 'Position', [0.5, 0.115, 0.01, 0.805])
    xlim([-0.5 0.5])
    ylim([1 length(Wilcox_sorted)])   
    testMeanSorted = testMean(SortIDX);
    for qq = 1:length(testMean)
        if Wilcox_sorted(qq) == 1
            if testMeanSorted(qq) < -g.SignificantZscore
                response_dir(qq) = -1;
            elseif testMeanSorted(qq) >= g.SignificantZscore
                response_dir(qq) = 1;
            end
        else
            response_dir(qq) = 0;
        end
    end        
    annotation('textbox', [0.55, 0.9, 1, 0], 'string',...
        ['Significantly modulated units:' num2str(numel(find(response_dir == 1))+numel(find(response_dir == -1)))], 'LineStyle','none', 'FontSize', 35)
    annotation('textbox', [0.55, 0.7, 1, 0], 'string',...
        ['Excited: ' num2str(numel(find(response_dir == 1))) '  ('...
        num2str(numel(find(response_dir == 1))/length(testMean)*100) '%)'], 'LineStyle','none', 'FontSize', 35)
    annotation('textbox', [0.55, 0.5, 1, 0], 'string',...
        ['Inhibited: ' num2str(numel(find(response_dir == -1))) '  ('...
        num2str(numel(find(response_dir == -1))/length(testMean)*100) '%)'], 'LineStyle','none', 'FontSize', 35)
    
    set(gcf,'Position',[500 50 1600 1300])
        switch g.Sorted
            case 0
%                save ([g.mainFolder '\BAparams.mat'], 'SortIDX','MyNewOrder', 'response_dir')
            case 1
                load ([g.mainFolder '\BAparams.mat'], 'response_dir') 
        end   
    % Plot the averages
    if g.average == 1   
        PSTHallSorted = PSTHall(SortIDX,:);
        Clustmeans(1,:) = mean(PSTHallSorted((find(response_dir==1)),:),1); % excited avg
        Clustmeans(2,:) = mean(PSTHallSorted((find(response_dir==-1)),:),1); % inhibited avg
        Clustmeans(1,:) = Clustmeans(1,:) - mean(Clustmeans(1,1:abs(g.int(1)*0.8/(g.psth_bin/g.fs)))); %offset, 0.8: PV start excluded
        Clustmeans(2,:) = Clustmeans(2,:) - mean(Clustmeans(2,1:abs(g.int(1)*0.8/(g.psth_bin/g.fs)))); %offset
        for tt = 1:2
            figure(tt+1)
            hold on
            plot(timeline, Clustmeans(tt,:), 'LineWidth', 2.5) % time shifted (different timing is necessary from heatmap's)
            hold off
            ylabel('Z-scored firing rate');
            xlabel('Time (s)');
            set(gca,'FontSize',35);
            set(gcf,'Position',[500 50 1600 1300])
            if tt == 1
                title('Average response of excited neurons')
            elseif tt == 2
                title('Average response of inhibited neurons')
            end    
        end
    end
    case 'Lineplot'
        plot (timeline,PSTHall) % plot individual neurons
        hold on
        plot(timeline, mean(PSTHall,1), 'Color', [.4 .4 .4], 'LineWidth', 3) % plot mean
        hold off
        xlabel('Time (s)');
        xlim([timeline(1) timeline(end)])
        if g.norm == 1
            ylabel('Z-scored firing rate');
        elseif g.norm == 0
            ylabel('non Z-scored firing rate');
        end
    case 'PCA'
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
        plot([-0.00125 -0.00125],[1 size(PSTHall,1)],'r');
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

end
end

