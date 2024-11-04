function Zscore_analysis_CellExp(Recordings, Stim, varargin)
% Input
% -Recordings: vertical list of included recordings in cell format, first 5 char
% -Stim: variable to use from TTLsKS.mat, string
% Example:
% -Zscore_analysis_KS({'MD130', 'MD131', 'MD132'}, 'BA_25_5Hz')

% Default params ----------------------------------------------------------
fs       = 30000; % Firing rate
int      = [-2 2]; % psth interval, even numbers are better for plot
norm     = 1; % Normalise data
psth_bin = 600; % 600 = 20ms; 1500=50ms
testBins = 25; % bins included for Z-score ordering !!!!!! WCX direction based on this, Must match WCX window
Wcx_win = [-0.5 0.5]; % Wilcoxon window in seconds
SignificantZscore = 0; % absolute value, if 0, not included.
average = 0; % plot averages, 1 if yes
mainFolder = 'C:\Users\dmagyar\Desktop\M2_shock_response';
minFR = 0;
disp(mainFolder)
PSTHall=[];
Wilcoxon_results = [];
% Load cell_metrics structure
for ii = 1:length(Recordings)
    basepaths(ii) = {[mainFolder '\' Recordings{ii} '_kilosort\kilosort25preprocess']};
    basenames(ii) = {'temp_wh'};
end
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

% Select neurons for PSTH; row number in cell_metrics, not the cellID

cellIdx_selected = 1:numel(cell_metrics.cellID); % all neurons

%cellIdx_selected = find(strcmp(cell_metrics.putativeCellType, 'Wide Interneuron'));
%cellIdx_selected = find(strcmp(cell_metrics.putativeCellType, 'Pyramidal Cell'));
% cellIdx_selected = find(strcmp(cell_metrics.putativeCellType, 'Narrow Interneuron'));
% firingRates = cell_metrics.firingRate(cellIdx_selected);
% cellIdx_selected(find(firingRates<10)) = [];



% Load PSTH matrix
prevBatch = 0;
loadBar = waitbar(0,'Loading neurons...');
for hh = cellIdx_selected
    AP = cell_metrics.spikes.times{hh};
    if cell_metrics.batchIDs(hh) ~= prevBatch
        load([cell_metrics.general.basepaths{cell_metrics.batchIDs(hh)} '\TTLsKS.mat']);
        TTL = eval(Stim);
        prevBatch = cell_metrics.batchIDs(hh);
    end
    waitbar(hh/numel(cellIdx_selected),loadBar);
    % PSTH matrix
    bin_time = psth_bin/fs;     
    pre_time = abs(int(1));      
    post_time = int(2);      
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
    
    
     %% paired Wilcoxon signed rank test: equal length of pre and post data. 
     % Each post point cointain all spikes during the whole illumination period. (Like one 2 secs long bin)
     for wcx = 1:numel(TTL) %Each TTL is a a column. 
         preAP_wcx{:,wcx} = AP(AP>=(TTL(wcx)-abs(Wcx_win(1))) & AP<TTL(wcx)); % spikes before each TTL separately
         postAP_wcx{:,wcx} = AP(AP>TTL(wcx) & AP<(TTL(wcx)+Wcx_win(2))); % spikes after each TTL separately
     end
     preAP_wcx_num = cellfun(@numel, preAP_wcx);
     postAP_wcx_num = cellfun(@numel, postAP_wcx);
     preAP_wcx_freq = preAP_wcx_num/abs(Wcx_win(1));
     postAP_wcx_freq = postAP_wcx_num/Wcx_win(2);
     [~,h] = signrank(preAP_wcx_freq, postAP_wcx_freq, 'alpha', 0.05);
     Wilcoxon_results{find(cellIdx_selected==hh),1} = hh; 
     Wilcoxon_results{find(cellIdx_selected==hh),2} = num2str(h);
     clear preAP_wcx postAP_wcx preAP_wcx_num postAP_wcx_num preAP_wcx_freq postAP_wcx_freq

            if norm == 1
                newpsth = zscore(psth_spx_dani);
            else
                newpsth = psth_spx_dani;
            end
            PSTHall(find(cellIdx_selected==hh),:) = newpsth';
            PSTHallIDs(find(cellIdx_selected==hh),1) = {['Batch ' num2str(cell_metrics.batchIDs(hh)) ' Cell ' num2str(cell_metrics.cellID(hh))]};

end
close(loadBar)
time = (int(1)*fs:psth_bin:int(2)*fs)';
time = time(1:end-1);

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

%% Calculating plotting order based on Z-score change
testWindow_firstBin = abs(int(1))*fs/psth_bin+1;
testWindow_lastBin = testWindow_firstBin + (testBins-1);    
testWindow = PSTHall(:, testWindow_firstBin:testWindow_lastBin);
testMean = mean(testWindow,2);
switch Stim
    case {'BA_50_5Hz', 'BA_50_5Hz_Resting','BA_50_10Hz_Resting', 'ChETA_50_20Hz', 'BA_25_5Hz', 'BA_250_5Hz', 'BA_25_10Hz', 'BA_250_10Hz', ...
            'shock_only', 'shock_only_Resting', 'BA_25_all_Resting', 'BA_25_5Hz_Resting', 'shocks', 'shock_motor_rest', ...
            'BA_5Hz_motor_rest_first', 'BA_20Hz_motor_rest_first', 'BA_5Hz_motor_rest_all', 'BA_20Hz_motor_rest_all', 'stim'}
        [~,SortIDX] = sort(testMean, 'descend');
        MyNewOrder = PSTHallIDs(SortIDX);
    case {'BA_50_5Hz_Running', 'BA_50_10Hz_Running', 'TO_25_5Hz', 'TO_250_5Hz', 'TO_25_10Hz','shock_inh', 'shock_only_Running' ,...
            'shock_inh_Running', 'shock_inh_Resting', 'BA_25_all_Running', 'shock_motor_run',...
            'BA_5Hz_motor_run_first', 'BA_20Hz_motor_run_first', 'BA_5Hz_motor_run_all', 'BA_20Hz_motor_run_all'}
        load ([mainFolder '\BAparams.mat'], 'SortIDX','MyNewOrder')
end
%% Common part

figure; 
subplot(1,4,1:2)
imagesc(time/fs, 1:size(PSTHall,1), PSTHall(SortIDX,:)); % plotting sorted psth matrix 
clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
colormap(mycolormap); 
cb = colorbar('FontSize',15);
cb.Ticks = linspace(-20,20,41);
cb.Label.String = 'Z-score change';
cb.Label.FontSize = 25;
cb.Position = [0.55, 0.115, 0.01, 0.255];
hold on;
plot([-0.00125 -0.00125],[1 size(PSTHall,1)],'r', 'LineWidth', 2);
% plot(int,[min(find(testMean(SortIDX)<0))-0.5 min(find(testMean(SortIDX)<0))-0.5], 'r') % finding the first negative Z-socre mean
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
        if testMeanSorted(qq) < -SignificantZscore
            response_dir(qq) = -1;
        elseif testMeanSorted(qq) >= SignificantZscore
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
    switch Stim
        case {'BA_50_5Hz', 'BA_50_5Hz_Resting', 'BA_50_10Hz_Resting', 'BA_25_5Hz', 'BA_250_5Hz', 'BA_25_10Hz', 'BA_250_10Hz',...
                'shock_only', 'shock_only_Resting', 'BA_25_all_Resting', 'BA_25_5Hz_Resting', 'shocks', 'shock_motor_rest',...
                'BA_5Hz_motor_rest_first', 'BA_20Hz_motor_rest_first', 'stim'}
            save ([mainFolder '\BAparams.mat'], 'SortIDX','MyNewOrder', 'response_dir')
        case {'BA_50_5Hz_Running', 'BA_50_10Hz_Running', 'TO_25_5Hz', 'TO_250_5Hz', 'TO_25_10Hz','shock_inh', 'shock_only_Running', ...
                'shock_inh_Running', 'shock_inh_Resting', 'BA_25_all_Running', 'shock_motor_run',...
                'BA_5Hz_motor_run_first', 'BA_20Hz_motor_run_first'}
            load ([mainFolder '\BAparams.mat'], 'response_dir') % for average plots
    end

%% Plot the averages
if average == 1

    PSTHallSorted = PSTHall(SortIDX,:);
    Clustmeans(1,:) = mean(PSTHallSorted((find(response_dir==1)),:),1); % excited avg
    Clustmeans(2,:) = mean(PSTHallSorted((find(response_dir==-1)),:),1); % inhibited avg
    Clustmeans(1,:) = Clustmeans(1,:) - mean(Clustmeans(1,1:abs(int(1)*0.8/(psth_bin/fs)))); %offset, 0.8: PV start excluded
    Clustmeans(2,:) = Clustmeans(2,:) - mean(Clustmeans(2,1:abs(int(1)*0.8/(psth_bin/fs)))); %offset
    for tt = 1:2
        figure(tt+1)
        hold on
        plot(time/fs+psth_bin/fs, Clustmeans(tt,:), 'LineWidth', 2.5) % time shifted (different timing is necessary from heatmap's)
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
end

