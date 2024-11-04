function Zscore_analysis_KS(Recordings, Stim)
% Input
% -Recordings: vertical list of included recordings in cell format, first 5 char
% -Stim: variable to use from TTLsKS.mat, string

% Zscore_analysis_KS({'MD130', 'MD131', 'MD132'}, 'BA_25_5Hz')
% Zscore_analysis_KS({'MD133', 'MD134', 'MD135', 'MD136', 'MD137'}, 'ChETA_50_20Hz')
% Zscore_analysis_KS({'MD138', 'MD139','MD140', 'MD141', 'MD142', 'MD143'}, 'shock_only')
% Zscore_analysis_KS({'MD138', 'MD139','MD140', 'MD141', 'MD142', 'MD143'}, 'shock_inh')
% Zscore_analysis_KS({'MD127', 'MD128', 'MD129'}, 'shocks')

% Default params ----------------------------------------------------------
fs       = 30000; % Firing rate
int      = [-0.025 0.025]; % psth interval, even numbers are better for plot
norm     = 1; % Normalise data
psth_bin = 60; % 600 = 20ms; 1500=50ms
testBins = 5; % bins included for Z-score ordering !!!!!! WCX direction based on this, Must match WCX window
Wcx_win = [-0.01 0.01]; % Wilcoxon window in seconds
SignificantZscore = 0; % absolute value, if 0, not included.
average = 0; % plot averages, 1 if yes
mainFolder = 'C:\Users\dmagyar\Desktop\M2_shock_response';
minFR = 0;
disp(mainFolder)
PSTHall=[];
Wilcoxon_results = [];

for ii = 1:length(Recordings)
    cd([mainFolder '\' Recordings{ii} '_kilosort\kilosort25preprocess'])

    cluster_group = tdfread('cluster_group.tsv');
    cluster_info = tdfread('cluster_info.tsv');
%     if length(cluster_group.group) ~= length(cluster_info.group)
%         error('Not all clusters have group label in Phy')
%     else
%         disp('All groups labeled')
%     end
    
    spike_times = readNPY('spike_times.npy');
    spike_clusters = readNPY('spike_clusters.npy');
    
    groupCell = cellstr(cluster_info.group);
    goodIdx = find(strcmp(groupCell, 'good')==1); % to include MUA: find(strcmp(groupCell, 'noise')==0)

    for lowFR = 1:length(goodIdx)
        totalTime = double((spike_times(length(spike_times)) - spike_times(1))/30000);
        spikeNum  = numel(find(spike_clusters(spike_clusters == goodIdx(lowFR))));
        if spikeNum/totalTime < minFR
            goodIdx(lowFR) = NaN;
        end
    end
    goodIdx = goodIdx(~isnan(goodIdx));

    goodClusterIDs = cluster_info.id(goodIdx); %not used in this script KS 25
    %goodClusterIDs = cluster_info.cluster_id(goodIdx); %KS3

    
    for kk = 1:length(goodIdx)
        % Load TTLs & Choose variable to use as stimulus
        load([cd '\TTLsKS.mat']); % #ok<LOAD>
        switch Stim
            case 'BA_50_5Hz'
                ttl = BA_50_5Hz;
            case 'BA_50_5Hz_Running'
                ttl = BA_50_5Hz_Running;   
            case 'BA_50_5Hz_Resting'
                ttl = BA_50_5Hz_Resting;
                ttl = BA_50_5Hz_Resting(1:length(BA_50_5Hz_Running));  
            case 'BA_50_10Hz_Running'
                ttl = BA_50_10Hz_Running;   
            case 'BA_50_10Hz_Resting'
                ttl = BA_50_10Hz_Resting;
                ttl = BA_50_10Hz_Resting(1:length(BA_50_10Hz_Running));  
            case 'ChETA_50_20Hz'
                ttl = ChETA_50_20Hz;
            case 'BA_500'
                ttl = BA_500;
            case 'BA_500_5Hz'
                ttl = BA_500_5Hz;
            case 'BA_50'
                ttl = BA_50;
            case 'TTL_500'
                ttl = TTL;
            case 'TTL_50'
                ttl = TTL(1:10:end);
            case 'BA_25'
                ttl = BA_25;
            case 'BA_25_5Hz'
                ttl = BA_25_5Hz;
    %            ttl = BA_25_5Hz(11:25);
            case 'BA_25_10Hz'
                ttl = BA_25_10Hz;
            case 'TO_25'
                ttl = TO_25;
            case 'TO_25_5Hz'
                ttl = TO_25_5Hz;
    %            ttl = TO_25_5Hz(1:10);
            case 'TO_25_10Hz'
                ttl = TO_25_10Hz;
            case 'BA_250'
                ttl = BA_250;
            case 'BA_250_5Hz'
                %ttl = BA_250_5Hz;
                %ttl = vertcat(BA_250_5Hz(1:10:end),[]);
                ttl = vertcat(BA_250_5Hz(1:10:end), BA_250_5Hz(2:10:end));
                %ttl = vertcat(BA_250_5Hz(1:10:end), BA_250_5Hz(2:10:end), BA_250_5Hz(3:10:end));
                ttl = ttl;
                %ttl = BA_250_5Hz(1:125);
%                 ttl = [BA_250_5Hz(6:10:end)];
%                 ttl = [ttl; BA_250_5Hz(7:10:end)];
%                 ttl = [ttl; BA_250_5Hz(8:10:end)];
%                 ttl = [ttl; BA_250_5Hz(9:10:end)];
                ttl = sort(ttl);
            case 'BA_250_10Hz'
                ttl = BA_250_10Hz;
            case 'TO_250'
                ttl = TO_250;
            case 'TO_250_5Hz'
                ttl = TO_250_5Hz;
    %            ttl = TO_250_5Hz(1:125);
            case 'TO_250_10Hz'
                ttl = TO_250_10Hz;
            case 'BA_500_5Hz_10Hz'
                ttl = BA_250_5Hz;
                ttl = [ttl; BA_250_10Hz];
                ttl = sort(ttl);
            case 'SK'
                ttl = TTL; % The SK TTL is already called TTL.
            case 'shock_only'
                ttl = shock_only;
            case 'shock_only_Running'
               % ttl = shock_only_Running;
               if length(shock_only_Resting) < length(shock_only_Running)
                   ttl = shock_only_Running(1:length(shock_only_Resting));
               else 
                   ttl = shock_only_Running;
               end
            case 'shock_only_Resting'
               % ttl = shock_only_Resting;
               if length(shock_only_Running) < length(shock_only_Resting)
                   ttl = shock_only_Resting(1:length(shock_only_Running));
               else 
                   ttl = shock_only_Resting;
               end
            case 'shock_inh_Running'
                ttl = shock_inh_Running;
            case 'shock_inh_Resting'
                ttl = shock_inh_Resting;
            case 'shock_inh'
                ttl = shock_inh;
            case 'shocks'
                ttl = shocks;
            case 'BA_25_all_Resting'
                ttl = BA_25_all_Resting;
            case 'BA_25_all_Running'
                ttl = BA_25_all_Running;
            case 'BA_25_5Hz_Resting'
                ttl = BA_25_5Hz_Resting;
                %ttl = sort(vertcat(ttl, ttl+0.2, ttl+0.4, ttl+0.6, ttl+0.8));
                %ttl = sort(vertcat(ttl, ttl+0.2, ttl+0.4, ttl+0.6, ttl+0.8, ttl+1, ttl+1.2, ttl+1.4, ttl+1.6, ttl+1.8));
            case 'shock_motor_rest'
                ttl = shock_motor_rest;
            case 'shock_motor_run'
                ttl = shock_motor_run;
            case 'BA_5Hz_motor_rest_first'
                ttl = BA_5Hz_motor_rest_first;
            case 'BA_20Hz_motor_rest_first'
                ttl = BA_20Hz_motor_rest_first;
            case 'BA_5Hz_motor_run_first' 
                ttl = BA_5Hz_motor_run_first;
            case 'BA_20Hz_motor_run_first'
                ttl = BA_20Hz_motor_run_first;
            case 'BA_5Hz_motor_rest_all'
                ttl = BA_5Hz_motor_rest_all;
            case 'BA_20Hz_motor_rest_all'
                ttl = BA_20Hz_motor_rest_all;
            case 'BA_5Hz_motor_run_all' 
                ttl = BA_5Hz_motor_run_all;
            case 'BA_20Hz_motor_run_all'
                ttl = BA_20Hz_motor_run_all;
            case 'stim'
                ttl = stim;
        end
    
        %ttl = ttl(1:10);
        %spk_idx = find(spike_clusters==cluster_info.cluster_id(goodIdx(kk))); %KS 3
        spk_idx = find(spike_clusters==cluster_info.id(goodIdx(kk))); %KS 2.5
        spk_t = double(spike_times(spk_idx));
    
        TT = spk_t/fs;
    
        % PSTH matrix
        bin_time = psth_bin/fs;     
        pre_time = abs(int(1));      
        post_time = int(2);      
        AP = TT;
        TTL = ttl;
        for jj = 1:numel(TTL) %Each TTL is a a column. 
             preAP{:,jj} = AP(AP>=(TTL(jj)-pre_time) & AP<TTL(jj)); % spikes before each TTL separately
             postAP{:,jj} = AP(AP>TTL(jj) & AP<(TTL(jj)+post_time)); % spikes after each TTL separately
             %postAP{:,jj} = AP(AP>(TTL(jj)+0.006) & AP<(TTL(jj)+0.006+post_time)); %BA_250
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
             %postAP_wcx{:,wcx} = AP(AP>(TTL(wcx)+0.006) & AP<(TTL(wcx)+0.006+Wcx_win(2))); %BA_250
         end
         preAP_wcx_num = cellfun(@numel, preAP_wcx);
         postAP_wcx_num = cellfun(@numel, postAP_wcx);
         preAP_wcx_freq = preAP_wcx_num/abs(Wcx_win(1));
         postAP_wcx_freq = postAP_wcx_num/Wcx_win(2);
         [~,h] = signrank(preAP_wcx_freq, postAP_wcx_freq, 'alpha', 0.05);
         %Wilcoxon_resCurr{kk,1} = cluster_info.cluster_id(kk); %KS 3
         Wilcoxon_resCurr{kk,1} = cluster_info.id(kk); %KS 2.5
         Wilcoxon_resCurr{kk,2} = num2str(h);
         clear preAP_wcx postAP_wcx preAP_wcx_num postAP_wcx_num preAP_wcx_freq postAP_wcx_freq
         % real-time scatter plot
% close all hidden;
% x1 = (0.01:0.01:2.5);
% x2 = (4.01:0.01:6.5);
% y1 = preAP_wcx_freq;
% y2 = postAP_wcx_freq;
% 
% figure
% hold on
% scatter(x1,y1, 30, 'r', 'filled');
% scatter(x2,y2, 30, 'g', 'filled');
% 
% plot([x1(:)';x2(:)'], [y1(:)';y2(:)'], 'k-')
         
    %    % Normality test included
    %      SWresults{kk,1} = NeuronID;
    %      [SWresults{kk,2},SWresults{kk,3},~] = swtest(preAP_wcx_freq, 0.05, true);
    %      [SWresults{kk,4},SWresults{kk,5},~] = swtest(postAP_wcx_freq, 0.05, true);
    %      if SWresults{kk,2} == 1
    %          if SWresults{kk,4} == 1
    %              h = ttest(preAP_wcx_freq, postAP_wcx_freq);
    %              Wilcoxon_results{kk,2} = num2str(h);
    %          end
    %      end
    
         
    
        % remove laser artefact bin (psth_spx length match PSTHall length)
%             psth_spx_dani(numel(psth_spx_dani)+1)=mean(psth_spx_dani); % add extra mean bin  
%             psth_spx_dani(abs(int(1)*fs/psth_bin)+1)=[]; %to remove first bin after TTL
            
        
            if norm == 1
                newpsth = zscore(psth_spx_dani);
            else
                newpsth = psth_spx_dani;
            end
            PSTHcurr(kk,:) = newpsth';
    end
    PSTHall = vertcat(PSTHall, PSTHcurr);
    Wilcoxon_results =  vertcat(Wilcoxon_results, Wilcoxon_resCurr);
    clear PSTHcurr Wilcoxon_resCurr;

end

% rZ = 1;
% while rZ <= length(PSTHall)
%     if PSTHall(rZ,:) == 0
%         PSTHall(rZ,:) = [];
%         Wilcoxon_results(rZ,:) = [];
%     end
%     rZ = rZ +1;
% end

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
        MyNewOrder = [];
%       MyNewOrder = goodClusterIDs(SortIDX);
%             MyNewOrder.Recording = Tab.Recording(SortIDX);
%             MyNewOrder.Group = Tab.Group(SortIDX);
%             MyNewOrder.Neuron = Tab.Neuron(SortIDX);
%             MyNewOrder.Type = Tab.Type(SortIDX);
%            MyNewOrder(:,5) = num2cell(Tab.ContactSite(SortIDX));
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
plot([-0.0125 -0.0125],[1 size(PSTHall,1)],'r', 'LineWidth', 2);
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
% saveas(gcf, [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) '.png'])


%     % Plot the NeuronIDs
%     MyNewOrderStrings=cellfun(@num2str,MyNewOrder,'un',0);
%     for pp = 1:length(MyNewOrderStrings)
%         CurrentID=strjoin(MyNewOrderStrings(pp,:));
%         text(1,pp,[CurrentID(1:5) CurrentID(17:end)], 'FontSize', 8);
%     end
    switch Stim
        case {'BA_50_5Hz', 'BA_50_5Hz_Resting', 'BA_50_10Hz_Resting', 'BA_25_5Hz', 'BA_250_5Hz', 'BA_25_10Hz', 'BA_250_10Hz',...
                'shock_only', 'shock_only_Resting', 'BA_25_all_Resting', 'BA_25_5Hz_Resting', 'shocks', 'shock_motor_rest',...
                'BA_5Hz_motor_rest_first', 'BA_20Hz_motor_rest_first'}
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
%         switch Stim
%             case {'BA_25_5Hz', 'shock_only'}
%                 if tt == 1
%                     saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) Stim '_AVG_exc.png'])
%                 elseif tt == 2
%                     saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) Stim '_AVG_inh.png'])
%                 end
%             case {'TO_25_5Hz', 'shock_inh'}
%                 if tt == 1
%                     saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) Stim '_AVG_exc.png'])
%                 elseif tt == 2
%                     saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) Stim '_AVG_inh.png'])
%                 end
%         end
    end
end
end
        
