function Zscore_analysis_Dani(File, nucleus, Stim, varargin)
% Default params ----------------------------------------------------------
fs       = 30000; % Firing rate
int      = [-5 5]; % psth interval, even numbers are better for plot
norm     = 1; % Normalise data
psth_bin = 6000; % 600 = 20ms
testBins = 10; % bins included for Z-score ordering 
Wcx_win = [-5 2]; % Wilcoxon window in seconds
average = 1; % plot averages, 1 if yes
mainFolder = 'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december';

% lines to change BA_25 -> 250 : line138-130, line189-190 

% User defined params -----------------------------------------------------
if nargin
    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'fs'
                fs       = varargin{ii+1};
            case 'interval'
                int      = varargin{ii+1};
            case 'numPCA'
                PCA_num  = varargin{ii+1};
            case 'numClusers'
                clustnum = varargin{ii+1};
            case 'normality'
                norm     = varargin{ii+1};
            case 'method'
                LKmethod = varargin{ii+1};
            case 'metric'
                LKmetric = varargin{ii+1};
            case 'bin'
                psth_bin = varargin{ii+1};
        end
    end
end

% Select neurons assigned to the selected NUCLEUS and Set the folders -----
Tab = File;
for ii = size(File,1):-1:1
    if startsWith(File.Location(ii), nucleus) == 0
        Tab(ii,:) = [];
    end
end


folder     = Tab.Recording;
nNeurons   = size(folder,1);

%% Generate the PSTHall matrix --------------------------------------------
for kk = 1:nNeurons
    % Go to folder
    cd([mainFolder filesep folder{kk,:}]);
    
    % Load TTLs & Choose variable to use as stimulus
    load('TTLs.mat'); %#ok<LOAD>
    switch Stim
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
            ttl = BA_250_5Hz;
%             ttl = BA_250_5Hz(1:125);
%             ttl = [ttl; BA_250_5Hz(1:10:end)];
%             ttl = [ttl; BA_250_5Hz(2:10:end)];
%             ttl = [ttl; BA_250_5Hz(3:10:end)];
%             ttl = [ttl; BA_250_5Hz(4:10:end)];
%             ttl = sort(ttl);
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
    end
    % Select the neuron
    NeuronID = dir(['GR',Tab.Group{kk},'_',num2str(Tab.Neuron(kk)),...
        '.mat']);
    NeuronID = NeuronID.name;    
    load(NeuronID,'TS')
    TT = TS(:,1)/10000;

    % PSTH matrix
    bin_time = psth_bin/fs;     
    pre_time = abs(int(1));      
    post_time = int(2);      
    AP = TT;
    TTL = ttl;
    for jj = 1:numel(TTL) %Each TTL is a a column. 
         preAP{:,jj} = AP(AP>=(TTL(jj)-pre_time) & AP<TTL(jj)); % spikes before each TTL separately
         postAP{:,jj} = AP(AP>TTL(jj) & AP<(TTL(jj)+post_time)); % spikes after each TTL separately
         %postAP{:,jj} = AP(AP>(TTL(jj)+0.005) & AP<(TTL(jj)+0.005+post_time)); %BA_250
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

     %% paired Wilcoxon signed rank test: equal length of pre and post data. 
     % Each post point cointain all spikes during the whole illumination period. (Like one 2 secs long bin)
     for wcx = 1:numel(TTL) %Each TTL is a a column. 
         preAP_wcx{:,wcx} = AP(AP>=(TTL(wcx)-abs(Wcx_win(1))) & AP<TTL(wcx)); % spikes before each TTL separately
         postAP_wcx{:,wcx} = AP(AP>TTL(wcx) & AP<(TTL(wcx)+Wcx_win(2))); % spikes after each TTL separately
         %postAP{:,wcx} = AP(AP>(TTL(wcx)+0.005) & AP<(TTL(wcx)+0.005+Wcx_win)); %BA_250
     end
     preAP_wcx_num = cellfun(@numel, preAP_wcx);
     postAP_wcx_num = cellfun(@numel, postAP_wcx);
     preAP_wcx_freq = preAP_wcx_num/abs(Wcx_win(1));
     postAP_wcx_freq = postAP_wcx_num/Wcx_win(2);
     [~,h] = signrank(preAP_wcx_freq, postAP_wcx_freq);
     Wilcoxon_results{kk,1} = NeuronID;
     Wilcoxon_results{kk,2} = num2str(h);

    % remove laser artefact bin (psth_spx length match PSTHall length)
        %psth_spx(abs(int(1)*fs/psth_bin)+1)=[]; %to remove first bin after TTL
        %psth_spx_dani(end)=[]; % to remove last bin  
    
        if norm == 1
            newpsth = zscore(psth_spx_dani);
        else
            newpsth = psth_spx_dani;
        end
        PSTHall(kk,:) = newpsth';
    end

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
    
    testWindow_firstBin = abs(int(1))*fs/psth_bin+1;
    testWindow_lastBin = testWindow_firstBin + (testBins-1);    
    testWindow = PSTHall(:, testWindow_firstBin:testWindow_lastBin);
    testMean = mean(testWindow,2);

    %% Calculating plotting order based on Z-score change
    switch Stim(1:2)
        case 'BA'
            [~,SortIDX] = sort(testMean, 'descend');
            MyNewOrder = table;
            MyNewOrder.Recording = Tab.Recording(SortIDX);
            MyNewOrder.Group = Tab.Group(SortIDX);
            MyNewOrder.Neuron = Tab.Neuron(SortIDX);
            MyNewOrder.Type = Tab.Type(SortIDX);
%            MyNewOrder(:,5) = num2cell(Tab.ContactSite(SortIDX));
        case 'TO'
            load ([mainFolder '\BAparams.mat'], 'SortIDX','MyNewOrder')
    end
    %% Common part
    
    figure; 
    subplot(1,4,1:2)
    imagesc(time/fs, 1:size(PSTHall,1), PSTHall(SortIDX,:)); % plotting sorted psth matrix 
    clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
    colormap(mycolormap); 
    hold on;
    plot([-0.125 -0.125],[1 size(PSTHall,1)],'r', 'LineWidth', 2);
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
            if testMeanSorted(qq) < 0
                response_dir(qq) = -1;
            elseif testMeanSorted(qq) >= 0
                response_dir(qq) = 1;
            end
        else
            response_dir(qq) = 0;
        end
    end    

    annotation('textbox', [0.55, 0.9, 1, 0], 'string',...
        ['Significantly modulated units:' num2str(numel(significant_idx))], 'LineStyle','none', 'FontSize', 35)
    annotation('textbox', [0.55, 0.7, 1, 0], 'string',...
        ['Excited: ' num2str(numel(find(response_dir == 1))) '  ('...
        num2str(numel(find(response_dir == 1))/length(testMean)*100) '%)'], 'LineStyle','none', 'FontSize', 35)
    annotation('textbox', [0.55, 0.5, 1, 0], 'string',...
        ['Inhibited: ' num2str(numel(find(response_dir == -1))) '  ('...
        num2str(numel(find(response_dir == -1))/length(testMean)*100) '%)'], 'LineStyle','none', 'FontSize', 35)

    set(gcf,'Position',[500 50 1600 1300])
    saveas(gcf, [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) '.png'])
   

%     % Plot the NeuronIDs
%     MyNewOrderStrings=cellfun(@num2str,MyNewOrder,'un',0);
%     for pp = 1:length(MyNewOrderStrings)
%         CurrentID=strjoin(MyNewOrderStrings(pp,:));
%         text(1,pp,[CurrentID(1:5) CurrentID(17:end)], 'FontSize', 8);
%     end
        switch Stim(1:2)
            case 'BA'
                save ([mainFolder '\BAparams.mat'], 'SortIDX','MyNewOrder', 'response_dir')
            case 'TO'
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
            if Stim(1:2) == 'BA'
                if tt == 1
                    saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) '_BA_AVG_exc.png'])
                elseif tt == 2
                    saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) '_BA_AVG_inh.png'])
                end
            elseif Stim(1:2) == 'TO'
                if tt == 1
                    saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) '_TO_AVG_exc.png'])
                elseif tt == 2
                    saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) '_TO_AVG_inh.png'])
                end
            end
        end
    end
end


