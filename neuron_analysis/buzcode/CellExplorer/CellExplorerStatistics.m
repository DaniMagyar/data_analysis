function CellExplorerStatistics

load('temp_wh.cell_metrics.cellinfo.mat')
load('TTLsKS.mat', 'BA_25_5Hz')
% Default params ----------------------------------------------------------
fs       = 30000; % Firing rate
int      = [-5 5]; % psth interval, even numbers are better for plot
norm     = 1; % Normalise data
psth_bin = 6000; % 600 = 20ms
testBins = 10; % bins included for Z-score ordering 
Wcx_win = [-5 2]; % Wilcoxon window in seconds

ttl = BA_25_5Hz+8;%%%%%%%%%%%%%

for ii = 1: length(cell_metrics.cellID)
    TT = cell_metrics.spikes.times{ii};

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
     [p,h] = signrank(preAP_wcx_freq, postAP_wcx_freq);
     cell_metrics.Wilcoxon_resH_TO(ii) = double(h);  %%%%%%%%%%
     cell_metrics.Wilcoxon_resP_TO(ii) = p; %%%%%%%%%%%%

     if norm == 1
        newpsth = zscore(psth_spx_dani);
     else
        newpsth = psth_spx_dani;
     end
     PSTHcurr(ii,:) = newpsth';
end

testWindow_firstBin = abs(int(1))*fs/psth_bin+1;
testWindow_lastBin = testWindow_firstBin + (testBins-1);    
testWindow = PSTHcurr(:, testWindow_firstBin:testWindow_lastBin);
cell_metrics.testMeanTO = mean(testWindow,2)'; %%%%%%%%%%%%%%

save('temp_wh.cell_metrics.cellinfo.mat', 'cell_metrics')


