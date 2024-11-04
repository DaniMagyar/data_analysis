cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'neutral'};
% Using Gaussian filter with tone_habit_all\5, so each pulse is one trial. NOT GOOD
for ii = 1:numel(cell_metrics.cellID)
    disp(ii)
    TTL = cell_metrics.general.TTL_tone_habit_first{ii};
    AP = cell_metrics.spikes.times{ii};
    pre_time = 10; % in sec     
    post_time = 0;    
    bin_time = 0.001; %in seconds
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
    psth_spx_dani1(ii,:) =  vertcat(sum(preAP_bin,2));
    clearvars -except psth_spx_dani1 cell_metrics ii
end

for ii = 1:numel(cell_metrics.cellID)
    disp(ii)
    TTL = cell_metrics.general.TTL_tone_habit_all{ii};
    AP = cell_metrics.spikes.times{ii};
    pre_time = 0; % in sec     
    post_time = 1;    
    bin_time = 0.001; %in seconds
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
    psth_spx_dani2(ii,:) =  vertcat(sum(postAP_bin,2))/5;
    clearvars -except psth_spx_dani1 psth_spx_dani2 cell_metrics ii
end

psth_spx = horzcat(psth_spx_dani1, psth_spx_dani2);




sigma = 20; % standard deviation of the Gaussian kernel (in ms)
binWidth = 1; % bin width for the spike train (in ms)
edges = -3*sigma:binWidth:3*sigma; % create kernel edges
gaussKernel = exp(-edges.^2 / (2*sigma^2)); % Gaussian formula
gaussKernel = gaussKernel / sum(gaussKernel); % normalize kernel

for ii = 1: size(psth_spx,1)
    smoothed_spx(ii,:) = conv(psth_spx(ii,:), gaussKernel, 'same'); % convolve
    smoothed_spx_zscore(ii,:) = (smoothed_spx(ii,:) - mean(smoothed_spx(ii,1:10000)))/std(smoothed_spx(ii,1:10000));
    cell_metrics.general.smoothDataHabit{ii} = smoothed_spx_zscore(ii,:);
end

