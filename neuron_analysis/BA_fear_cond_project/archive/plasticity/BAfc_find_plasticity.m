function [pRanksum, resRanksum] = BAfc_find_plasticity

cell_metrics = BAfc_load_neurons;
% Load Habit psth_spx
for ii = 1:numel(cell_metrics.cellID)
    TTL = cell_metrics.general.TTL_tone_habit_all{ii};
    AP = cell_metrics.spikes.times{ii};
    pre_time = 1; % in sec     
    post_time = 5;    
    bin_time = 0.05; %in seconds
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
    post_spx_habit(ii,:) = sum(postAP_bin,2);
    clearvars -except  cell_metrics ii post_spx_habit post_spx_recall
end

% Load Habit psth_spx
for ii = 1:numel(cell_metrics.cellID)
    TTL = cell_metrics.general.TTL_tone_recall_all{ii};
    AP = cell_metrics.spikes.times{ii};
    pre_time = 1; % in sec     
    post_time = 5;    
    bin_time = 0.05; %in seconds
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
    post_spx_recall(ii,:) = sum(postAP_bin,2);
    clearvars -except  cell_metrics ii post_spx_habit post_spx_recall 
end




window200_habit = horzcat(...
    post_spx_habit(:,1:4), ...
    post_spx_habit(:,21:24), ...
    post_spx_habit(:,41:44),...
    post_spx_habit(:,61:64),...
    post_spx_habit(:,81:84));

window200_recall = horzcat(...
    post_spx_recall(:,1:4), ...
    post_spx_recall(:,21:24), ...
    post_spx_recall(:,41:44),...
    post_spx_recall(:,61:64),...
    post_spx_recall(:,81:84));


for ii = 1:numel(cell_metrics.cellID)
    [pRanksum(ii), resRanksum(ii)] = signrank(window200_habit(ii,:), window200_recall(ii,:), 'alpha', 0.05);
end