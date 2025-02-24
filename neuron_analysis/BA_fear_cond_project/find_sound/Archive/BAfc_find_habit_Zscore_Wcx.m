clear all
psth_bin = 7500; % 7500  = 250ms bin time
int = [-20 5];
fs = 30000;
trial_length = 5; %trial length is seconds
n_trial_bin = (trial_length*fs)/psth_bin;
n_baseline_bin = (abs(int(1))*fs)/psth_bin;


[PSTHall_habit, psth_spx_bin] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_habit_first'},...
    'TTLselect', [1:3], 'psth_bin', psth_bin, 'int', int);

% Find excitation
for pp = 1:size(PSTHall_habit,1)
    data_1 = psth_spx_bin{pp,1};
    %data_2 = psth_spx_bin{pp,2};
    % only the first bin after stim (500ms binsize)
    %data_2 = vertcat(psth_spx_bin{pp,2}(1:(fs/psth_bin):end,:));
    % first 2 bin after stim (250 ms binsize)
    data_2 = vertcat(psth_spx_bin{pp,2}(1:(fs/psth_bin):end,:), psth_spx_bin{pp,2}(2:(fs/psth_bin):end,:));    
    [pRanksum(pp,1), resRanksum(pp,1)] = ranksum(data_1(:), data_2(:));
    post_zscore(pp,:) = (sum(data_2,2)-(mean(sum(data_1,2))))/std(sum(data_1,2));
    if any(post_zscore(pp,:)>=3)
        Significant(pp,1) = 1;
    elseif any(post_zscore(pp,:)>=2)
        trues = post_zscore(pp,:)>=2;
        lookfor = [1 1];
        out = strfind(trues, lookfor);
        out_odd = mod(out,2);
        if any(out_odd == 1)
            Significant(pp,1) = 1;
        end
    else
        Significant(pp,1) = 0;
    end
end


[PSTHall_habit, psth_spx_bin] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_habit_first'},...
    'TTLselect', [1:3], 'psth_bin', 7500, 'int', [-20 5]);

% Find inhibition
for pp = 1:size(PSTHall_habit,1)
    data_1 = psth_spx_bin{pp,1};
    data_2 = vertcat(psth_spx_bin{pp,2}(1:(fs/psth_bin):end,:), psth_spx_bin{pp,2}(2:(fs/psth_bin):end,:)); %first two bins
    data_3 = psth_spx_bin{pp,2}; % all bins
    [pRanksum(pp,2), resRanksum(pp,2)] = ranksum(data_1(:), data_3(:,1)); % Wcx calculated from first trial, but all bins
    post_zscore2(pp,:) = (sum(data_2(:,1),2)-(mean(sum(data_1,2))))/std(sum(data_1,2)); % Zscore caluclated from first two bins after cs 
    post_zscore3(pp,:) = (sum(data_3(:,1),2)-(mean(sum(data_1,2))))/std(sum(data_1,2)); %all bins zscore
%     [pRanksum(pp,2), resRanksum(pp,2)] = ranksum(data_1(:), data_2(:)); % all trial
%     post_zscore2(pp,:) = (sum(data_2(:),2)-(mean(sum(data_1,2))))/std(sum(data_1,2)); 
    preAP_bin = psth_spx_bin{pp,1};
    postAP_bin = psth_spx_bin{pp,2}(:,1);

    psth_spx_zscore(pp,:) = (postAP_bin - mean(preAP_bin(:)))/std(preAP_bin(:));


    if any(post_zscore2(pp,:)<=-2)
        Significant(pp,2) = -1;
    elseif any(post_zscore2(pp,:)<=-1)
        trues = post_zscore2(pp,:)<=-1;
        lookfor = [1 1];
        out = strfind(trues, lookfor);
        out_odd = mod(out,2);
        if any(out_odd == 1)
            Significant(pp,2) = -1;
        end
    else
        Significant(pp,2) = 0;
    end    
end

mean_post_zscore2 = mean(post_zscore2,2);
mean_post_zscore3 = mean(post_zscore3,2);
mean_post_zscore4 = mean(psth_spx_zscore,2);






idx = find(resRanksum(:,1)==1)';
idx2 = find(Significant(:,1)==1);

idx4 = intersect(idx, idx2); % excited neurons
idx5 = setdiff(idx, idx2); % only Wcx excit
idx6 = setdiff(idx2, idx); % only Zscore excit

idx7 = find(resRanksum(:,2)==1)';
idx8 = find(Significant(:,2)==-1); 
idx9 = intersect(idx7, idx8); % inhibited neurons
idx10 = setdiff(idx7, idx8); % only Wcx inhib
idx11 = setdiff(idx8, idx7); % only Zscore inhib
idx12 = find(mean_post_zscore2<=-2);
idx13 = find(mean_post_zscore4<=-0.2);

idx14 = vertcat(idx12, idx13);
idx15 = intersect(idx4, idx13); % both responses
idx16 = intersect(idx7, idx13); % inhibited and wcx diff.

cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'0'};
cell_metrics.labels(idx4) = {'excited_Both'};
%cell_metrics.labels(idx5) = {'excited_Wcx'};
cell_metrics.labels(idx6) = {'excited_Zscore'};
%cell_metrics.labels(idx9) = {'inhibited_Both'};
%cell_metrics.labels(idx10) = {'inhibited_Wcx'};
cell_metrics.labels(idx16) = {'inhibited_Zscore'};
cell_metrics.labels(idx15) = {'Biphasic'};


