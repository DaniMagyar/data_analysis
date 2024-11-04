% Archive


clear all
[~, baseline_data] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_cond_first'},...
    'TTLselect', [1:5], 'psth_bin', 600, 'int', [-20 0]); % each ttl is a column

[~, response_data] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_shocks'},...
   'TTLselect', [1:5], 'psth_bin', 600, 'int', [-1 2]); % each ttl is a column


% wilcoxon
for qq = 1:size(baseline_data,1)
    respdat1 = reshape(response_data{qq,2}(1:5,:),1,[]); % edit here if binsize changed
    [pRanksum_short(qq,1), resRanksum_short(qq,1)] = ranksum(baseline_data{qq,1}(:), respdat1, 'alpha', 0.001); 
    if mean(baseline_data{qq,1}(:)) < mean(respdat1)
        wcx_comp_short{qq} = 'exc';
    elseif mean(baseline_data{qq,1}(:)) > mean(respdat1)
        wcx_comp_short{qq} = 'inh';
    else
        wcx_comp_short{qq} = 'non';
    end

    if wcx_comp_short{qq} == 'non'
        respdat2 = reshape(response_data{qq,2}(6:10,:),1,[]); % edit here if binsize changed
        [pRanksum_short(qq,1), resRanksum_short(qq,1)] = ranksum(baseline_data{qq,1}(:), respdat2, 'alpha', 0.001); 
        if mean(baseline_data{qq,1}(:)) < mean(respdat2)
            wcx_comp_short{qq} = 'exc';
        elseif mean(baseline_data{qq,1}(:)) > mean(respdat2)
            wcx_comp_short{qq} = 'inh';
        else
            wcx_comp_short{qq} = 'non';
        end
    end

    respdat3 = reshape(response_data{qq,2}(11:end,:),1,[]); % edit here if binsize changed
    [pRanksum_late(qq,1), resRanksum_late(qq,1)] = ranksum(baseline_data{qq,1}(:), respdat3, 'alpha', 0.001); 
    if mean(baseline_data{qq,1}(:)) < mean(respdat3)
        wcx_comp_late{qq} = 'exc';
    elseif mean(baseline_data{qq,1}(:)) > mean(respdat3)
        wcx_comp_late{qq} = 'inh';
    else
        wcx_comp_late{qq} = 'non';
    end

end

% %Calculate Zscore separately
% [~, baseline_data2] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_cond_first'},...
%     'TTLselect', [1:5], 'psth_bin', 12000, 'int', [-20 0]); % each ttl is a column
% 
% [~, response_data2] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_shocks'},...
%    'TTLselect', [1:5], 'psth_bin', 12000, 'int', [-1 2]); % each ttl is a column
% 
% for rr = 1:size(baseline_data2,1)
%     respdat4 = sum(response_data2{rr,2},2);
%     respdat4 = respdat4(1:end); % edit here if binsize changed
%     basedat4 = sum(baseline_data2{rr,1},2);
%     zscore_base(rr,:) = (basedat4 - mean(basedat4))/std(basedat4);
%     zscore_resp(rr,:) = (respdat4 - mean(basedat4))/std(basedat4);
%     if any(isnan(zscore_resp(rr,:))) || any(isnan(zscore_base(rr,:)))
%         zscore_comp_wcx(rr,1) = NaN;
%     else
%         zscore_comp_wcx(rr,1) = ranksum(zscore_base(rr,:), zscore_resp(rr,:), 'alpha', 0.05);
%     end
% end
% 
% [rows_exc,~] = find(zscore_resp>=3);
% [rows_inh,~] = find(zscore_resp<=-2);
% 
% rows_exc = sort(rows_exc);
% repeats_exc = diff(rows_exc) == 0; 
% rows_exc(repeats_exc) = [];
% 
% rows_inh = sort(rows_inh);
% repeats_inh = diff(rows_inh) == 0; 
% rows_inh(repeats_inh) = [];
% 
% rows_biphasic = intersect(rows_exc, rows_inh);
% cell_metrics.labels(rows_exc) = {'exc'};
% cell_metrics.labels(rows_inh) = {'inh'};
% cell_metrics.labels(rows_biphasic) = {'biphasic'};

%Calculate Zscore separately
[~, baseline_data2] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_cond_first'},...
    'TTLselect', [1:3], 'psth_bin', 3000, 'int', [-20 0]); % each ttl is a column

[~, response_data2] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_shocks'},...
   'TTLselect', [1:3], 'psth_bin', 3000, 'int', [-1 2.1]); % each ttl is a column

for rr = 1:size(baseline_data2,1)
    respdat100ms = sum(response_data2{rr,2}(1,:),2); % 0-100 ms
    respdat_other = sum(response_data2{rr,2}(2:end,:),2); % 100-2000ms
    basedat100ms = sum(baseline_data2{rr,1},2);
    zscore_base(rr,:) = (basedat100ms - mean(basedat100ms))/std(basedat100ms);
    zscore_100(rr,:) = (respdat100ms - mean(basedat100ms))/std(basedat100ms);
    basedat400ms = sum(reshape(basedat100ms, [4 50])',2);
    respdat400ms = sum(reshape(respdat_other, [4 5])',2);
    zscore_400(rr,:) = (respdat400ms - mean(basedat400ms))/std(basedat400ms);
    avgFR(rr,1) = sum(respdat400ms)/(2*5); % (2sec*5TTL)
    resp100spikenum(rr,:) = respdat100ms/5; %avg. spikenum in first 100 ms / 5TTL
end


[idx_exc_100,~] = find(zscore_100>=3);
[idx_inh_100,~] = find(zscore_100<=-2);
[idx_exc_400,~] = find(zscore_400>=3);
[idx_inh_400,~] = find(zscore_400<=-2);

% Excitatory response curation

idx_exc_100_filt = setdiff(idx_exc_100, find(resp100spikenum<1)); %filtering short exc. responses where spikenum<1/TTL
idx_exc_100_orig = setdiff(idx_exc_100, idx_exc_100_filt);

idx_exc_400_filt = setdiff(idx_exc_400, find(avgFR<1)); % Filtering late excitatory responses where FR<1Hz  after US
idx_exc_400_orig = setdiff(idx_exc_400, idx_exc_400_filt);

idx_exc_filt = sort([idx_exc_100_filt;idx_exc_400_filt]);
repeats_exc = diff(idx_exc_filt) == 0; 
idx_exc_filt(repeats_exc) = [];


% Inhibitory response curation

idx_inh = sort([idx_inh_100;idx_inh_400]);
repeats_inh = diff(idx_inh) == 0; 
idx_inh(repeats_inh) = [];

idx_biphasic = intersect(idx_exc_filt, idx_inh);

lookfor = [1 1];
trues = zscore_400<=-1;
for ss = 1:size(baseline_data2,1)
    out = strfind(trues(ss,:), lookfor);
    if any(out>=1)
        inh_400_lowZscore(ss,1) = 1;
    end
end
idx_inh_400_lowZscore = find(inh_400_lowZscore==1);
idx_inh_400_onlylowZscore = setdiff(idx_inh_400_lowZscore, idx_inh);





cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'neutral'};

cell_metrics.labels(idx_exc_100_orig) = {'exc_100_orig'};
cell_metrics.labels(idx_exc_400_orig) = {'exc_400_orig'};
cell_metrics.labels(idx_exc_filt) = {'exc_filt'};
cell_metrics.labels(idx_inh) = {'inh'};
cell_metrics.labels(idx_inh_400_onlylowZscore) = {'idx_inh_400_onlylowZscore'};
cell_metrics.labels(idx_biphasic) = {'biphasic'};



cell_metrics = CellExplorer('metrics',cell_metrics);









idx_wcx_short_exc = intersect(find(resRanksum_short == 1), find(contains(wcx_comp_short, 'exc')));
idx_wcx_short_inh = intersect(find(resRanksum_short == 1), find(contains(wcx_comp_short, 'inh')));
idx_wcx_late_exc = intersect(find(resRanksum_late == 1), find(contains(wcx_comp_late, 'exc')));
idx_wcx_late_inh = intersect(find(resRanksum_late == 1), find(contains(wcx_comp_late, 'inh')));

idx_wcx_both_exc = intersect(idx_wcx_short_exc, idx_wcx_late_exc);
idx_wcx_only_short_exc = setdiff(idx_wcx_short_exc,idx_wcx_both_exc);
idx_wcx_only_late_exc = setdiff(idx_wcx_late_exc,idx_wcx_both_exc);

idx_wcx_both_inh = intersect(idx_wcx_short_inh, idx_wcx_late_inh);
idx_wcx_only_short_inh = setdiff(idx_wcx_short_inh,idx_wcx_both_inh);
idx_wcx_only_late_inh = setdiff(idx_wcx_late_inh,idx_wcx_both_inh);

idx_wcx_short_exc_late_inh = intersect(idx_wcx_only_short_exc,idx_wcx_only_late_inh);
idx_wcx_short_inh_late_exc = intersect(idx_wcx_only_short_inh,idx_wcx_only_late_exc);



cell_metrics.labels(idx_wcx_only_short_exc) = {'only_short_exc'}; % this must be the first to not overwrite SELE/SELI
cell_metrics.labels(idx_wcx_only_short_inh) = {'only_short_inh'}; %this must be the first to not overwrite SILE/SILI

cell_metrics.labels(idx_wcx_only_late_exc) = {'only_late_exc'}; % this must be the second to not overwrite SELE/SELI
cell_metrics.labels(idx_wcx_only_late_inh) = {'only_late_inh'}; % this must be the second to not overwrite SELE/SELI

cell_metrics.labels(idx_wcx_both_exc) = {'short_exc_late_exc'};
cell_metrics.labels(idx_wcx_short_exc_late_inh) = {'short_exc_late_inh'};
cell_metrics.labels(idx_wcx_short_inh_late_exc) = {'short_inh_late_exc'};
cell_metrics.labels(idx_wcx_both_inh) = {'short_inh_late_inh'};


