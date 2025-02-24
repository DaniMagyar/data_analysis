function [responses_sound] =  BAfc_find_recall_v01(inputTTLcell)

[~, baseline_data] = BAfc_calc_Zscore_prenorm('Stims', inputTTLcell,...
    'psth_bin', 3000, 'int', [-20 0]); % each ttl is a column
[~, response_data] = BAfc_calc_Zscore_prenorm('Stims', inputTTLcell,...
    'psth_bin', 3000, 'int', [-1 1]); % each ttl is a column
cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'neutral'};

% zscore
for pp = 1:size(baseline_data,1)
    baseline_spx{pp} = sum(baseline_data{pp,1},2)/10;
    response_spx{pp} = sum(response_data{pp,2},2)/50;
    response_zscore{pp} = (response_spx{pp} - mean(baseline_spx{pp}))/std(baseline_spx{pp});
    meanFR(pp,1) = (sum(baseline_spx{pp}))/20;

    zscore_exc_high = 1.5;
    zscore_exc_low = 0.75;
    zscore_inh_high = -1;
    zscore_inh_low = -0.5;

    % find short responses + classification
    if response_zscore{pp}(1)>=zscore_exc_high
        short_window{pp} = 'excited';
    elseif response_zscore{pp}(1)>=zscore_exc_low & response_zscore{pp}(2)>=zscore_exc_low
        short_window{pp} = 'excited';
    elseif response_zscore{pp}(1)<=zscore_inh_high
        short_window{pp} = 'inhibited';
    elseif response_zscore{pp}(1)<=zscore_inh_low & response_zscore{pp}(2)<=zscore_inh_low
        short_window{pp} = 'inhibited';
    else
        short_window{pp} = 'neutral';
    end

    % find late responses + classification
    if response_zscore{pp}(2)>=zscore_exc_high
        late_window{pp} = 'excited';
    elseif response_zscore{pp}(2)>=zscore_exc_low & response_zscore{pp}(3)>=zscore_exc_low
        late_window{pp} = 'excited';
    elseif response_zscore{pp}(2)<=zscore_inh_high
        late_window{pp} = 'inhibited';
    elseif response_zscore{pp}(2)<=zscore_inh_low & response_zscore{pp}(3)<=zscore_inh_low
        late_window{pp} = 'inhibited';
    else
        late_window{pp} = 'neutral';
    end

end

% wilcoxon
for qq = 1:size(baseline_data,1)
    [pRanksum_short(qq,1), resRanksum_short(qq,1)] = ranksum(baseline_data{qq,1}(:), response_data{qq,2}(1,:), 'alpha', 0.01); 
    if mean(baseline_data{qq,1}(:)) < mean(response_data{qq,2}(1,:))
        wcx_comp_short{qq} = 'exc';
    elseif mean(baseline_data{qq,1}(:)) > mean(response_data{qq,2}(1,:))
        wcx_comp_short{qq} = 'inh';
    else
        wcx_comp_short{qq} = 'non';
    end
    
    [pRanksum_late(qq,1), resRanksum_late(qq,1)] = ranksum(baseline_data{qq,1}(:), response_data{qq,2}(2,:), 'alpha', 0.01); 
    if mean(baseline_data{qq,1}(:)) < mean(response_data{qq,2}(2,:))
        wcx_comp_late{qq} = 'exc';
    elseif mean(baseline_data{qq,1}(:)) > mean(response_data{qq,2}(2,:))
        wcx_comp_late{qq} = 'inh';
    else
        wcx_comp_late{qq} = 'non';
    end
end

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

idx_short_exc = find(contains(short_window,'excited'));
idx_short_inh = find(contains(short_window,'inhibited'));
idx_late_exc = find(contains(late_window,'excited'));
idx_late_inh = find(contains(late_window,'inhibited'));
idx_only_late_exc = setdiff(idx_late_exc, [idx_short_exc idx_short_inh]);
idx_only_late_inh = setdiff(idx_late_inh, [idx_short_exc idx_short_inh]);

idx_SELE = intersect(idx_short_exc, idx_late_exc); %short exc, late exc
idx_SELI = intersect(idx_short_exc, idx_late_inh); %short exc, late inh
idx_SILE = intersect(idx_short_inh, idx_late_exc); %short inh, late exc
idx_SILI = intersect(idx_short_inh, idx_late_inh); %short inh, late inh

% Detailed classification
% cell_metrics.labels(idx_wcx_only_short_exc) = {'only_short_exc'}; % this must be the first to not overwrite SELE/SELI
% cell_metrics.labels(idx_wcx_only_short_inh) = {'only_short_inh'}; %this must be the first to not overwrite SILE/SILI
% cell_metrics.labels(idx_wcx_only_late_inh) = {'only_late_inh'}; % this must be the second to not overwrite SELE/SELI
% cell_metrics.labels(idx_wcx_only_late_exc) = {'only_late_exc'}; % this must be the second to not overwrite SELE/SELI
% cell_metrics.labels(idx_wcx_both_exc) = {'short_exc_late_exc'};
% cell_metrics.labels(idx_wcx_short_exc_late_inh) = {'short_exc_late_inh'};
% cell_metrics.labels(idx_wcx_short_inh_late_exc) = {'short_inh_late_exc'};
% cell_metrics.labels(idx_wcx_both_inh) = {'short_inh_late_inh'};


% Simplified classification (exc, inh, neut)
cell_metrics.labels(idx_wcx_only_short_exc) = {'exc'}; % this must be the first to not overwrite SELE/SELI
cell_metrics.labels(idx_wcx_only_short_inh) = {'inh'}; %this must be the first to not overwrite SILE/SILI
cell_metrics.labels(idx_wcx_only_late_inh) = {'inh'}; % this must be the second to not overwrite SELE/SELI
cell_metrics.labels(idx_wcx_only_late_exc) = {'exc'}; % this must be the second to not overwrite SELE/SELI
cell_metrics.labels(idx_wcx_both_exc) = {'exc'};
cell_metrics.labels(idx_wcx_short_exc_late_inh) = {'exc'};
cell_metrics.labels(idx_wcx_short_inh_late_exc) = {'inh'};
cell_metrics.labels(idx_wcx_both_inh) = {'inh'};


responses_sound = cell_metrics.labels;