function [responses_shock] =  BAfc_find_shocks_v01      %10 ms bintime

TTLinclude = 10;

[~, baseline_data2] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_habit_first'},...
    'TTLselect', [1:TTLinclude], 'psth_bin', 3000, 'int', [-20 0]); % each ttl is a column
[~, response_data2] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_habit_all'},...
   'TTLselect', [1:TTLinclude*5], 'psth_bin', 3000, 'int', [0 1]); % each ttl is a column.
cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'neutral'};

for rr = 1:size(baseline_data2,1)
    respdat100ms = sum(response_data2{rr,2}(1,:),2)/(TTLinclude*5); % 0-100 ms , divided by ttl*5
    basedat100ms = sum(baseline_data2{rr,1},2)/(TTLinclude); % divided by 3
    zscore_base100ms(rr,:) = (basedat100ms - mean(basedat100ms))/std(basedat100ms);
    zscore_100(rr,:) = (respdat100ms - mean(basedat100ms))/std(basedat100ms);
    resp100spikenum(rr,:) = respdat100ms; %avg. spikenum in first 100 ms !!! ALREADY DIVIDED BY 'TTLinclude'
end
[idx_exc_100,~] = find(zscore_100>=3);
[idx_inh_100,~] = find(zscore_100<=-2);
% Excitatory response curation
idx_exc_100_filt = setdiff(idx_exc_100, find(resp100spikenum<0.5)); %filtering short exc. responses 
idx_exc_100_orig = setdiff(idx_exc_100, idx_exc_100_filt);
idx_exc_100_curated = sort([idx_exc_100_filt]);

repeats_exc_100 = diff(idx_exc_100_curated) == 0; 
idx_exc_100_curated(repeats_exc_100) = [];
cell_metrics.labels(idx_exc_100_curated) = {'exc_100'};
cell_metrics.labels(idx_exc_100_orig) = {'exc_100_lowFR'};





[~, baseline_data3] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_habit_first'},...
    'TTLselect', [1:TTLinclude], 'psth_bin', 30000, 'int', [-20 0]); % each ttl is a column
[~, response_data3] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_habit_first'},...
   'TTLselect', [1:TTLinclude], 'psth_bin', 30000, 'int', [0 5]); % each ttl is a column.

neutrals = find(contains(cell_metrics.labels, 'neutral'));
for rr = 1:size(baseline_data3,1)
    if ismember(rr, neutrals)
        respdat400ms = sum(response_data3{rr,2},2)/(TTLinclude); 
        basedat400ms = sum(baseline_data3{rr,1},2)/(TTLinclude); 
        zscore_base400ms(rr,:) = (basedat400ms - mean(basedat400ms))/std(basedat400ms);
        zscore_400(rr,:) = (respdat400ms - mean(basedat400ms))/std(basedat400ms);
        resp400spikenum(rr,:) = mean(respdat400ms); %mean FR in 400 ms bins, ALREADY divided by ttl
    end
end

[idx_exc_400,~] = find(zscore_400>=3);
[idx_inh_400,~] = find(zscore_400<=-2);
% Excitatory response curation
idx_exc_400_filt = setdiff(idx_exc_400, find(resp400spikenum<1)); %filtering short exc. responses 
idx_exc_400_orig = setdiff(idx_exc_400, idx_exc_400_filt);
idx_exc_400_curated = sort([idx_exc_400_filt]);

repeats_exc_400 = diff(idx_exc_400_curated) == 0; 
idx_exc_400_curated(repeats_exc_400) = [];
cell_metrics.labels(idx_exc_400_curated) = {'exc_400'};
cell_metrics.labels(idx_exc_400_orig) = {'exc_400_lowFR'};





cell_metrics = CellExplorer('metrics',cell_metrics);