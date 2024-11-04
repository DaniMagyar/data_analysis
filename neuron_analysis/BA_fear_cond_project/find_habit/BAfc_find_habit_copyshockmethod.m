function [responses_shock] =  BAfc_find_shocks_v01

TTLinclude = 3;

[~, baseline_data2] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_habit_first'},...
    'TTLselect', [1:TTLinclude], 'psth_bin', 3000, 'int', [-15 0]); % each ttl is a column
[~, response_data2] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_habit_all'},...
   'TTLselect', [1:TTLinclude*5], 'psth_bin', 3000, 'int', [0 1]); % each ttl is a column.
cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'neutral'};

for rr = 1:size(baseline_data2,1)
    respdat100ms = sum(response_data2{rr,2}(1,:),2)/(TTLinclude*5); % 0-100 ms , divided by 15
    respdat_other = sum(response_data2{rr,2}(2:end,:),2)/(TTLinclude*5); % 100-2000ms divided by 15
    basedat100ms = sum(baseline_data2{rr,1},2)/(TTLinclude); % divided by 3
    zscore_base(rr,:) = (basedat100ms - mean(basedat100ms))/std(basedat100ms);
    zscore_100(rr,:) = (respdat100ms - mean(basedat100ms))/std(basedat100ms);
    basedat400ms = sum(reshape(basedat100ms, [3 50])',2);
    respdat400ms = sum(reshape(respdat_other, [3 3])',2);
    zscore_400(rr,:) = (respdat400ms - mean(basedat400ms))/std(basedat400ms);
    avgFR(rr,1) = sum(respdat400ms)/(0.9*TTLinclude); % (0.9sec*TTL)
    resp100spikenum(rr,:) = respdat100ms/TTLinclude; %avg. spikenum in first 100 ms / 5TTL
end
[idx_exc_100,~] = find(zscore_100>=1.5);
[idx_inh_100,~] = find(zscore_100<=-1);
[idx_exc_400,~] = find(zscore_400>=1.5);
[idx_inh_400,~] = find(zscore_400<=-1);

% Excitatory response curation
idx_exc_100_filt = setdiff(idx_exc_100, find(resp100spikenum<1)); %filtering short exc. responses where spikenum<1/TTL
idx_exc_100_orig = setdiff(idx_exc_100, idx_exc_100_filt);

idx_exc_400_filt = setdiff(idx_exc_400, find(avgFR<1)); % Filtering late excitatory responses where FR<1Hz  after US
idx_exc_400_orig = setdiff(idx_exc_400, idx_exc_400_filt);

idx_exc_curated = sort([idx_exc_100_filt;idx_exc_400_filt]);
repeats_exc = diff(idx_exc_curated) == 0; 
idx_exc_curated(repeats_exc) = [];
cell_metrics.labels(idx_exc_curated) = {'exc'};

% Inhibitory response curation
idx_inh = sort([idx_inh_100;idx_inh_400]);
repeats_inh = diff(idx_inh) == 0; 
idx_inh(repeats_inh) = [];

neutrals = find(contains(cell_metrics.labels, 'neutral'));
lookfor = [1 1 1];
trues = zscore_400<=-0.5;
for ss = 1:size(baseline_data2,1)
    if ismember(ss, neutrals)
        out = strfind(trues(ss,:), lookfor);
        if any(out>=1)
            inh_400_lowZscore(ss,1) = 1;
        else 
            inh_400_lowZscore(ss,1) = 0;
        end
    else
        inh_400_lowZscore(ss,1) = 0;
    end
end
idx_inh_400_2zscore = find(inh_400_lowZscore==1);

idx_inh_curated = sort([idx_inh;idx_inh_400_2zscore]);
cell_metrics.labels(idx_inh_curated) = {'inh'};

% Sort biphasic neurons based on first 100 ms
idx_biphasic = intersect(idx_exc_curated, idx_inh_curated);
idx_exc_biph_100 = intersect(idx_exc_100_filt, idx_biphasic);
cell_metrics.labels(idx_exc_biph_100) = {'exc'};
idx_inh_biph_100 = intersect(idx_inh_100, idx_biphasic);
cell_metrics.labels(idx_inh_biph_100) = {'inh'};

idx_biphasic2 = setdiff(idx_biphasic, [idx_exc_biph_100; idx_inh_biph_100]);
cell_metrics.labels(idx_biphasic2) = {'biphasic400'};
neutrals_final = find(contains(cell_metrics.labels, 'neutral'));

responses_shock = cell_metrics.labels;