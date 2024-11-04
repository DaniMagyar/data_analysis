%Archive
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

cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'neutral'};

cell_metrics.labels(idx_exc_100_orig) = {'exc_100_orig'};
cell_metrics.labels(idx_exc_400_orig) = {'exc_400_orig'};
cell_metrics.labels(idx_exc_filt) = {'exc_filt'};
cell_metrics.labels(idx_inh) = {'inh'};
cell_metrics.labels(idx_biphasic) = {'biphasic'};



cell_metrics = CellExplorer('metrics',cell_metrics);
