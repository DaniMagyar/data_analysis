function [responses_shock, responses_shock_class] =  BAfc_find_shocks_v01(TTLinclude)

[~, baseline_data2] = BAfc_calc_Zscore_v01('Stims', {'TTL_tone_cond_first'},...
    'TTLselect', TTLinclude, 'psth_bin', 3000, 'int', [-20 0]); % each ttl is a column
[~, response_data2] = BAfc_calc_Zscore_v01('Stims', {'TTL_shocks'},...
   'TTLselect', TTLinclude, 'psth_bin', 3000, 'int', [-1 2.1]); % each ttl is a column.
cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'neutral'};

for rr = 1:size(baseline_data2,1)
    respdat100ms = sum(response_data2{rr,2}(1,:),2); % 0-100 ms
    respdat_other = sum(response_data2{rr,2}(2:end,:),2); % 100-2000ms
    basedat100ms = sum(baseline_data2{rr,1},2);
    zscore_base(rr,:) = (basedat100ms - mean(basedat100ms))/std(basedat100ms);
    zscore_100(rr,:) = (respdat100ms - mean(basedat100ms))/std(basedat100ms);
    basedat400ms = sum(reshape(basedat100ms, [4 50])',2);
    respdat400ms = sum(reshape(respdat_other, [4 5])',2);
    zscore_400(rr,:) = (respdat400ms - mean(basedat400ms))/std(basedat400ms);
    avgFR(rr,1) = sum(respdat400ms)/(2*numel(TTLinclude)); % (2sec*5TTL)
    resp100spikenum(rr,:) = respdat100ms/numel(TTLinclude); %avg. spikenum in first 100 ms / numTTL
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

idx_exc_curated = sort([idx_exc_100_filt;idx_exc_400_filt]);
repeats_exc = diff(idx_exc_curated) == 0; 
idx_exc_curated(repeats_exc) = [];
cell_metrics.labels(idx_exc_curated) = {'exc'};

% Inhibitory response curation
idx_inh = sort([idx_inh_100;idx_inh_400]);
repeats_inh = diff(idx_inh) == 0; 
idx_inh(repeats_inh) = [];

neutrals = find(contains(cell_metrics.labels, 'neutral'));
lookfor = [1 1 1 1];
trues = zscore_400<=-1;
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
%cell_metrics.labels(idx_biphasic2) = {'biphasic400'};
neutrals_final = find(contains(cell_metrics.labels, 'neutral'));

responses_shock = cell_metrics.labels;

% Shock responses with classification
cell_metrics.shockRespClass(1:numel(responses_shock)) = {'neutral'};


cell_metrics.shockRespClass(idx_exc_curated) = {'excLate'}; %all excited receives "excLate" 
cell_metrics.shockRespClass(idx_exc_100_filt) = {'excEarly'}; % overwrite excLate values if also receives early excit
cell_metrics.shockRespClass(idx_inh_curated) = {'inhLate'}; %all excited receives "excLate" 
cell_metrics.shockRespClass(idx_inh_100) = {'inhEarly'}; % overwrite excLate values if also receives early excit

responses_shock_class = cell_metrics.shockRespClass;
