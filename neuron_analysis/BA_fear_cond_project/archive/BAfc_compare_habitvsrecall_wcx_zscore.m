psth_bin = 15000;
int = [-20 5];
fs = 30000;
trial_length = 5; %trial length is seconds
n_trial_bin = (trial_length*fs)/psth_bin;
n_baseline_bin = (abs(int(1))*fs)/psth_bin;


idx1 = BAfc_compare_habitvsrecall_Wcx_prenorm;
%idx2 = BAfc_compare_habitvsrecall_Zscore_prenorm;
PSTHall_recall = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_recall_first'}, 'TTLselect', [1:5],'psth_bin', psth_bin, 'int', int);
significant = any(PSTHall_recall(:,n_baseline_bin+1:n_baseline_bin+n_trial_bin)>=3,2);
idx2 = find(significant==1);

for ii = idx2'
    test_period = PSTHall_recall(ii,n_baseline_bin+1:n_baseline_bin+n_trial_bin);
    num = numel(test_period(test_period>=3));
    if num >= 2
        idxpass(find(idx2==ii)) = 1;
    else
        idxpass(find(idx2==ii)) = 0;
    end
end
idx2(find(idxpass==0)) = [];
idx = intersect(idx1, idx2);
cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'0'};
cell_metrics.labels(idx) = {'1'};