function [idx] = BAfc_compare_habitvsrecall_Wcx_prenorm

psth_bin = 7500;
int = [-20 5];
fs = 30000;
trial_length = 5; %trial length is seconds
n_trial_bin = (trial_length*fs)/psth_bin;
n_baseline_bin = (abs(int(1))*fs)/psth_bin;


PSTHall_habit = BAfc_calc_wcx_NOTnorm('Stims', {'TTL_tone_habit_first'}, 'TTLselect', [6:10], 'psth_bin', psth_bin, 'int', int);
PSTHall_recall = BAfc_calc_wcx_NOTnorm('Stims', {'TTL_tone_recall_first'}, 'TTLselect', [1:5],'psth_bin', psth_bin, 'int', int);

for pp = 1:size(PSTHall_habit,1)
    if any(isnan(PSTHall_habit(pp,:)))
        pRanksum(pp) = 0;
        resRanksum(pp) = 0;   
    elseif any(isnan(PSTHall_recall(pp,:)))
        pRanksum(pp) = 0;
        resRanksum(pp) = 0;   
    else
        [pRanksum(pp), resRanksum(pp)] = ranksum(PSTHall_habit(pp,n_baseline_bin+1:n_baseline_bin+n_trial_bin),...
        PSTHall_recall(pp,n_baseline_bin+1:n_baseline_bin+n_trial_bin), 'alpha', 0.05);
    end
end

idx = find(resRanksum==1);
% cell_metrics = BAfc_load_neurons;
% cell_metrics.labels(1:end) = {'0'};
% cell_metrics.labels(idx) = {'1'};