function [idx] = BAfc_compare_habitvsrecall_Zscore_prenorm

psth_bin = 30000;
int = [-20 5];
fs = 30000;
trial_length = 5; %trial length is seconds
n_trial_bin = (trial_length*fs)/psth_bin;
n_baseline_bin = (abs(int(1))*fs)/psth_bin;


PSTHall_habit = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_habit_first'}, 'TTLselect', [6:10], 'psth_bin', psth_bin, 'int', int);
PSTHall_recall = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_recall_first'}, 'TTLselect', [1:5],'psth_bin', psth_bin, 'int', int);
PSTHall_diff = PSTHall_recall-PSTHall_habit;
significant = any(PSTHall_diff(:,n_baseline_bin+1:n_baseline_bin+n_trial_bin)>=3,2);
idx = find(significant==1);
% cell_metrics = BAfc_load_neurons;
% cell_metrics.labels(1:end) = {'0'};
% cell_metrics.labels(idx) = {'1'};