
psth_bin = 7500; % 7500  = 250ms bin time
int = [-20 5];
fs = 30000;
trial_length = 5; %trial length is seconds
n_trial_bin = (trial_length*fs)/psth_bin;
n_baseline_bin = (abs(int(1))*fs)/psth_bin;


PSTHall_habit = BAfc_calc_Zscore_prenorm('Stims', {'TTL_tone_habit_first'},...
    'TTLselect', [1:3], 'psth_bin', psth_bin, 'int', int);

for pp = 1:size(PSTHall_habit,1)
    if any(horzcat(PSTHall_habit(pp,n_baseline_bin+1:4:n_baseline_bin+n_trial_bin), ...
            PSTHall_habit(pp,n_baseline_bin+2:4:n_baseline_bin+n_trial_bin))>=3,2);
        Significant(pp) = 1;
    else 
        Significant(pp) = 0;
    end
end

idx = find(Significant==1);
cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'0'};
cell_metrics.labels(idx) = {'1'};