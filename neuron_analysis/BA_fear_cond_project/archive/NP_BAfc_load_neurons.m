function [cell_metrics] = NP_BAfc_load_neurons

basenames = {'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD288_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD289_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD290_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD291_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

cell_metrics.general.TTL_shocks = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks;
    cell_metrics.general.TTL_shocks(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.general.TTL_shocks_nonpredicted = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks(1:20);
    cell_metrics.general.TTL_shocks_nonpredicted(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.general.TTL_shocks_predicted = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks(21:40);
    cell_metrics.general.TTL_shocks_predicted(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.general.TTL_tone_habit_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_habit_first;
    cell_metrics.general.TTL_tone_habit_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.general.TTL_tone_cond_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_cond_first;
    cell_metrics.general.TTL_tone_cond_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.general.TTL_tone_recall_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_recall_first;
    cell_metrics.general.TTL_tone_recall_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.general.TTL_tone_habit_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_habit_all;
    cell_metrics.general.TTL_tone_habit_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.general.TTL_tone_cond_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_cond_all;
    cell_metrics.general.TTL_tone_cond_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.general.TTL_tone_recall_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_recall_all;
    cell_metrics.general.TTL_tone_recall_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

