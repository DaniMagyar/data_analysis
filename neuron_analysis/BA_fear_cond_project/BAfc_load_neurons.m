function [cell_metrics] = BAfc_load_neurons

basenames = {'temp_wh', 'temp_wh','temp_wh','temp_wh', 'temp_wh', 'temp_wh', 'temp_wh','temp_wh','temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD243_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD250_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD251_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD252_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD253_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD254_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD266_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD267_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD268_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD269_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD275_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD276_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD277_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD278_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);




cell_metrics.general.TTL_shocks = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks;
    cell_metrics.general.TTL_shocks(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.general.TTL_tone_habit_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_habit_first;
    cell_metrics.general.TTL_tone_habit_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_habit_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_habit_first;
    cell_metrics.general.TTL_noise_habit_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_tone_cond_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_cond_first;
    cell_metrics.general.TTL_tone_cond_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_cond_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_cond_first;
    cell_metrics.general.TTL_noise_cond_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_tone_recall_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_recall_first;
    cell_metrics.general.TTL_tone_recall_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_recall_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_recall_first;
    cell_metrics.general.TTL_noise_recall_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_tone_habit_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_habit_all;
    cell_metrics.general.TTL_tone_habit_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_habit_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_habit_all;
    cell_metrics.general.TTL_noise_habit_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_tone_cond_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_cond_all;
    cell_metrics.general.TTL_tone_cond_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_cond_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_cond_all;
    cell_metrics.general.TTL_noise_cond_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_tone_recall_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_recall_all;
    cell_metrics.general.TTL_tone_recall_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_recall_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_recall_all;
    cell_metrics.general.TTL_noise_recall_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
