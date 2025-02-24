function [cell_metrics] = BAfc_load_neurons(varargin)

% Default params
prs = inputParser;
addRequired(prs,'experiment',@ischar)
parse(prs, varargin{:})
g = prs.Results;

switch g.experiment
    case 'BAfc'
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

    case 'NP_BAfc'
        basenames = {'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh'};
        basepaths = {...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD288_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD289_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD290_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD291_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD293_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD294_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD295_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD296_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD297_kilosort\kilosort25preprocess'};
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

    case 'NP_BAfc_triptest'
        basenames = {'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh'};
        basepaths = {...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD292_002_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD293_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD294_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD295_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD296_kilosort\kilosort25preprocess', ...
            'C:\Users\dmagyar\Desktop\BA_fear_cond\MD297_kilosort\kilosort25preprocess'};
        cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

              
        cell_metrics.general.TTL_triptest_shocks_only = {};
        for ii = 1:max(cell_metrics.batchIDs)
            allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
            TTL = allTTL.triptest_shocks_only;
            cell_metrics.general.TTL_triptest_shocks_only(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
        end
        clear TTL AP ii allTTL


        cell_metrics.general.TTL_triptest_sound_only = {};
        for ii = 1:max(cell_metrics.batchIDs)
            allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
            TTL = allTTL.triptest_sound_only;
            cell_metrics.general.TTL_triptest_sound_only(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
        end
        clear TTL AP ii allTTL


        cell_metrics.general.TTL_triptest_both = {};
        for ii = 1:max(cell_metrics.batchIDs)
            allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
            TTL = allTTL.triptest_both;
            cell_metrics.general.TTL_triptest_both(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
        end
        clear TTL AP ii allTTL


end
