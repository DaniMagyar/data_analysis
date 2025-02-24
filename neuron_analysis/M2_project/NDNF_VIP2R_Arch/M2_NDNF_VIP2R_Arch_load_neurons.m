function [cell_metrics] = M2_NDNF_VIP2R_Arch_load_neurons

basenames = {'temp_wh', 'temp_wh','temp_wh','temp_wh','temp_wh','temp_wh', 'temp_wh','temp_wh','temp_wh','temp_wh', 'temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD281_001_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD281_002_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD281_003_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD281_004_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD281_005_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD282_001_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD282_002_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD283_002_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_NDNF_VIP2R_Arch\MD283_003_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);



cell_metrics.general.TTL_shocks_only = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks_only;
    cell_metrics.general.TTL_shocks_only(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_shocks_light = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks_light;
    cell_metrics.general.TTL_shocks_light(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL


