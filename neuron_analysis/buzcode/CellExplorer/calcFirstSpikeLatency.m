function calcFirstSpikeLatency 

search_time = 0.25; % in seconds

basenames = {'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD138_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD139_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD140_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD141_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD142_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD143_kilosort\kilosort3preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);


for ii = 1:length(cell_metrics.cellID)
    batch  = cell_metrics.batchIDs(ii);
    allTTL = load([cell_metrics.general.basepaths{1} '\TTLsKS.mat']);
    TTL = allTTL.shock_only; % select TTL 
    AP = cell_metrics.spikes.times{ii};
    for jj = 1:numel(TTL) %Each TTL is a a column. 
         postAP{jj} = min(AP(AP>TTL(jj) & AP<(TTL(jj)+search_time))); % spikes after each TTL separately
         postAP_norm{jj} = postAP{jj}-TTL(jj);
    end
    boxchart(cell2mat(postAP_norm(find(~cellfun(@isempty,postAP_norm)))), 'orientation', 'horizontal')
end

disp('done')