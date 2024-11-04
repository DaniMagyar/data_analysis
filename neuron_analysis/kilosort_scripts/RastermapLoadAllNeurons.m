function RastermapLoadAllNeurons(varargin)
% Default params
prs =  inputParser;
addParameter(prs,'mainFolder','C:\Users\dmagyar\Desktop\PFC_shock_run_rest',@ischar) % Path of main folder. (e.g. 'C:\Users\dmagyar\Desktop\M2_shock_response')
addParameter(prs,'Stim','none',@ischar) % TTLs to use from TTLsKS.mat. (e.g. 'shocks')
addParameter(prs,'loadBatch','none',@ischar) % Use batch file instead of original data, default 'none'. (e.g. 'cell_metrics_batch.mat')
addParameter(prs,'Recordings',@iscell) % Use original files; list of included recordings in cell format. (e.g. {'MD127', 'MD128', 'MD129'}).
addParameter(prs,'selVariable','none',@ischar) % Variable used for neuron selection from cell_metrics. (e.g. 'labels')
addParameter(prs,'selValue','none',@ischar) % Value used for neuron selection from 'selectVariable'. (e.g. '2')
addParameter(prs,'neuronType','none',@ischar) % Options: 'Narrow Interneuron', 'Wide Interneuron', 'Pyramidal Cell'. Default 'none' loads all types. 
addParameter(prs,'plotType','Heatmap',@ischar) % Options: 'Heatmap', 'Lineplot'
addParameter(prs,'Sorted',0,@isnumeric) % If 1, Load 'SortIDX' from 'BAparams.mat'. Default 0.
addParameter(prs,'fs',30000,@isnumeric) % Sampling rate, default 30000.
addParameter(prs,'psth_bin',6000,@isnumeric) % Binning window calculated from fs. (e.g. 1ms=30, 20ms=600(20x30), 50ms=1500(50x30))
addParameter(prs,'int',[0 2400],@isnumeric) % Plotting window [-pre post] in seconds.  
addParameter(prs,'Wcx_win',[0 0],@isnumeric) % Time window [-pre post] for Wilcoxon in seconds. Two sides should be equal.
addParameter(prs,'Wcx_alpha',0.05,@isnumeric) % alpha value for Wilcoxon, default 0.05
addParameter(prs,'norm',0,@isnumeric) % Normalise data (Z-score). Options: 1 or 0.
addParameter(prs,'average',0,@isnumeric) % Plot excited/inhibited group averages. Options: 1 or 0. 
addParameter(prs,'SignificantZscore',0,@isnumeric) % Defines significance based on provided Z-score value. Default 0.
parse(prs,varargin{:})
g = prs.Results;
disp(g.mainFolder)
PSTHall=[]; 
Wilcoxon_results = []; 
% Load cell_metrics structure
if strcmp(g.loadBatch, 'none')
    for ii = 1:length(g.Recordings)
        basepaths(ii) = {[g.mainFolder '\' g.Recordings{ii} '_kilosort\kilosort25preprocess']};
        basenames(ii) = {'temp_wh'};
    end
    cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
else
    cell_metrics = importdata(g.loadBatch);
end
% Select neurons from cell_metrics based on structure variable
if strcmp(g.selVariable, 'none')
    cellIdx_selVariable = 1:numel(cell_metrics.cellID); % all neurons
else
    cellIdx_selVariable = find(strcmp(cell_metrics.(g.selVariable), g.selValue)); 
end
% Select neuron type
if strcmp(g.neuronType, 'none')
    cellIdx_type = 1:numel(cell_metrics.cellID); % all neurons
else
    cellIdx_type = find(strcmp(cell_metrics.putativeCellType, g.neuronType)); 
end
% Find neurons that match selected variable and celltype
cellIdx_selected = intersect(cellIdx_selVariable, cellIdx_type);
% Load PSTH matrix
prevBatch = 0;
loadBar = waitbar(0,'Loading neurons...');
for hh = cellIdx_selected
    AP = cell_metrics.spikes.times{hh};
    if cell_metrics.batchIDs(hh) ~= prevBatch
        load([cell_metrics.general.basepaths{cell_metrics.batchIDs(hh)} '\TTLsKS.mat']);
        TTL = 1;
        prevBatch = cell_metrics.batchIDs(hh);
    end
    waitbar(hh/numel(cellIdx_selected),loadBar);
    % PSTH matrix
    bin_time = g.psth_bin/g.fs;     
    pre_time = abs(g.int(1));      
    post_time = g.int(2);      
    for jj = 1:numel(TTL) %Each TTL is a a column. 
         preAP{:,jj} = AP(AP>=(TTL(jj)-pre_time) & AP<TTL(jj)); % spikes before each TTL separately
         postAP{:,jj} = AP(AP>TTL(jj) & AP<(TTL(jj)+post_time)); % spikes after each TTL separately
     end
     for ll = 1:numel(TTL) %Each TTL is a a column. 
         preAP_norm{ll} = preAP{ll}-TTL(ll); % spikes relative to their own TTL
         postAP_norm{ll} = postAP{ll}-TTL(ll);
     end
     for tt = 1:numel(TTL) %Each TTL is a a column. 
         for nn = 1:(pre_time/bin_time) % number of timestamps in each bin.
             preAP_bin(nn,tt) = sum(preAP_norm{tt}>=(-pre_time+(nn-1)*bin_time) & preAP_norm{tt}<(-pre_time+nn*bin_time));
         end
         for oo = 1:(post_time/bin_time)
             postAP_bin(oo,tt) = sum(postAP_norm{tt}>=((oo-1)*bin_time) & postAP_norm{tt}<(oo*bin_time));
         end
     end  
     psth_spx_dani = vertcat(sum(postAP_bin,2));
     clear preAP_bin postAP_bin preAP postAP preAP_norm postAP_norm 

     PSTHall(find(cellIdx_selected==hh),:) = psth_spx_dani;
end

s=1