function [PSTHall, psth_spx_bin] = BAfc_calc_Zscore_v01(varargin)

prs =  inputParser;
addParameter(prs,'mainFolder','C:\Users\dmagyar\Desktop\BA_fear_cond',@ischar) % Path of main folder. (e.g. 'C:\Users\dmagyar\Desktop\M2_shock_response')
addParameter(prs,'Stims',@iscell) % TTLs to use from TTLsKS.mat. (e.g. 'shocks')
addParameter(prs,'TTLselect',[0],@isnumeric) % Selected TTLs (e.g. [1:3])
addParameter(prs,'selVariable','none',@ischar) % Variable used for neuron selection from cell_metrics. (e.g. 'labels')
addParameter(prs,'selValue','none',@ischar) % Value used for neuron selection from 'selectVariable'. (e.g. '2')
addParameter(prs,'neuronType','none',@ischar) % Options: 'Narrow Interneuron', 'Wide Interneuron', 'Pyramidal Cell'. Default 'none' loads all types. 
addParameter(prs,'fs',30000,@isnumeric) % Sampling rate, default 30000.
addParameter(prs,'psth_bin',15000,@isnumeric) % Binning window calculated from fs. (e.g. 1ms=30, 20ms=600(20x30), 50ms=1500(50x30))
addParameter(prs,'int',[-20 4.5],@isnumeric) % Plotting window [-pre post] in seconds.  
addParameter(prs,'Wcx_win',[-5 5],@isnumeric) % Time window [-pre post] for Wilcoxon in seconds. Two sides should be equal.
addParameter(prs,'Wcx_alpha',0.05,@isnumeric) % alpha value for Wilcoxon, default 0.05
addParameter(prs,'norm',1,@isnumeric) % Normalise data (Z-score). Options: 1 or 0.
addParameter(prs,'offset',1,@isnumeric) % Offset to pre stimulus firing rate ( int(1)).
addParameter(prs,'Sorted',0,@isnumeric) % If 1, Load 'SortIDX' from 'BAparams.mat'. Default 0.
addParameter(prs,'SignificantZscore',0,@isnumeric) % Defines significance based on provided Z-score value. Default 0.

parse(prs,varargin{:})
g = prs.Results;

cell_metrics = BAfc_load_neurons;

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
% Remove 'junk' 
cellIdx_keep =  1:numel(cell_metrics.cellID); % all neurons
cellIdx_keep(find(strcmp(cell_metrics.putativeCellType, 'junk'))) = []; % neurons labeled as 'junk'
% Find neurons that match selected variable and celltype and not junk
cellIdx_selected = intersect(intersect(cellIdx_selVariable, cellIdx_type), cellIdx_keep);

% Load PSTH matrix
prevBatch = 0;
loadBar = waitbar(0,'Loading neurons...');
Significance(1:numel(cellIdx_selected)) = 0;

for hh = cellIdx_selected
    preAP_bin = [];
    postAP_bin = [];
    TTL = [];
    AP = cell_metrics.spikes.times{hh};
    for ii = 1:length(g.Stims)
        ttl = cell_metrics.general.(g.Stims{ii}){hh};
        TTL = vertcat(TTL,ttl);
    end
    TTL = sort(TTL);
    if g.TTLselect ~= 0 
        TTL = TTL(g.TTLselect);
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
    psth_spx_dani = vertcat(sum(preAP_bin,2), sum(postAP_bin,2));
    psth_spx_bin(hh,1) = {preAP_bin};
    psth_spx_bin(hh,2) = {postAP_bin};
    clear preAP_bin postAP_bin preAP postAP preAP_norm postAP_norm  % reset matrices, because different TTL num    
    % Joshua Johansen method

    psth_spx_dani = psth_spx_dani/numel(TTL);
    n_baseline_bin = (abs(g.int(1))*g.fs)/g.psth_bin;
    % μ and σ are the mean and s.d.,respectively, of all Si values in a set of baseline bins
    if g.norm == 1
        psth_spx_zscore = (psth_spx_dani - mean(psth_spx_dani(1:n_baseline_bin)))/std(psth_spx_dani(1:n_baseline_bin));
        PSTHall(find(cellIdx_selected==hh),:) = psth_spx_zscore';
    elseif g.norm == 0 
        PSTHall(find(cellIdx_selected==hh),:) = psth_spx_dani';
    end
    PSTHallIDs(find(cellIdx_selected==hh),1) = {['Batch ' num2str(cell_metrics.batchIDs(hh)) ' Cell ' num2str(cell_metrics.cellID(hh))]};
    PSTHallIDs(find(cellIdx_selected==hh),2) = cell_metrics.brainRegion(hh);
end
close(loadBar)