function [response_dir] = BAfc_calc_Zscore_JJohansen_trial(varargin)

prs =  inputParser;
addParameter(prs,'mainFolder','C:\Users\dmagyar\Desktop\BA_fear_cond',@ischar) % Path of main folder. (e.g. 'C:\Users\dmagyar\Desktop\M2_shock_response')
addParameter(prs,'Stims',@iscell) % TTLs to use from TTLsKS.mat. (e.g. 'shocks')
addParameter(prs,'TTLselect',[0],@isnumeric) % Selected TTLs (e.g. [1:3])
addParameter(prs,'selVariable','none',@ischar) % Variable used for neuron selection from cell_metrics. (e.g. 'labels')
addParameter(prs,'selValue','none',@ischar) % Value used for neuron selection from 'selectVariable'. (e.g. '2')
addParameter(prs,'neuronType','none',@ischar) % Options: 'Narrow Interneuron', 'Wide Interneuron', 'Pyramidal Cell'. Default 'none' loads all types. 
addParameter(prs,'fs',30000,@isnumeric) % Sampling rate, default 30000.
addParameter(prs,'psth_bin',3000,@isnumeric) % Binning window calculated from fs. (e.g. 1ms=30, 20ms=600(20x30), 50ms=1500(50x30))
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
    clear preAP_bin postAP_bin preAP postAP preAP_norm postAP_norm  % reset matrices, because different TTL num    
    % Joshua Johansen method

    psth_spx_dani = psth_spx_dani/numel(TTL);
    n_baseline_bin = (abs(g.int(1))*g.fs)/g.psth_bin;
    % μ and σ are the mean and s.d.,respectively, of all Si values in a set of baseline bins
    psth_spx_zscore = (psth_spx_dani - mean(psth_spx_dani(1:n_baseline_bin)))/std(psth_spx_dani(1:n_baseline_bin));
    PSTHall(find(cellIdx_selected==hh),:) = psth_spx_zscore';
    PSTHallIDs(find(cellIdx_selected==hh),1) = {['Batch ' num2str(cell_metrics.batchIDs(hh)) ' Cell ' num2str(cell_metrics.cellID(hh))]};
    PSTHallIDs(find(cellIdx_selected==hh),2) = cell_metrics.brainRegion(hh);
    
    trial_length = 4.5; %trial length is seconds
    n_trial_bin = (trial_length*g.fs)/g.psth_bin;
    test_bin = psth_spx_zscore(n_baseline_bin+1:n_baseline_bin+n_trial_bin);
    
    % Significant #1: if 1 bin is above 3 zscore 
    if max(test_bin) > 3
        Significance(hh) = 3;
    end
    % Significant #2: if 2 consecutive bins are above 2 zscore
    idx1 = find((test_bin)>2);
    idx1_diff = diff(idx1);
    if any(idx1_diff == 1)
       Significance(hh) = 2;
    end
    idx2 = find((test_bin)>1);
    idx2_diff = diff(idx2);
    lookfor = [1 1 1];
    out = strfind(idx2_diff', lookfor);
    if any(out >= 1)
        Significance(hh) = 1;
    end

    if min(test_bin) <= -3
        Significance(hh) = -3;
    end
    % Significant #2: if 2 consecutive bins are above 2 zscore
    idx3 = find((test_bin)<=-2);
    idx3_diff = diff(idx3);
    if any(idx3_diff == 1)
       Significance(hh) = -2;
    end
    idx4 = find((test_bin)<=-1);
    idx4_diff = diff(idx4);
    lookfor = [1 1 1];
    out = strfind(idx4_diff', lookfor);
    if any(out >= 1)
        Significance(hh) = -1;
    end




end














close(loadBar)
time = (g.int(1)*g.fs:g.psth_bin:g.int(2)*g.fs)';
time = time(1:end-1);
timeline = time/g.fs+g.psth_bin/g.fs;
% Colormap blue:white:red -------------------------------------------------
c1 = 1/255*[0,115,185];
c2 = 1/255*[239,239,239];
c3 = 1/255*[217,84,26];
mycolormap(1:20,1)  = linspace(c1(1),c2(1),20);
mycolormap(21:64,1) = linspace(c2(1),c3(1),44);
mycolormap(1:20,2)  = linspace(c1(2),c2(2),20);
mycolormap(21:64,2) = linspace(c2(2),c3(2),44);
mycolormap(1:20,3)  = linspace(c1(3),c2(3),20);
mycolormap(21:64,3) = linspace(c2(3),c3(3),44);  
% Calculating plotting order based on Z-score change
testWindow_firstBin = abs(g.int(1))*g.fs/g.psth_bin+1;
testBins = g.Wcx_win(2)*g.fs/g.psth_bin; % number of bins in Wcx test window. Z-score ordering, Wcx direction based on this
testWindow_lastBin = testWindow_firstBin + (testBins-1);    
testWindow = PSTHall(:, testWindow_firstBin:testWindow_lastBin);
testMean = mean(testWindow,2);
switch g.Sorted
    case 0
        [~,SortIDX] = sort(testMean, 'descend');
        MyNewOrder = PSTHallIDs(SortIDX,:);
        save ([g.mainFolder '\BAparams_' g.selValue '.mat'], 'SortIDX','MyNewOrder')
    case 1
        load ([g.mainFolder '\BAparams_' g.selValue '.mat'], 'SortIDX','MyNewOrder')
end
 
y = imagesc(time/g.fs, 1:size(PSTHall,1), PSTHall(SortIDX,:)); % plotting sorted psth matrix 
clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
colormap(mycolormap); 
% cb = colorbar('FontSize',10);
% cb.Ticks = linspace(-20,20,41);
% cb.Label.String = 'Z-score change';
% cb.Label.FontSize = 15;
% cb.Position = [0.75, 0.215, 0.01, 0.255];
hold on;
plot([-0.0125 -0.0125],[1 size(PSTHall,1)],'r', 'LineWidth', 1.5);
hold off;
ylabel('# Cell');
xlabel('Time (s)');
set(gca,'FontSize',8);   





















