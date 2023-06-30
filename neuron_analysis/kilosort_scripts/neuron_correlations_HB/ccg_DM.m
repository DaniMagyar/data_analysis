function ccg_DM(experiment,recording,window,varargin)

% based on ccg.m
% https://github.com/hangyabalazs/CellBase/blob/master/CellBase_R2013a/Functions/analysis_functions/spike_correlations/ccg.m

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   14-June-2013


prs = inputParser;
addRequired(prs,'experiment')
addRequired(prs,'recording')
addRequired(prs,'window',@isscalar)  % CCG window, in seconds
addParameter(prs,'issave',false,@islogical)   % control saving behavior
addParameter(prs,'resdir',[],@(s)isdir(s)|isempty(s))   % results directory
addParameter(prs,'minspikeno',100,@isnumeric)   % calculate CCG above minimal spike number
addParameter(prs,'maxspikeno',5000,@isnumeric)   % maximize included spikes
parse(prs,experiment,recording,window,varargin{:})
g = prs.Results;
mainFolder = 'C:\Users\dmagyar\Desktop\PFC_project_base\Experiments';
cd([mainFolder '\' experiment '\' recording '_kilosort\kilosort25preprocess'])

cluster_group = tdfread('cluster_group.tsv');
cluster_info = tdfread('cluster_info.tsv');
spike_times = readNPY('spike_times.npy');
spike_clusters = readNPY('spike_clusters.npy');
groupCell = cellstr(cluster_info.group);
goodIdx = find(strcmp(groupCell, 'good')==1);
goodCluIDs = cluster_info.id(goodIdx);

PairOfCells = nchoosek(goodCluIDs,2);

% Determine time window
sr = 30000;      % sampling rate
wn = g.window * sr;    % CCG window in data points

% Cell pair loop for CCG
wb = waitbar(0,'Please wait...','Name','Running CCG...');  % progress indicator
global WB
WB(end+1) = wb;
limit_spikes = [g.minspikeno g.maxspikeno];

numPairs = size(PairOfCells,1);

for iP = 1:numPairs   % loop through pairs of cells
    cell1 = PairOfCells(iP,1);
    cell2 = PairOfCells(iP,2);
    idx1 = find(spike_clusters==cell1);
    idx2 = find(spike_clusters==cell2);
    ncc1 = double(spike_times(idx1));
    ncc2 = double(spike_times(idx2));
    ncc1 = ncc1/30000000;
    ncc2 = ncc2/30000000;

    % Implement upper spike number limits
    if length(ncc1) > limit_spikes(2)     % crop if too long to avoid out of memory
        ncc1 = ncc1(1:limit_spikes(2));
    end
    if length(ncc2) > limit_spikes(2)
        ncc2 = ncc2(1:limit_spikes(2));
    end

    if length(ncc1) > limit_spikes(1) && length(ncc2) > limit_spikes(1)     % implement minimum number of spikes
        [H1, ccr, lwr, upr, rccg] = somccg_conf_filter(ncc1,ncc2,wn,1100);    % 1->2
        xl = xlim;   % put the cell IDs in the plot
        yl = ylim;
        %text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.9,regexprep(cell1,'_',' '))
        %text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.85,regexprep(cell2,'_',' '))
        if g.issave   % save figure
            ncl1 = regexprep(cell1,'\.','_');
            ncl2 = regexprep(cell2,'\.','_');
            fnm = ['CCG_' ncl1 '_' ncl2 '.fig'];
            saveas(H1,fullfile(g.resdir,fnm))   % save CCG plot
            close(H1)
        end
        CCR(iP,:) = ccr;   % cross-correlogram
        LCCR(iP,:) = lwr;  % lower significance limit
        UCCR(iP,:) = upr;  % upper significance limit
        MeanH0(iP) = mean(mean(rccg,2),1);   % surrogate mean
        SDH0(iP) = mean(std(rccg,[],2),1);   % surrogate SD
        disp(['Pair #' num2str(iP) ' / ' num2str(numPairs) ' done......'])
    end
    
    % Save
    if g.issave
        if isequal(mod(iP,50),0)   % save after every 50 pairs to prevent data loss
            save(fullfile(g.resdir,fnm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength')
            disp('Autosave done.')
        end
    end
    
    waitbar(iP/numPairs)
end
close(wb)   % eliminate progress indicator

% Save
if g.issave
    save(fullfile(g.resdir,fnmm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength')
end