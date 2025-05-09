function [cell_metrics] = BAfc_load_neurons(varargin)

% Load cell_metrics, stimulation and LFP data
% Daniel Magyar, 3/12/2025

% Default params
prs = inputParser;
addParameter(prs,'recordings',{},@iscell)
addParameter(prs,'ttl',[],@ischar)
addParameter(prs,'LFP',false,@islogical)
addParameter(prs,'notch', false,@islogical)
addParameter(prs,'twin',[],@isnumeric) %twin for LFP around stimulus, e.g. [0.1 0.5];
addParameter(prs,'clean_metrics',false,@islogical)
addParameter(prs,'mainFolder', 'C:\Users\dmagyar\Desktop\BA_fear_cond', @ischar)
parse(prs, varargin{:})
g = prs.Results;

basenames = repmat({'temp_wh'}, 1, size(g.recordings,2));
basepaths(1:size(g.recordings,2)) = cell(1,size(g.recordings,2));
for rc = 1:size(g.recordings,2)
    basepaths{rc} = [g.mainFolder '\' g.recordings{rc} '\kilosort25preprocess'];
end
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
if ~isempty(g.ttl)
    cell_metrics.general.(g.ttl) = {};
    for ii = 1:max(cell_metrics.batchIDs)
        allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
        TTL = allTTL.(g.ttl);
        cell_metrics.general.(g.ttl)(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
    end
end
if g.LFP
    if isempty(g.ttl)
        error ('Message: No stimulus selected for LFP!')
    end
    cell_metrics = DM_load_LFP(g,cell_metrics,'notch', g.notch);    
end
% Rename cell_metrics.animal
for ii = 1:size(cell_metrics.animal,2)
    cell_metrics.animal(ii) = regexp(cell_metrics.general.batch{1, ...
        cell_metrics.batchIDs(ii)}.basepath, 'MD.{3}', 'match');
end