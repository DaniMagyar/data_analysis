function [cell_metrics] = BAfc_load_neurons(varargin)

% Load cell_metrics, stimulation and LFP data
% Daniel Magyar, 3/12/2025

% Default params
prs = inputParser;
addParameter(prs,'recordings',{},@iscell)
addParameter(prs,'ttl',{},@iscell)
addParameter(prs,'LFP',false,@islogical)
addParameter(prs,'notch', false,@islogical)
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
    for tt = 1:size(g.ttl,2)
        cell_metrics.general.(g.ttl{tt}) = {};
        for ii = 1:max(cell_metrics.batchIDs)
            allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
            TTL = allTTL.(g.ttl{tt});
            cell_metrics.general.(g.ttl{tt})(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
        end
    end
end

% % Rename cell_metrics.animal
% for ii = 1:size(cell_metrics.animal,2)
%     recID = regexp(cell_metrics.general.batch{1, ...
%         cell_metrics.batchIDs(ii)}.basepath, 'MD.{27}', 'match');
%     cell_metrics.animal(ii) = {recID{1}([1:5, 26:29])};
%     % cell_metrics.animal(ii) = regexp(cell_metrics.general.batch{1, ...
%     %     cell_metrics.batchIDs(ii)}.basepath, 'MD.{27}', 'match');
% end
animalIDs = {};
for ii = 1:size(g.recordings,2)
    animalIDs = [animalIDs repmat({g.recordings{ii}(1:9)}, 1, cell_metrics.general.batch_benchmark.file_cell_count(ii))];
end
cell_metrics.animal = animalIDs;
cell_metrics.general.mainFolder = g.mainFolder;