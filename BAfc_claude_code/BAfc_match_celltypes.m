function [cell_metrics] = BAfc_match_celltypes(varargin)

% Hardcoded script for BAfc project's waweform plot.

prs = inputParser;
addParameter(prs,'cell_metrics',[],@isstruct)
addParameter(prs, 'filePath', 'C:\Users\dmagyar\Desktop\BA_fear_cond\sorting_results.mat');
parse(prs,varargin{:})

load(prs.Results.filePath);
cell_metrics = prs.Results.cell_metrics;
%% Match neurons and copy putativeCellType

fprintf('\nMatching neurons between cell_metrics and putative_celltypes...\n');

% Initialize putativeCellType field
cell_metrics.putativeCellType = cell(1,length(cell_metrics.cellID));

% Match neurons based on: animal, brainRegion, cellID, cluID
num_matched = 0;
for i = 1:length(cell_metrics.cellID)
    matched = false;
    for j = 1:length(cell_ids)
        if strcmp(cell_metrics.animal{i}, animals{j}) && ...
           strcmp(cell_metrics.brainRegion{i}, brain_regions{j}) && ...
           cell_metrics.cellID(i) == cell_ids(j) && ...
           cell_metrics.cluID(i) == clu_ids(j)
           if trough_to_peak_corrected_ms(j) < 0.4
               cell_metrics.putativeCellType{1,i} = 'IN';
           else
               cell_metrics.putativeCellType{1,i} = 'PN';
           end
           num_matched = num_matched + 1;
           matched = true;
           break;
        end
    end

    if ~matched
        cell_metrics.putativeCellType{1,i} = 'Unmatched';
    end
end

fprintf('Matched %d / %d neurons (%.1f%%)\n', num_matched, length(cell_metrics.cellID), ...
    100 * num_matched / length(cell_metrics.cellID));