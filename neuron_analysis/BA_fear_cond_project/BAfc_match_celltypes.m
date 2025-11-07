function [cell_metrics] = BAfc_match_celltypes(varargin)

% Hardcoded script for BAfc project's waweform plot.

prs = inputParser;
addParameter(prs,'cell_metrics',[],@isstruct)
addParameter(prs, 'filePath', 'C:\Users\dmagyar\Desktop\BA_fear_cond\sorting_results.mat');
addParameter(prs, 'filePath2', 'C:\Users\dmagyar\Desktop\BA_fear_cond\BAfc_ccg_doublets.mat');
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
           cell_metrics.cellID(i) == cell_ids(j) && ...
           cell_metrics.cluID(i) == clu_ids(j)
           if trough_to_peak_corrected_ms(j) < 0.4 && cell_metrics.firingRate(i) > 10
               cell_metrics.putativeCellType{1,i} = 'IN';
           else
               cell_metrics.putativeCellType{1,i} = 'PN';
           end
           cell_metrics.spikes.ttp(i,1) = trough_to_peak_corrected_ms(j);
           cell_metrics.spikes.half_width(i,1) = half_width_corrected_ms(j);
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

load(prs.Results.filePath2);
num_matched = 0;
for i = 1:length(cell_metrics.cellID)
    for j = 1:size(potential_duplicates_phy,1)
        if strcmp(cell_metrics.animal{i}, potential_duplicates_phy{j,1}) && ...
           cell_metrics.cluID(i) == potential_duplicates_phy{j,3}
                cell_metrics.putativeCellType{1,i} = 'Doublet';
                cell_metrics.brainRegion{1,i} = 'SKIP';
                num_matched = num_matched + 1;
        end
    end
end

fprintf('removed dobulets %d / %d neurons (%.1f%%)\n', num_matched, length(cell_metrics.cellID), ...
    100 * num_matched / length(cell_metrics.cellID));
