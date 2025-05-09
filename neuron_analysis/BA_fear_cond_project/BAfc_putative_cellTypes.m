function [cell_metrics] = BAfc_putative_cellTypes(varargin)

% Hardcoded script for BAfc project's waweform plot.

prs = inputParser;
addParameter(prs,'cell_metrics',[],@isstruct)
addParameter(prs,'plot',false,@islogical) % plot results
addParameter(prs,'width_critical',0.475,@isnumeric) % border between narrow and wide waveforms
addParameter(prs,'fr_critical',6,@isnumerci) % border between fast and regular firing rate
addParameter(prs,'ab_ratio',0.1,@isnumeric)
parse(prs,varargin{:})
g = prs.Results;

cell_metrics = g.cell_metrics;

narrow_idx = find(cell_metrics.troughToPeak<=g.width_critical);
wide_idx = find(cell_metrics.troughToPeak>g.width_critical);
% idx_fs = find(cell_metrics.firingRate>=g.fr_critical);
% idx_rs = find(cell_metrics.firingRate<g.fr_critical);
idx_sym = find(cell_metrics.ab_ratio>=g.ab_ratio);
idx_asym = find(cell_metrics.ab_ratio<g.ab_ratio);
idx_LA = find(contains(cell_metrics.brainRegion, 'LA'));
idx_BA = find(contains(cell_metrics.brainRegion, 'BA'));

idx_IN = intersect(intersect(narrow_idx, idx_sym),[idx_LA idx_BA]);
idx_PN = intersect(intersect(wide_idx, idx_asym),[idx_LA idx_BA]);
idx_unknown = setdiff([idx_LA idx_BA], [idx_IN idx_PN]);

cell_metrics.putativeCellType(idx_IN) = {'IN'};
cell_metrics.putativeCellType(idx_PN) = {'PN'};
cell_metrics.putativeCellType(idx_unknown) = {'unknown'};
