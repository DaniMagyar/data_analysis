function [cell_metrics] = BAfc_putative_cellTypes(varargin)

% Hardcoded script for BAfc project's waweform plot.

prs = inputParser;
addParameter(prs,'cell_metrics',[],@isstruct)
addParameter(prs,'plot',false,@islogical) % plot results
addParameter(prs,'width_critical',0.4,@isnumeric) % border between narrow and wide waveforms
addParameter(prs,'fr_critical',10,@isnumerci) % border between fast and regular firing rate
% Bienvenu cikkben a 10Hz jo tampont a PV IN-ekre
addParameter(prs,'ttp_d_crit',0.16,@isnumerci) % border between narrow and wide
addParameter(prs,'ab_ratio',0.1,@isnumeric)
parse(prs,varargin{:})
g = prs.Results;

cell_metrics = g.cell_metrics;

% Based on TTP
narrow_idx = find(cell_metrics.troughToPeak<=g.width_critical);
wide_idx = find(cell_metrics.troughToPeak>g.width_critical);

% % Based on TTP-derivative
% narrow_idx = find(cell_metrics.troughtoPeakDerivative>g.ttp_d_crit);
% % wide_idx = find(cell_metrics.troughtoPeakDerivative<=g.ttp_d_crit);
% 
% % Based on FR
% idx_fs = find(cell_metrics.firingRate>=g.fr_critical);
% idx_rs = find(cell_metrics.firingRate<g.fr_critical);
% idx_IN = intersect(narrow_idx, idx_fs);
% idx_PN = intersect(wide_idx, idx_rs);

% Based on assymetry
idx_sym = find(cell_metrics.ab_ratio>=g.ab_ratio);
idx_asym = find(cell_metrics.ab_ratio<g.ab_ratio);
idx_IN = intersect(narrow_idx, idx_sym);
idx_PN = intersect(wide_idx, idx_asym);

cell_metrics.putativeCellType(1:end) = {'unknown'};
cell_metrics.putativeCellType(idx_IN) = {'IN'};
cell_metrics.putativeCellType(idx_PN) = {'PN'};


% 
% 
% % cell_metrics.putativeCellType(1:end) = {'unknown'};
% % cell_metrics.putativeCellType(narrow_idx) = {'IN'};
% % cell_metrics.putativeCellType(wide_idx) = {'PN'};
% 
% % Manually fixed:
% manual{1,1} = 'MD307'; manual{1,2} = 'PN';
% 
% for ii = 1:size(manual,1)
%     if any(strcmp(cell_metrics.animal, manual{ii,1}))
%         idx = find(strcmp(cell_metrics.animal, manual{ii,1}) & cell_metrics.cellID == 13);
% %         cell_metrics.putativeCellType{idx} = manual{ii,2};
% %     end
% % end
% 
% filePath = 'C:\Users\dmagyar\Documents\data_analysis\BAfc_claude_code\uni_classifier\sorting_results_002.mat';
% data = load(filePath);
% timebin = mean(diff(data.timeaxis_super));
% idx_IN = data.trough_to_peak*timebin < 0.4 & data.half_width*timebin<0.2;
% idx_PN = data.trough_to_peak*timebin >= 0.4 & data.half_width*timebin>=0.2;
% idx_IN = data.trough_to_peak_corrected_ms < 0.4;
% idx_PN = data.trough_to_peak_corrected_ms >= 0.4;
% 
% 
% % cell_metrics.putativeCellType(1:end) = {'unknown'};
% cell_metrics.putativeCellType(1:end) = {'PN'};
% cell_metrics.putativeCellType(idx_IN) = {'IN'};
