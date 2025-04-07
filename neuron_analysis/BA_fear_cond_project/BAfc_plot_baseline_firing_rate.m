function BAfc_plot_baseline_firing_rate(varargin)

prs = inputParser;
addParameter(prs,'cell_metrics',[],@isstruct) % input cell_metrics
addParameter(prs,'plot',false,@islogical)
parse(prs,varargin{:})
g = prs.Results;

cell_metrics = g.cell_metrics;

idx_IN = find(contains(cell_metrics.putativeCellType, 'IN'));
idx_PN = find(contains(cell_metrics.putativeCellType, 'PN'));

idx_LA = find(contains(cell_metrics.brainRegion, 'LA'));
idx_BA = find(contains(cell_metrics.brainRegion, 'BA'));

matrix_data = {...
    cell_metrics.firingRate(intersect(idx_IN,idx_LA))',...
    cell_metrics.firingRate(intersect(idx_PN,idx_LA))',...
    cell_metrics.firingRate(intersect(idx_IN,idx_BA))',...
    cell_metrics.firingRate(intersect(idx_PN,idx_BA))',...
    };

group_labels = [repmat(1, size(matrix_data{1},1), 1); 
                repmat(2, size(matrix_data{2},1), 1); 
                repmat(3, size(matrix_data{3},1), 1); 
                repmat(4, size(matrix_data{4},1), 1)];

all_data = [matrix_data{1}; matrix_data{2}; matrix_data{3}; matrix_data{4}];

if g.plot
    % Plot the boxplot
    figure;
    boxplot(all_data, group_labels, 'Labels', {'LA IN', 'LA PN', 'BA IN', 'BA PN'});
    
    title('Firing rate across regions');
    ylabel('Firing rate');
end