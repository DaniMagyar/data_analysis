function [potential_duplicates_phy, idx_remove] = BAfc_xcorr_neurons(varargin)

% Finds double detected neurons. The brainRegion of these neurons should be
% labeled as 'SKIP_double' or something.
% Daniel Magyar 4/2/2025.

prs = inputParser;
addParameter(prs,'cell_metrics',[],@isstruct) 
addParameter(prs,'bin_size',0.005,@isnumeric) % bin size. Default 5ms
addParameter(prs,'th',0.8,@isnumeric) % correlation threshold. Default 0.8
addParameter(prs,'max_time',1500,@isnumerci) % duration of checking. Default 1500 sec. (30 min)
parse(prs,varargin{:})
g = prs.Results;

cell_metrics = g.cell_metrics;

num_batches = max(cell_metrics.batchIDs);
potential_duplicates = cell(num_batches, 1); % Preallocate storage for batch results
loadBar = waitbar(0,'Looking for double detected neurons...');
for batch = 1:num_batches
    waitbar(batch/num_batches, loadBar);
    neurons = cell_metrics.spikes.times(cell_metrics.batchIDs == batch);
    num_neurons = length(neurons);
    bins = 0:g.bin_size:g.max_time; % 5 ms bin size
    binned_spikes = zeros(num_neurons, length(bins)-1, 'uint8');
    % Convert spike times to binned spike trains
    for ii = 1:num_neurons
        binned_spikes(ii, :) = histcounts(neurons{ii}, bins);
    end
    % Temporary storage for parallel results
    batch_results = cell(num_neurons, 1);
    % Compute cross-correlations within this batch
    parfor ii = 1:num_neurons
        local_results = {};
        for jj = ii+1:num_neurons
            [xc, ~] = xcorr(binned_spikes(ii, :), binned_spikes(jj, :), 'coeff');
            max_corr = max(xc);

            if max_corr > g.th
                local_results = [local_results; {batch, ii, jj, max_corr}]; 
            end
        end
        batch_results{ii} = local_results;
    end   
    % Store results for this batch
    potential_duplicates{batch} = vertcat(batch_results{:});
end
% Merge all results into a single array
potential_duplicates = vertcat(potential_duplicates{:});
idx_cellexp = zeros(size(potential_duplicates,1),2);
for ii = 1:size(potential_duplicates,1)
    idx_cellexp(ii,1) = intersect(find(cell_metrics.batchIDs == potential_duplicates{ii,1}), ...
        find(cell_metrics.cellID == potential_duplicates{ii,2}));
    idx_cellexp(ii,2) = intersect(find(cell_metrics.batchIDs == potential_duplicates{ii,1}), ...
        find(cell_metrics.cellID == potential_duplicates{ii,3}));
end
close(loadBar)
potential_duplicates_phy = cell(size(potential_duplicates));
for ii = 1:size(potential_duplicates,1)
    potential_duplicates_phy{ii,1} = cell_metrics.animal{find(cell_metrics.batchIDs == potential_duplicates{ii,1},1)};
    potential_duplicates_phy{ii,2} = cell_metrics.cluID(idx_cellexp(ii,1));
    potential_duplicates_phy{ii,3} = cell_metrics.cluID(idx_cellexp(ii,2));
    potential_duplicates_phy{ii,4} = potential_duplicates{ii,4};
end
% removing unique elements from the second column. Usually the first index
% is better quality because it was found first in Phy. Higher indexes were
% probably detected with lower amplitude on an other channel.
idx_remove = unique(idx_cellexp(:,2));