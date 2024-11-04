function [responses_shock, responses_shock_class] =  BAfc_find_shocks_v01(TTLinclude)

% v02: Uses Gaussian Kernel filtering, not matching the cell_metrics figure
% method, so use v3 instead of this

TTLinclude = 1:10;

[psth_spx] = BAfc_calc_Zscore_v02('Stims', {'TTL_shocks'},...
    'TTLselect', TTLinclude, 'psth_bin', 30, 'int', [-14.5 2], 'norm', 0); % each ttl is a column
% [response_data2] = BAfc_calc_Zscore_v02('Stims', {'TTL_shocks'},...
%    'TTLselect', TTLinclude, 'psth_bin', 30, 'int', [-0 2], 'norm', 0); % each ttl is a column.
cell_metrics = BAfc_load_neurons;
cell_metrics.labels(1:end) = {'neutral'};


sigma = 20; % standard deviation of the Gaussian kernel (in ms)
binWidth = 1; % bin width for the spike train (in ms)
edges = -3*sigma:binWidth:3*sigma; % create kernel edges
gaussKernel = exp(-edges.^2 / (2*sigma^2)); % Gaussian formula
gaussKernel = gaussKernel / sum(gaussKernel); % normalize kernel


for ii = 1: size(psth_spx,1)
    smoothed_spx(ii,:) = conv(psth_spx(ii,:), gaussKernel, 'same'); % convolve
    smoothed_spx_zscore(ii,:) = (smoothed_spx(ii,:) - mean(smoothed_spx(ii,1:10000)))/std(smoothed_spx(ii,1:10000));
end

[idx_exc,~] = find(smoothed_spx_zscore(:,14501:16500)>=5);
idx_exc = sort(idx_exc);
repeats_exc = diff(idx_exc) == 0; 
idx_exc(repeats_exc) = [];

[idx_exc100,~] = find(smoothed_spx_zscore(:,14501:14600)>=5);
idx_exc100 = sort(idx_exc100);
repeats_exc = diff(idx_exc100) == 0; 
idx_exc100(repeats_exc) = [];


cell_metrics.labels(idx_exc) = {'exc'};
cell_metrics.labels(idx_exc100) = {'exc100'};

cell_metrics = CellExplorer('metrics',cell_metrics);