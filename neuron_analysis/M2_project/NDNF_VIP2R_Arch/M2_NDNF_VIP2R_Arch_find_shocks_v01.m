function [psth_spx_1, psth_spx_2, idx] =  M2_NDNF_VIP2R_Arch_find_shocks_v01(pre_time, post_time, bin_time)

cell_metrics = M2_NDNF_VIP2R_Arch_load_neurons;

cell_metrics.labels(1:end) = {'neutral'};
% pre_time = 5; % in sec
% post_time = 5;  
% bin_time = 0.001; % in seconds
illum_start = pre_time - 1;
timewin1 = pre_time/bin_time+1:pre_time/bin_time+500;
zscore_exc = 5;
zscore_inh = -3;
% first TTL
% Preallocate arrays for results
num_cells = numel(cell_metrics.cellID);
psth_spx_1 = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly
smoothed_spike_train_zscore_1 = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly
inputTTL = 'TTL_shocks_only';
parfor ii = 1:num_cells
    disp(ii)
    TTL = cell_metrics.general.(inputTTL){ii};
    AP = cell_metrics.spikes.times{ii};

    % Preallocate cell arrays for spikes
    preAP = cell(numel(TTL), 1);
    postAP = cell(numel(TTL), 1);
    preAP_norm = cell(numel(TTL), 1);
    postAP_norm = cell(numel(TTL), 1);

    % Binning variables
    preAP_bin = zeros(pre_time/bin_time, numel(TTL));
    postAP_bin = zeros(post_time/bin_time, numel(TTL));

    % Loop over TTL events
    for jj = 1:numel(TTL)
        % Spikes before and after each TTL
        preAP{jj} = AP(AP >= (TTL(jj) - pre_time) & AP < TTL(jj));
        postAP{jj} = AP(AP > TTL(jj) & AP < (TTL(jj) + post_time));

        % Normalize spike times to TTL
        preAP_norm{jj} = preAP{jj} - TTL(jj);
        postAP_norm{jj} = postAP{jj} - TTL(jj);
    end

    % Bin spike times
    for tt = 1:numel(TTL)
        preAP_bin(:, tt) = histcounts(preAP_norm{tt}, linspace(-pre_time, 0, pre_time/bin_time + 1));
        postAP_bin(:, tt) = histcounts(postAP_norm{tt}, linspace(0, post_time, post_time/bin_time + 1));
    end

    % Sum bins across TTLs
    psth_spx_1(ii,:) = sum([preAP_bin; postAP_bin], 2);

    % Smoothing and Z-scoring
    sigma = 20; % standard deviation of Gaussian kernel (ms)
    binWidth = 1; % bin width (ms)
    edges = -3*sigma:binWidth:3*sigma; % kernel edges
    gaussKernel = exp(-edges.^2 / (2*sigma^2)); % Gaussian kernel
    gaussKernel = gaussKernel / sum(gaussKernel); % Normalize
    smoothed_spike_train_1 = conv(psth_spx_1(ii,:), gaussKernel, 'same'); % Convolve
    smoothed_spike_train_zscore_1(ii,:) = (smoothed_spike_train_1 - mean(smoothed_spike_train_1(1:illum_start/bin_time))) / std(smoothed_spike_train_1(1:illum_start/bin_time));
end
% Find responses in 500ms time window
window200_1 = smoothed_spike_train_zscore_1(:,timewin1);
[idx_exc1,~] = find(window200_1>=zscore_exc);
idx_exc1 = unique(idx_exc1);
cell_metrics.labels(idx_exc1) = {'exc'};

% Find inhibitory responses 500ms time window
[idx_inh1,~] = find(window200_1<=zscore_inh);
idx_inh1 = unique(idx_inh1);
cell_metrics.labels(idx_inh1) = {'inh'};




% second TTL
% Preallocate arrays for results
num_cells = numel(cell_metrics.cellID);
psth_spx_2 = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly
smoothed_spike_train_zscore_2 = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly
inputTTL = 'TTL_shocks_light';
parfor ii = 1:num_cells
    disp(ii)
    TTL = cell_metrics.general.(inputTTL){ii};
    AP = cell_metrics.spikes.times{ii};

    % Preallocate cell arrays for spikes
    preAP = cell(numel(TTL), 1);
    postAP = cell(numel(TTL), 1);
    preAP_norm = cell(numel(TTL), 1);
    postAP_norm = cell(numel(TTL), 1);

    % Binning variables
    preAP_bin = zeros(pre_time/bin_time, numel(TTL));
    postAP_bin = zeros(post_time/bin_time, numel(TTL));

    % Loop over TTL events
    for jj = 1:numel(TTL)
        % Spikes before and after each TTL
        preAP{jj} = AP(AP >= (TTL(jj) - pre_time) & AP < TTL(jj));
        postAP{jj} = AP(AP > TTL(jj) & AP < (TTL(jj) + post_time));

        % Normalize spike times to TTL
        preAP_norm{jj} = preAP{jj} - TTL(jj);
        postAP_norm{jj} = postAP{jj} - TTL(jj);
    end

    % Bin spike times
    for tt = 1:numel(TTL)
        preAP_bin(:, tt) = histcounts(preAP_norm{tt}, linspace(-pre_time, 0, pre_time/bin_time + 1));
        postAP_bin(:, tt) = histcounts(postAP_norm{tt}, linspace(0, post_time, post_time/bin_time + 1));
    end

    % Sum bins across TTLs
    psth_spx_2(ii,:) = sum([preAP_bin; postAP_bin], 2);

    % Smoothing and Z-scoring
    sigma = 20; % standard deviation of Gaussian kernel (ms)
    binWidth = 1; % bin width (ms)
    edges = -3*sigma:binWidth:3*sigma; % kernel edges
    gaussKernel = exp(-edges.^2 / (2*sigma^2)); % Gaussian kernel
    gaussKernel = gaussKernel / sum(gaussKernel); % Normalize
    smoothed_spike_train_2 = conv(psth_spx_2(ii,:), gaussKernel, 'same'); % Convolve
    smoothed_spike_train_zscore_2(ii,:) = (smoothed_spike_train_2 - mean(smoothed_spike_train_2(1:illum_start/bin_time))) / std(smoothed_spike_train_2(1:illum_start/bin_time));
end
% Find responses in 500ms time window
window200_2 = smoothed_spike_train_zscore_2(:,timewin1);
[idx_exc2,~] = find(window200_2>=zscore_exc);
idx_exc2 = unique(idx_exc2);
cell_metrics.labels(idx_exc2) = {'exc'};

% Find inhibitory responses 500ms time window
[idx_inh2,~] = find(window200_2<=zscore_inh);
idx_inh2 = unique(idx_inh2);
cell_metrics.labels(idx_inh2) = {'inh'};


idx = unique(sort([idx_exc1; idx_exc2; idx_inh1; idx_inh2]));