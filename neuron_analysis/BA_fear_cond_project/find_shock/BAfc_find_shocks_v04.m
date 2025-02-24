function [psth_spx_dani, idx] =  BAfc_find_shocks_v04(experiment, time_window)

% v4: using parfor loop, faster
% experiment: 'BAfc', or 'NP_BAfc'
% time_window: time window for spike detection in ms (eg. 200)

switch experiment
    case 'BAfc'
        cell_metrics = BAfc_load_neurons;
        inputTTL = 'TTL_shocks';
    case 'NP_BAfc'
        cell_metrics = NP_BAfc_load_neurons;
        inputTTL = 'TTL_shocks_nonpredicted';
end


pre_time = 14.5; % in sec, MUST BE 14.5
post_time = 5;    
baseline_time = 10;
bin_time = 0.001; % in seconds
timewin1 = pre_time/bin_time+1:pre_time/bin_time+time_window;

% Preallocate arrays for results
num_cells = numel(cell_metrics.cellID);
psth_spx_dani = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly
smoothed_spike_train_zscore = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly

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
    psth_spx_dani(ii,:) = sum([preAP_bin; postAP_bin], 2);

    % Smoothing and Z-scoring
    sigma = 20; % standard deviation of Gaussian kernel (ms)
    binWidth = 1; % bin width (ms)
    edges = -3*sigma:binWidth:3*sigma; % kernel edges
    gaussKernel = exp(-edges.^2 / (2*sigma^2)); % Gaussian kernel
    gaussKernel = gaussKernel / sum(gaussKernel); % Normalize
    smoothed_spike_train = conv(psth_spx_dani(ii,:), gaussKernel, 'same'); % Convolve
    smoothed_spike_train_zscore(ii,:) = (smoothed_spike_train - mean(smoothed_spike_train(1:baseline_time/bin_time))) / std(smoothed_spike_train(1:baseline_time/bin_time));
end

% Find lowFR neurons in 200ms time window
window_spx = psth_spx_dani(:,timewin1);
[idx_lowFR,~] = find(sum(window_spx,2)<10);
idx_lowFR = unique(idx_lowFR);

% Find responses in 200ms time window
window = smoothed_spike_train_zscore(:,timewin1);
[idx_exc, ~] = find(window >= 5);
idx_exc = unique(idx_exc);
idx_exc = setdiff(idx_exc, idx_lowFR); % remove lowFr only from excitatiory responses.

% Find inhibitory responses 200ms time window
[idx_inh,~] = find(window<=-2);
idx_inh = unique(idx_inh);


idx = sort([idx_exc; idx_inh]);
