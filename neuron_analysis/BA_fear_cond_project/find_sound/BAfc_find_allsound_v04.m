function [psth_spx_dani, idx_exc200, idx_inh200] =  BAfc_find_allsound_v04(inputTTL)

% v4: using parfor loop, faster
%inputTTL = 'TTL_tone_recall_all';
cell_metrics = BAfc_load_neurons;

pre_time = 14.5; % in sec
post_time = 2;    
bin_time = 0.001; % in seconds

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
    smoothed_spike_train_zscore(ii,:) = (smoothed_spike_train - mean(smoothed_spike_train(4501:14500))) / std(smoothed_spike_train(4501:14500));
end

timewin1 = pre_time/bin_time+1:pre_time/bin_time+200;
% Find lowFR neurons in 200ms time window
window200_spx = psth_spx_dani(:,timewin1);
[idx_lowFR,~] = find(window200_spx<10);
idx_lowFR = unique(idx_lowFR);

% Find responses in 200ms time window
window200 = smoothed_spike_train_zscore(:,timewin1);
[idx_exc200, ~] = find(window200 >= 5);
idx_exc200 = unique(idx_exc200);

% Find inhibitory responses 200ms time window
[idx_inh200,~] = find(window200<=-3);
idx_inh200 = unique(idx_inh200);


