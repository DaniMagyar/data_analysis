function [psth_spx_dani, idx] =  BAfc_find_sound_v04(experiment, inputTTL, TTLinclude)

% v4: using parfor loop, faster
%inputTTL = 'TTL_tone_habit_first';
switch experiment
    case 'BAfc'
        cell_metrics = BAfc_load_neurons;
    case 'NP_BAfc'
        cell_metrics = NP_BAfc_load_neurons;
end

pre_time = 14.5; % in sec, MUST BE 14.5, HARDCODED
post_time = 14.5;   % MUST BE 14.5
baseline_time = 10;
bin_time = 0.001; % in seconds

% Preallocate arrays for results
num_cells = numel(cell_metrics.cellID);
psth_spx_dani = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly
smoothed_zscore_short = zeros(num_cells, (pre_time + post_time) / bin_time); % adjust size accordingly

parfor ii = 1:num_cells
    TTL = cell_metrics.general.(inputTTL){ii}(TTLinclude);
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

    %% Smoothing and Z-scoring short time window
    sigma1 = 20; % standard deviation of Gaussian kernel (ms)
    binWidth1 = 1; % bin width (ms)
    edges1 = -3*sigma1:binWidth1:3*sigma1; % kernel edges
    gaussKernel1 = exp(-edges1.^2 / (2*sigma1^2)); % Gaussian kernel
    gaussKernel1 = gaussKernel1 / sum(gaussKernel1); % Normalize
    smoothed_spike_train1 = conv(psth_spx_dani(ii,:), gaussKernel1, 'same'); % Convolve
    baseline1 = smoothed_spike_train1(pre_time/bin_time-baseline_time/bin_time+1:pre_time/bin_time);
    smoothed_zscore_short(ii,:) = (smoothed_spike_train1 - mean(baseline1))/std(baseline1);

    %% Smoothing and Z-scoring large time window
    sigma2 = 1000; % standard deviation of Gaussian kernel (ms)
    binWidth2 = 1; % bin width (ms)
    edges2 = -3*sigma2:binWidth2:3*sigma2; % kernel edges
    gaussKernel2 = exp(-edges2.^2 / (2*sigma2^2)); % Gaussian kernel
    gaussKernel2 = gaussKernel2 / sum(gaussKernel2); % Normalize
    smoothed_spike_train2 = conv(psth_spx_dani(ii,:), gaussKernel2, 'same'); % Convolve
    baseline2 = smoothed_spike_train2(pre_time/bin_time-baseline_time/bin_time+1:pre_time/bin_time);
    smoothed_zscore_long(ii,:) = (smoothed_spike_train2 - mean(baseline2))/std(baseline2);

end

%% Short timescale (200ms) calculations
% Find lowFR neurons in 200ms time window
% window200_spx = horzcat(...
%     psth_spx_dani(:,4501:4700), ...
%     psth_spx_dani(:,5501:5700), ...
%     psth_spx_dani(:,6501:6700),...
%     psth_spx_dani(:,7501:7700),...
%     psth_spx_dani(:,8501:8700));
% 
% window200_spx = sum(window200_spx,2);
% [idx_lowFR_short,~] = find(window200_spx<10);
% idx_lowFR_short = unique(idx_lowFR_short);

window200 = horzcat(...
    smoothed_zscore_short(:,pre_time/bin_time+1:pre_time/bin_time+200), ...
    smoothed_zscore_short(:,pre_time/bin_time+1001:pre_time/bin_time+1200), ...
    smoothed_zscore_short(:,pre_time/bin_time+2001:pre_time/bin_time+2200),...
    smoothed_zscore_short(:,pre_time/bin_time+3001:pre_time/bin_time+3200),...
    smoothed_zscore_short(:,pre_time/bin_time+4001:pre_time/bin_time+4200));

% Find responses in 200ms time window
[idx_exc200,~] = find(window200>=5);
idx_exc200 = unique(idx_exc200);

% Find inhibitory responses 200ms time window
[idx_inh200,~] = find(window200<=-3);
idx_inh200 = unique(idx_inh200);

%% Large timescale calculations
% Find lowFR neurons in 4.5s time window
window4500_spx = psth_spx_dani(:,pre_time/bin_time-baseline_time/bin_time+1:pre_time/bin_time);


minSpkN = numel(cell_metrics.general.(inputTTL){1}) * baseline_time * 0.1; % below 0.1Hz on average filtered
window4500_spx = sum(window4500_spx,2);
[idx_lowFR_long,~] = find(window4500_spx<minSpkN); 
idx_lowFR_long = unique(idx_lowFR_long);

window4500 = smoothed_zscore_long(:,pre_time/bin_time+1:pre_time/bin_time+4500);

% Find responses in 200ms time window
[idx_exc4500,~] = find(window4500>=3);
idx_exc4500 = unique(idx_exc4500);

% Find inhibitory responses 200ms time window
[idx_inh4500,~] = find(window4500<=-2);
idx_inh4500 = unique(idx_inh4500);

%% Remove LowFR neurons

idx_exc200 = setdiff(idx_exc200, idx_lowFR_long); 
idx_inh200 = setdiff(idx_inh200, idx_lowFR_long); 
idx_exc4500 = setdiff(idx_exc4500, idx_lowFR_long);
idx_inh4500 = setdiff(idx_inh4500, idx_lowFR_long);

% idx_exc = unique([idx_exc200; idx_exc4500]);
% idx_inh = unique([idx_inh200; idx_inh4500]);

idx = [idx_exc200; idx_inh200];
 


%% Optional cell_explorer
% 
% cell_metrics.Excited(1:numel(cell_metrics.cellID)) = 0;
% cell_metrics.Excited(idx_exc200) = 200;
% cell_metrics.Excited(idx_exc4500) = 450;
% cell_metrics.Excited(intersect(idx_exc200, idx_exc4500)) = 900;
% 
% cell_metrics.Inhibited(1:numel(cell_metrics.cellID)) = 0;
% cell_metrics.Inhibited(idx_inh200) = 200;
% cell_metrics.Inhibited(idx_inh4500) = 450;
% cell_metrics.Inhibited(intersect(idx_inh200, idx_inh4500)) = 900;
% 
% disp(['w/oLFR_Excited = ' num2str(numel(cell_metrics.Excited) - numel(find(cell_metrics.Excited == 0)))])
% disp(['w/oLFR_Inhibited = ' num2str(numel(cell_metrics.Inhibited) - numel(find(cell_metrics.Inhibited == 0)))])
% 
% cell_metrics = CellExplorer('metrics',cell_metrics);
% 


% a Low-FR-t visszairja 900 ra, azt megkell csinalni. 