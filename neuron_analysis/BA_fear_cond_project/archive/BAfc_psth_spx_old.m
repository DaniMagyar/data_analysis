function [psth_spx] =  BAfc_psth_spx(inputTTL, TTLinclude, pre_time, post_time, bin_time)

% v4: using parfor loop, faster
cell_metrics = BAfc_load_neurons;

% Preallocate arrays for results
num_cells = numel(cell_metrics.cellID);
psth_spx = zeros(num_cells, round((pre_time + post_time) / bin_time)); % adjust size accordingly
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
        if pre_time == 0
            preAP_bin = []
        else
            preAP_bin(:, tt) = histcounts(preAP_norm{tt}, linspace(-pre_time, 0, pre_time/bin_time + 1));
        end
        postAP_bin(:, tt) = histcounts(postAP_norm{tt}, linspace(0, post_time, post_time/bin_time + 1));
    end
    % Sum bins across TTLs
    psth_spx(ii,:) = sum([preAP_bin; postAP_bin], 2);
end