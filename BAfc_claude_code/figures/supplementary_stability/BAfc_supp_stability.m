% BAfc_supp_stability
% Visualize stability of neuronal responses across trial blocks
% Shows that responses remained stable during the experiment (did not change over time)

clear all; close all

%% Setup
recordings = {...
    'MD292_002_kilosort',...
    'MD293_001_kilosort',...
    'MD294_001_kilosort',...
    'MD295_001_kilosort',...
    'MD296_001_kilosort',...
    'MD297_001_kilosort',...
    'MD298_001_kilosort',...
    'MD299_001_kilosort',...
    'MD300_001_kilosort',...
    'MD304_001_kilosort',...
    'MD305_001_kilosort',...
    'MD307_001_kilosort',...
    'MD309_001_kilosort',...
    'MD310_001_kilosort',...
    'MD311_002_kilosort',...
    'MD312_001_kilosort',...
    'MD313_001_kilosort',...
    'MD314_001_kilosort',...
    'MD315_001_kilosort',...
    'MD316_002_kilosort',...
    'MD317_001_kilosort',...
    'MD318_001_kilosort',...
    'MD318_002_kilosort',...
    'MD319_003_kilosort'};

ttl = {'triptest_sound_only','triptest_shocks_only','triptest_both'};
stim_names = {'CS', 'US', 'CS+US'};

g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
g.cell_metrics = BAfc_match_celltypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 10;
g.fontSize2 = 10;
g.xlinewidth = 2;
g.pre_time = 5;
g.post_time = 4;
g.test_time = 1;
g.bin_time = 0.001;

% Response detection parameters (same as BAfc_figure_3)
g.excitation_threshold = 2;
g.inhibition_fr_drop = 0.50;
g.onset_threshold = g.excitation_threshold;
g.min_consec_bins = max(1, round(0.020 / g.bin_time));
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;

% Brain regions to analyze
brain_regions = {'LA', 'BA', 'Astria'};

% Region colors
region_colors = struct();
region_colors.LA = [0.7 0.2 0.2];
region_colors.BA = [0.2 0.7 0.2];
region_colors.Astria = [0.2 0.4 0.7];

%% Calculate PSTHs and get trial-by-trial spikes
fprintf('Calculating PSTHs and extracting trial data...\n');
psthZ_full = cell(1, numel(ttl));
psthHz_full = cell(1, numel(ttl));
preAP_norm_all = cell(1, numel(ttl));
postAP_norm_all = cell(1, numel(ttl));

for hmp = 1:numel(ttl)
    [psth_spx_og, preAP_norm_all{hmp}, postAP_norm_all{hmp}] = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

    % Convert to Hz
    num_trials = size(g.cell_metrics.general.(ttl{hmp}){1}, 1);
    psth_hz = psth_spx_og / (num_trials * g.bin_time);

    % Smooth Hz data
    psth_hz_smooth = smoothdata(psth_hz, 2, 'sgolay', 201);
    psthHz_full{hmp} = psth_hz_smooth;

    % Z-score
    baseline_idx = 1:(g.pre_time / g.bin_time);
    baseline_mean = mean(psth_spx_og(:, baseline_idx), 2);
    baseline_std = std(psth_spx_og(:, baseline_idx), 0, 2);
    psth_spx = (psth_spx_og - baseline_mean) ./ baseline_std;
    psth_spx = smoothdata(psth_spx, 2, 'sgolay', 201);
    psthZ_full{hmp} = psth_spx;
end

%% Process each brain region and identify responsive neurons
results_all = cell(1, 3);

for br = 1:3
    fprintf('\nProcessing %s...\n', brain_regions{br});

    % Get neuron indices
    idx_neurons = strcmp(g.cell_metrics.brainRegion, brain_regions{br});
    n_neurons = sum(idx_neurons);
    fprintf('  %d neurons\n', n_neurons);

    if n_neurons == 0
        continue;
    end

    % Extract PSTHs for CS and US (for classification)
    psth_CS = psthZ_full{1}(idx_neurons, :);
    psth_US = psthZ_full{2}(idx_neurons, :);
    psth_CS_Hz = psthHz_full{1}(idx_neurons, :);
    psth_US_Hz = psthHz_full{2}(idx_neurons, :);

    % Calculate responses based on CS and US only
    CS_peak = max(psth_CS(:, g.roi), [], 2);
    US_peak = max(psth_US(:, g.roi), [], 2);

    baseline_idx = 1:(g.pre_time / g.bin_time);
    CS_baseline_fr = mean(psth_CS_Hz(:, baseline_idx), 2);
    US_baseline_fr = mean(psth_US_Hz(:, baseline_idx), 2);
    CS_test_fr = mean(psth_CS_Hz(:, g.roi), 2);
    US_test_fr = mean(psth_US_Hz(:, g.roi), 2);
    CS_fr_drop = (CS_baseline_fr - CS_test_fr) ./ (CS_baseline_fr + eps);
    US_fr_drop = (US_baseline_fr - US_test_fr) ./ (US_baseline_fr + eps);

    % Classify neurons (same as BAfc_figure_3)
    CS_excited = CS_peak >= g.excitation_threshold;
    US_excited = US_peak >= g.excitation_threshold;
    CS_inhibited = CS_fr_drop >= g.inhibition_fr_drop;
    US_inhibited = US_fr_drop >= g.inhibition_fr_drop;

    Clusters = zeros(n_neurons, 1);
    Clusters(CS_excited & ~US_excited) = 1;  % CS-selective
    Clusters(US_excited & ~CS_excited) = 2;  % US-selective
    Clusters(CS_excited & US_excited) = 3;   % Multisensory
    Clusters(~CS_excited & ~US_excited & (CS_inhibited | US_inhibited)) = 5;  % Inhibited
    Clusters(~CS_excited & ~US_excited & ~CS_inhibited & ~US_inhibited) = 4;  % Non-responsive

    % Get responsive neurons (clusters 1-3 only)
    responsive_mask = ismember(Clusters, [1, 2, 3]);
    responsive_idx = find(idx_neurons);
    responsive_idx = responsive_idx(responsive_mask);
    n_responsive = length(responsive_idx);

    fprintf('  %d responsive neurons (clusters 1-3)\n', n_responsive);

    % Store results
    results_all{br}.responsive_idx = responsive_idx;
    results_all{br}.n_responsive = n_responsive;
    results_all{br}.Clusters = Clusters(responsive_mask);
end

%% Calculate block-averaged firing rates for each stimulus
fprintf('\nCalculating block-averaged firing rates...\n');

% Storage for block data
block_data = cell(3, 3);  % brain_regions × stimuli

for br = 1:3
    if isempty(results_all{br})
        continue;
    end

    responsive_idx = results_all{br}.responsive_idx;
    n_responsive = results_all{br}.n_responsive;

    for stim = 1:3
        % Skip CS for BA (only show US and CS+US)
        if br == 2 && stim == 1
            continue;
        end

        fprintf('  %s, %s...\n', brain_regions{br}, stim_names{stim});

        % Storage for this region-stimulus combination
        all_block_means = [];

        for n = 1:n_responsive
            neuron_idx = responsive_idx(n);

            % Get trial-by-trial spike times
            preAP = preAP_norm_all{stim}{neuron_idx};
            postAP = postAP_norm_all{stim}{neuron_idx};

            % Number of trials
            n_trials = length(preAP);

            % Calculate block size (trials per block for 10 blocks)
            block_size = floor(n_trials / 10);

            if block_size == 0
                continue;  % Skip if too few trials
            end

            % Calculate baseline firing rate for this neuron
            baseline_spikes = 0;
            for trial = 1:n_trials
                % Spikes in baseline window (-5 to 0s)
                spikes_baseline = preAP{trial}(preAP{trial} >= -g.pre_time & preAP{trial} < 0);
                baseline_spikes = baseline_spikes + length(spikes_baseline);
            end
            baseline_fr = baseline_spikes / (n_trials * g.pre_time);

            % Calculate firing rate for each block in response window (0 to 1s)
            block_means = zeros(1, 10);

            for block = 1:10
                trial_start = (block - 1) * block_size + 1;
                trial_end = block * block_size;

                if trial_end > n_trials
                    trial_end = n_trials;
                end

                % Count spikes in response window for these trials
                total_spikes = 0;
                for trial = trial_start:trial_end
                    % Spikes in response window (0 to 1s)
                    spikes_in_window = postAP{trial}(postAP{trial} >= 0 & postAP{trial} < g.test_time);
                    total_spikes = total_spikes + length(spikes_in_window);
                end

                % Calculate mean firing rate (Hz) and delta FR
                n_trials_in_block = trial_end - trial_start + 1;
                block_fr = total_spikes / (n_trials_in_block * g.test_time);
                block_means(block) = block_fr - baseline_fr;  % Delta FR
            end

            % Store this neuron's block means
            all_block_means = [all_block_means; block_means];
        end

        % Store mean and SEM across neurons
        block_data{br, stim}.mean = mean(all_block_means, 1);
        block_data{br, stim}.sem = std(all_block_means, 0, 1) / sqrt(size(all_block_means, 1));
        block_data{br, stim}.n_neurons = size(all_block_means, 1);

        fprintf('    %d neurons with sufficient trials\n', size(all_block_means, 1));
    end
end

%% Create figure
fig = figure('Units', 'pixels', 'Position', [100, 100, 1500, 500], 'Visible', 'on');
t = tiledlayout(fig, 1, 3, 'TileSpacing', 'tight', 'Padding', 'compact');

% X-axis (block numbers)
x_blocks = 1:10;

for stim = 1:3
    ax = nexttile(t, stim);
    hold on;

    % Plot each brain region - first shaded areas, then lines
    legend_handles = [];
    legend_labels = {};

    % First pass: plot shaded SEM areas
    for br = 1:3
        % Skip CS for BA (only show US and CS+US)
        if br == 2 && stim == 1
            continue;
        end
        if isempty(block_data{br, stim})
            continue;
        end

        data = block_data{br, stim};
        color = region_colors.(brain_regions{br});

        % Add shaded SEM
        fill([x_blocks, fliplr(x_blocks)], ...
             [data.mean + data.sem, fliplr(data.mean - data.sem)], ...
             color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    % Second pass: plot mean lines (for legend)
    for br = 1:3
        % Skip CS for BA (only show US and CS+US)
        if br == 2 && stim == 1
            continue;
        end
        if isempty(block_data{br, stim})
            continue;
        end

        data = block_data{br, stim};
        color = region_colors.(brain_regions{br});

        % Plot line
        h = plot(x_blocks, data.mean, '-', 'Color', color, 'LineWidth', 3);
        legend_handles = [legend_handles; h];
        legend_labels{end+1} = sprintf('%s (n=%d)', brain_regions{br}, block_data{br, stim}.n_neurons);
    end

    hold off;

    % Formatting
    xlabel('Block #', 'FontSize', g.fontSize2);

    if stim == 1
        ylabel('\DeltaFiring Rate (Hz)', 'FontSize', g.fontSize2);
    end

    title(stim_names{stim}, 'FontSize', g.fontSize1, 'FontWeight', 'bold');

    xlim([0.5 10.5]);
    ylim([0 30]);
    xticks(1:10);

    set(gca, 'FontSize', g.fontSize2);
    box off;

    % Add legend only on CS+US panel
    if stim == 3 && ~isempty(legend_handles)
        lg = legend(legend_handles, legend_labels, 'Location', 'best');
        lg.FontSize = g.fontSize2;
    end
end

fprintf('\nDone.\n');
