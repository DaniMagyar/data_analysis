% BAfc_supplementary_figure_2_v2
% Visualize stability of neuronal responses across trial blocks
% Uses cluster assignments from figure_2_data.xlsx

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

ttl = {'triptest_sound_only','triptest_shocks_only'};
stim_names = {'CS', 'US'};

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

% Brain regions to analyze
brain_regions = {'LA', 'BA', 'Astria'};

% Region colors
region_colors = struct();
region_colors.LA = [0.7 0.2 0.2];
region_colors.BA = [0.2 0.7 0.2];
region_colors.Astria = [0.2 0.4 0.7];

%% Load cluster assignments from figure_2_data.xlsx
fprintf('Loading cluster assignments from figure_2_data.xlsx...\n');
figure_2_file = '..\figure_2\figure_2_data.xlsx';

% Read the raw data sheet
raw_data = readcell(figure_2_file, 'Sheet', 'RawData_AllNeurons');

% Find header row - search for it
header_row = [];
for r = 1:min(10, size(raw_data, 1))
    if any(strcmp(raw_data(r, :), 'Global Index'))
        header_row = r;
        break;
    end
end

if isempty(header_row)
    error('Could not find header row with Global Index in figure_2_data.xlsx');
end

global_idx_col = find(strcmp(raw_data(header_row, :), 'Global Index'));
cluster_col = find(strcmp(raw_data(header_row, :), 'Cluster'));

if isempty(global_idx_col) || isempty(cluster_col)
    error('Could not find Global Index or Cluster columns in figure_2_data.xlsx');
end

% Extract data
neuron_clusters = struct();
n_loaded = 0;
for row = (header_row+1):size(raw_data, 1)
    % Check if we have data in this row
    if row <= size(raw_data, 1) && global_idx_col <= size(raw_data, 2)
        global_idx_val = raw_data{row, global_idx_col};

        if ~isempty(global_idx_val) && isnumeric(global_idx_val)
            global_idx = global_idx_val;
            cluster_name = raw_data{row, cluster_col};

            % Store cluster assignment
            field_name = sprintf('n%d', global_idx);
            neuron_clusters.(field_name).cluster = cluster_name;
            neuron_clusters.(field_name).global_idx = global_idx;
            n_loaded = n_loaded + 1;
        end
    end
end

fprintf('  Loaded %d neuron cluster assignments\n', length(fieldnames(neuron_clusters)));

%% Get trial-by-trial spike data
fprintf('Extracting trial-by-trial spike data...\n');
preAP_norm_all = cell(1, numel(ttl));
postAP_norm_all = cell(1, numel(ttl));

for hmp = 1:numel(ttl)
    [~, preAP_norm_all{hmp}, postAP_norm_all{hmp}] = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
end

%% Organize neurons by region and cluster
fprintf('\nOrganizing neurons by region and cluster...\n');

% Fixed block size: 5 trials per block
trials_per_block = 5;
max_trials = 50;  % Only use first 50 trials for all neurons

% Storage: brain_regions × stimuli
results_all = cell(3, 3);

for br = 1:3
    fprintf('  Processing %s...\n', brain_regions{br});

    % Get all neurons in this region
    idx_neurons = strcmp(g.cell_metrics.brainRegion, brain_regions{br});
    neuron_indices = find(idx_neurons);

    % Filter neurons by stimulus type based on cluster
    for stim = 1:2  % Only CS and US
        % Skip CS for BA
        if br == 2 && stim == 1
            continue;
        end

        responsive_indices = [];

        for n = 1:length(neuron_indices)
            global_idx = neuron_indices(n);
            field_name = sprintf('n%d', global_idx);
            if isfield(neuron_clusters, field_name)
                cluster_name = neuron_clusters.(field_name).cluster;

                % Assign neurons to appropriate stimuli based on cluster
                include_neuron = false;
                if stim == 1  % CS
                    % CS-sel and Multi respond to CS
                    if ismember(cluster_name, {'CS-sel', 'Multi'})
                        include_neuron = true;
                    end
                elseif stim == 2  % US
                    % US-sel and Multi respond to US
                    if ismember(cluster_name, {'US-sel', 'Multi'})
                        include_neuron = true;
                    end
                end

                if include_neuron
                    responsive_indices = [responsive_indices; global_idx];
                end
            end
        end

        results_all{br, stim}.responsive_idx = responsive_indices;
        results_all{br, stim}.n_responsive = length(responsive_indices);

        fprintf('    %s: %d neurons\n', stim_names{stim}, length(responsive_indices));
    end
end

%% Calculate maximum blocks
fprintf('\nDetermining maximum number of blocks...\n');
max_blocks = floor(max_trials / trials_per_block);
fprintf('  Using first %d trials for all neurons (%d blocks)\n', max_trials, max_blocks);

%% Calculate block-averaged firing rates
fprintf('\nCalculating block-averaged firing rates...\n');

block_data = cell(3, 2);  % brain_regions × stimuli (CS, US only)

for br = 1:3
    for stim = 1:2  % Only CS and US
        % Skip CS for BA
        if br == 2 && stim == 1
            continue;
        end

        if isempty(results_all{br, stim})
            continue;
        end

        fprintf('  %s, %s...\n', brain_regions{br}, stim_names{stim});

        responsive_idx = results_all{br, stim}.responsive_idx;
        n_responsive = results_all{br, stim}.n_responsive;

        % Storage for this region-stimulus combination
        all_block_means = nan(n_responsive, max_blocks);
        neuron_global_idx = zeros(n_responsive, 1);

        for n = 1:n_responsive
            neuron_idx = responsive_idx(n);
            neuron_global_idx(n) = neuron_idx;

            % Get trial-by-trial spike times
            preAP = preAP_norm_all{stim}{neuron_idx};
            postAP = postAP_norm_all{stim}{neuron_idx};

            % Number of trials (limit to max_trials)
            n_trials_total = length(preAP);
            n_trials = min(n_trials_total, max_trials);

            % Skip if neuron has fewer than 5 trials
            if n_trials < trials_per_block
                continue;
            end

            % Calculate how many blocks this neuron can contribute (capped at max_blocks)
            n_blocks_this_neuron = floor(n_trials / trials_per_block);

            % Calculate firing rate for each block (5 trials per block)
            for block = 1:n_blocks_this_neuron
                trial_start = (block - 1) * trials_per_block + 1;
                trial_end = block * trials_per_block;

                % Count spikes in response window for these trials
                total_spikes = 0;
                for trial = trial_start:trial_end
                    % Spikes in response window (0 to 1s)
                    spikes_in_window = postAP{trial}(postAP{trial} >= 0 & postAP{trial} < g.test_time);
                    total_spikes = total_spikes + length(spikes_in_window);
                end

                % Calculate mean firing rate (Hz)
                block_fr = total_spikes / (trials_per_block * g.test_time);
                all_block_means(n, block) = block_fr;  % Simple FR, not delta
            end
        end

        % Calculate mean and SEM across neurons, ignoring NaN values
        block_data{br, stim}.mean = nanmean(all_block_means, 1);
        block_data{br, stim}.sem = nanstd(all_block_means, 0, 1) ./ sqrt(sum(~isnan(all_block_means), 1));
        block_data{br, stim}.n_neurons_per_block = sum(~isnan(all_block_means), 1);
        block_data{br, stim}.n_neurons = sum(any(~isnan(all_block_means), 2));

        % Store which neurons contributed to each block
        block_data{br, stim}.neurons_per_block = cell(1, max_blocks);
        for block = 1:max_blocks
            contributing_neurons = neuron_global_idx(~isnan(all_block_means(:, block)));
            block_data{br, stim}.neurons_per_block{block} = contributing_neurons;
        end

        fprintf('    %d neurons included\n', block_data{br, stim}.n_neurons);
    end
end

%% Create figure
fig = figure('Position', [100, 100, 1000, 500], 'Units', 'pixels');
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add main title
title(t, 'Response stability across trial blocks', 'FontSize', 12, 'FontWeight', 'bold');

% X-axis (block numbers)
x_blocks = 1:max_blocks;

for stim = 1:2  % Only CS and US
    % Create nested tiledlayout for each stimulus
    t_nest = tiledlayout(t, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    t_nest.Layout.Tile = stim;

    % Add stimulus-specific title
    if stim == 1
        title(t_nest, 'CS responses', 'FontSize', 12, 'FontWeight', 'bold');
    else
        title(t_nest, 'US responses', 'FontSize', 12, 'FontWeight', 'bold');
    end

    ax = nexttile(t_nest, 1);
    hold on;

    % Plot each brain region - first shaded areas, then lines
    legend_handles = [];
    legend_labels = {};

    % First pass: plot shaded SEM areas
    for br = 1:3
        % Skip CS for BA
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
        % Skip CS for BA
        if br == 2 && stim == 1
            continue;
        end
        if isempty(block_data{br, stim})
            continue;
        end

        data = block_data{br, stim};
        color = region_colors.(brain_regions{br});

        % Plot line
        h = plot(x_blocks, data.mean, '-', 'Color', color, 'LineWidth', 2);
        legend_handles = [legend_handles; h];
        legend_labels{end+1} = sprintf('%s (n=%d)', brain_regions{br}, block_data{br, stim}.n_neurons);
    end

    hold off;

    % Formatting
    xlabel('Block #', 'FontSize', 10);

    if stim == 1
        ylabel('Firing Rate (Hz)', 'FontSize', 10);
    end

    xlim([0.5 max_blocks+0.5]);
    ylim([0 30]);
    xticks(1:max_blocks);

    set(gca, 'FontSize', 10);
    box off;

    % Add legend only on US panel
    if stim == 2 && ~isempty(legend_handles)
        lg = legend(legend_handles, legend_labels, 'Location', 'best');
        lg.FontSize = 10;
    end
end

%% Export to Excel
export_supplementary_figure_2_to_excel(block_data, results_all, brain_regions, stim_names, max_blocks, trials_per_block, max_trials, g);
exportgraphics(gcf, 'supplementary_figure_2.png', 'Resolution', 300);
fprintf('\nDone.\n');

%% Export function
function export_supplementary_figure_2_to_excel(block_data, results_all, brain_regions, stim_names, max_blocks, trials_per_block, max_trials, g)
    output_filename = 'supplementary_figure_2_data.xlsx';

    if exist(output_filename, 'file')
        delete(output_filename);
    end

    %% Sheet 1: Summary
    summary_data = {};
    summary_data{1,1} = 'Supplementary Figure 2: Response Stability Across Trial Blocks';
    summary_data{2,1} = '';
    summary_data{3,1} = 'Cluster assignments loaded from figure_2_data.xlsx';
    summary_data{4,1} = '';
    summary_data{5,1} = 'Analysis Parameters';
    summary_data{6,1} = 'Trials per block';
    summary_data{6,2} = trials_per_block;
    summary_data{7,1} = 'Max trials used';
    summary_data{7,2} = max_trials;
    summary_data{8,1} = 'Number of blocks';
    summary_data{8,2} = max_blocks;
    summary_data{9,1} = 'Response window';
    summary_data{9,2} = '0-1s';
    summary_data{10,1} = 'Baseline window';
    summary_data{10,2} = '-5-0s';
    summary_data{11,1} = '';
    summary_data{12,1} = 'Neuron Counts by Region';
    summary_data{13,1} = 'Region';
    summary_data{13,2} = 'CS';
    summary_data{13,3} = 'US';

    row = 14;
    for br = 1:3
        summary_data{row,1} = brain_regions{br};
        for stim = 1:2  % Only CS and US
            if br == 2 && stim == 1
                summary_data{row, stim+1} = 'N/A';
            elseif ~isempty(block_data{br, stim})
                summary_data{row, stim+1} = block_data{br, stim}.n_neurons;
            else
                summary_data{row, stim+1} = 0;
            end
        end
        row = row + 1;
    end

    writecell(summary_data, output_filename, 'Sheet', 'Summary');

    %% Sheet 2: Block Statistics
    stats_data = {};
    stats_data{1,1} = 'Block-averaged Firing Rate Statistics';
    stats_data{2,1} = '';
    stats_data{3,1} = 'Mean ± SEM (Hz) for each block (response window 0-1s)';
    stats_data{4,1} = '';

    row = 5;
    for br = 1:3
        for stim = 1:2  % Only CS and US
            if br == 2 && stim == 1
                continue;
            end
            if isempty(block_data{br, stim})
                continue;
            end

            stats_data{row,1} = sprintf('%s %s', brain_regions{br}, stim_names{stim});
            stats_data{row+1,1} = 'Block';
            stats_data{row+2,1} = 'Mean (Hz)';
            stats_data{row+3,1} = 'SEM (Hz)';
            stats_data{row+4,1} = 'N neurons';

            for block = 1:max_blocks
                stats_data{row+1, block+1} = block;
                if block <= length(block_data{br, stim}.mean)
                    stats_data{row+2, block+1} = block_data{br, stim}.mean(block);
                    stats_data{row+3, block+1} = block_data{br, stim}.sem(block);
                    stats_data{row+4, block+1} = block_data{br, stim}.n_neurons_per_block(block);
                else
                    stats_data{row+2, block+1} = NaN;
                    stats_data{row+3, block+1} = NaN;
                    stats_data{row+4, block+1} = 0;
                end
            end

            row = row + 6;
        end
    end

    writecell(stats_data, output_filename, 'Sheet', 'Block_Statistics');

    %% Sheet 3: Neurons Per Block (Global IDs)
    neurons_per_block_data = {};
    neurons_per_block_data{1,1} = 'Neurons Contributing to Each Block';
    neurons_per_block_data{2,1} = '';
    neurons_per_block_data{3,1} = 'Region-Stimulus';
    neurons_per_block_data{3,2} = 'Block';
    neurons_per_block_data{3,3} = 'Global Index';

    row = 4;
    for br = 1:3
        for stim = 1:2  % Only CS and US
            if br == 2 && stim == 1
                continue;
            end
            if isempty(block_data{br, stim})
                continue;
            end

            for block = 1:max_blocks
                neuron_ids = block_data{br, stim}.neurons_per_block{block};
                if ~isempty(neuron_ids)
                    for idx = 1:length(neuron_ids)
                        neurons_per_block_data{row, 1} = sprintf('%s %s', brain_regions{br}, stim_names{stim});
                        neurons_per_block_data{row, 2} = block;
                        neurons_per_block_data{row, 3} = neuron_ids(idx);
                        row = row + 1;
                    end
                end
            end
        end
    end

    writecell(neurons_per_block_data, output_filename, 'Sheet', 'Neurons_Per_Block');

    %% Sheet 4: Neuron List
    neuron_data = {};
    neuron_data{1,1} = 'Neurons Included in Analysis';
    neuron_data{2,1} = '';
    neuron_data{3,1} = 'Global Index';
    neuron_data{3,2} = 'Animal ID';
    neuron_data{3,3} = 'Cell ID';
    neuron_data{3,4} = 'Region';
    neuron_data{3,5} = 'Cell Type';

    row = 4;
    for br = 1:3
        for stim = 1:2  % Only CS and US
            if isempty(results_all{br, stim})
                continue;
            end

            responsive_idx = results_all{br, stim}.responsive_idx;

            for n = 1:length(responsive_idx)
                neuron_idx = responsive_idx(n);

                % Check if already added
                already_added = false;
                for check_row = 4:(row-1)
                    if ~isempty(neuron_data{check_row, 1}) && neuron_data{check_row, 1} == neuron_idx
                        already_added = true;
                        break;
                    end
                end

                if ~already_added
                    neuron_data{row, 1} = neuron_idx;
                    neuron_data{row, 2} = g.cell_metrics.animal{neuron_idx};
                    neuron_data{row, 3} = g.cell_metrics.cellID(neuron_idx);
                    neuron_data{row, 4} = brain_regions{br};
                    neuron_data{row, 5} = g.cell_metrics.putativeCellType{neuron_idx};

                    row = row + 1;
                end
            end
        end
    end

    writecell(neuron_data, output_filename, 'Sheet', 'Neuron_List');

end
