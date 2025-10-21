% BAfc_figure_monosynaptic.m
clear all
%% ===== USER CONFIGURATION =====

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
    'MD307_001_kilosort',...
    'MD309_001_kilosort',...
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

% Stimulus selection
ttl = {'triptest_sound_only','triptest_shocks_only', 'triptest_both'};

%% ===== LOAD DATA =====

fprintf('\nLoading neural data...\n');
g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
g.cell_metrics = BAfc_match_celltypes('cell_metrics', g.cell_metrics);

%% ===== PARAMETERS =====

fprintf('\nLoading default parameters...\n');

% Get default parameters (modify BAfc_optogenetics_params.m to change defaults)
params = BAfc_optogenetics_params();
g.params = params;
% Override specific parameters here if needed for this figure
% params.zscore_threshold = 8;  % Example: use lower threshold
% params.bR = 'LA';             % Example: only analyze LA neurons

fprintf('\nIdentifying responsive neurons...\n');
responsive_results = BAfc_identify_responsive_neurons(g.cell_metrics, ttl, params);

% BAfc_monosyn_raster_ui(g, responsive_results, ttl);

%% ===== FIGURE: MONOSYNAPTIC RESPONSES BY STIMULUS TYPE =====

% Define brain regions and cell types to plot
brain_regions = {'LA', 'BA'};
cell_types = {'PN', 'IN'};
stimulus_names = {'CS (Sound)', 'US (Shock)', 'CS+US (Both)'};

% Create figure with 3 rows (one per stimulus) x 3 columns (magnitude, probability, latency)
fig = figure('Position', [50 50 1600 1200]);

for stim_idx = 1:length(ttl)
    ttl_fn = matlab.lang.makeValidName(ttl{stim_idx});

    % Get all responsive neurons for this stimulus
    responsive_idx = responsive_results.(ttl_fn).responsive_idx;

    if isempty(responsive_idx)
        fprintf('No responsive neurons found for %s\n', ttl{stim_idx});
        continue;
    end

    % Extract data for responsive neurons only
    all_brain_regions = g.cell_metrics.brainRegion(responsive_results.(ttl_fn).all_neuron_idx(responsive_idx));
    all_cell_types = g.cell_metrics.putativeCellType(responsive_results.(ttl_fn).all_neuron_idx(responsive_idx));
    magnitude = responsive_results.(ttl_fn).peak_zscore;
    probability = responsive_results.(ttl_fn).response_probability;
    latency = responsive_results.(ttl_fn).mean_latency;

    % --- SUBPLOT 1: Response Magnitude (Z-score) ---
    ax1 = subplot(3, 3, (stim_idx-1)*3 + 1);
    hold on;

    % Plot each group separately
    group_positions = [];
    group_labels = {};
    pos = 1;

    for br_idx = 1:length(brain_regions)
        for ct_idx = 1:length(cell_types)
            br = brain_regions{br_idx};
            ct = cell_types{ct_idx};

            % Find neurons matching this brain region and cell type
            idx = strcmp(all_brain_regions, br) & strcmp(all_cell_types, ct);

            if sum(idx) > 0
                data = magnitude(idx);

                % Plot individual points with jitter
                x_jitter = pos + 0.3 * (rand(sum(idx), 1) - 0.5);
                scatter(x_jitter, data, 20, 'o', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

                % Plot mean and SD
                mean_val = mean(data);
                sd_val = std(data);
                errorbar(pos, mean_val, sd_val, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'LineWidth', 2);

                group_positions(end+1) = pos;
                group_labels{end+1} = sprintf('%s-%s (n=%d)', br, ct, sum(idx));
                pos = pos + 1;
            end
        end
    end

    xlim([0.5 pos-0.5]);
    xticks(group_positions);
    xticklabels(group_labels);
    ylabel('Peak Z-score');
    title(sprintf('%s - Response Magnitude', stimulus_names{stim_idx}));
    grid on;
    hold off;

    % --- SUBPLOT 2: Response Probability ---
    ax2 = subplot(3, 3, (stim_idx-1)*3 + 2);
    hold on;

    pos = 1;
    group_positions = [];
    group_labels = {};

    for br_idx = 1:length(brain_regions)
        for ct_idx = 1:length(cell_types)
            br = brain_regions{br_idx};
            ct = cell_types{ct_idx};

            % Find neurons matching this brain region and cell type
            idx = strcmp(all_brain_regions, br) & strcmp(all_cell_types, ct);

            if sum(idx) > 0
                data = probability(idx);

                % Plot individual points with jitter
                x_jitter = pos + 0.3 * (rand(sum(idx), 1) - 0.5);
                scatter(x_jitter, data, 20, 'o', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

                % Plot mean and SD
                mean_val = mean(data);
                sd_val = std(data);
                errorbar(pos, mean_val, sd_val, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'LineWidth', 2);

                group_positions(end+1) = pos;
                group_labels{end+1} = sprintf('%s-%s (n=%d)', br, ct, sum(idx));
                pos = pos + 1;
            end
        end
    end

    xlim([0.5 pos-0.5]);
    ylim([0 1]);
    xticks(group_positions);
    xticklabels(group_labels);
    ylabel('Response Probability');
    title(sprintf('%s - Response Probability', stimulus_names{stim_idx}));
    grid on;
    hold off;

    % --- SUBPLOT 3: Response Latency ---
    ax3 = subplot(3, 3, (stim_idx-1)*3 + 3);
    hold on;

    pos = 1;
    group_positions = [];
    group_labels = {};

    for br_idx = 1:length(brain_regions)
        for ct_idx = 1:length(cell_types)
            br = brain_regions{br_idx};
            ct = cell_types{ct_idx};

            % Find neurons matching this brain region and cell type
            idx = strcmp(all_brain_regions, br) & strcmp(all_cell_types, ct);

            if sum(idx) > 0
                data = latency(idx);

                % Plot individual points with jitter
                x_jitter = pos + 0.3 * (rand(sum(idx), 1) - 0.5);
                scatter(x_jitter, data, 20, 'o', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

                % Plot mean and SD
                mean_val = mean(data);
                sd_val = std(data);
                errorbar(pos, mean_val, sd_val, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'LineWidth', 2);

                group_positions(end+1) = pos;
                group_labels{end+1} = sprintf('%s-%s (n=%d)', br, ct, sum(idx));
                pos = pos + 1;
            end
        end
    end

    xlim([0.5 pos-0.5]);
    xticks(group_positions);
    xticklabels(group_labels);
    ylabel('Mean Latency (ms)');
    title(sprintf('%s - Response Latency', stimulus_names{stim_idx}));
    grid on;
    hold off;
end

fprintf('\nFigure complete!\n');
