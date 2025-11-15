% BAfc_figure_2_example_rasters.m

% NOTE: Only first 50 trials plotted per condition

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

g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
g.cell_metrics = BAfc_match_celltypes('cell_metrics', g.cell_metrics);
g.colors = BAfc_colors;
g.fontSize1 = 10;
g.fontSize2 = 10;
g.xlinewidth = 2;
g.pre_time = 5;
g.post_time = 4;
g.bin_time = 0.001;

% Time window for rasters - same as heatmap
g.plotwin = [2 2];  % -2s to +2s

% Define example neurons: [animal, cellID, title]
example_neurons = {
    'MD307_001', 14, 'Example CS-selective neuron';
    'MD296_001', 32, 'Example US-selective neuron';
    'MD296_001', 31, 'Example Multisensory neuron';
    'MD296_001', 11, 'Example Inhibited neuron'
};

%% Create figure
fig = figure('Position', [100, 100, 1000, 300], 'Units', 'pixels');
t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Add panel label A
annotation(fig, 'textbox', [0.01 0.95 0.05 0.05], 'String', 'A', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Get spike times for both conditions once
fprintf('Loading spike times...\n');
[~, preAP_norm_cs, postAP_norm_cs] = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{1}, ...
    'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
[~, preAP_norm_us, postAP_norm_us] = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{2}, ...
    'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

%% Plot rasters
for ex = 1:4
    animal_name = example_neurons{ex, 1};
    target_cellID = example_neurons{ex, 2};
    neuron_title = example_neurons{ex, 3};

    % Find neuron index
    animal_match = strcmp(g.cell_metrics.animal, animal_name);
    cellID_match = g.cell_metrics.cellID == target_cellID;
    neuron_idx = find(animal_match & cellID_match);

    if isempty(neuron_idx)
        fprintf('Warning: Neuron %s/%d not found\n', animal_name, target_cellID);
        continue;
    end

    fprintf('Plotting %s/%d...\n', animal_name, target_cellID);

    % Calculate tile positions for nested layout
    % Main layout is 2x2, so tiles are: 1,2 (row 1), 3,4 (row 2)
    % ex=1 -> tile 1, ex=2 -> tile 2, ex=3 -> tile 3, ex=4 -> tile 4
    main_tile = ex;

    % Create nested tiledlayout for this neuron (1 row, 2 columns: CS and US)
    t_nested = tiledlayout(t, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');
    t_nested.Layout.Tile = main_tile;

    % Add title to nested layout
    title(t_nested, neuron_title, 'FontSize', g.fontSize1, 'FontWeight', 'bold');

    % CS raster (left column for this neuron)
    ax_cs = nexttile(t_nested, 1);
    preAP_cs = preAP_norm_cs{neuron_idx};
    postAP_cs = postAP_norm_cs{neuron_idx};

    if ~isempty(preAP_cs) || ~isempty(postAP_cs)
        hold(ax_cs, 'on');
        n_trials_total = length(preAP_cs);
        n_trials = min(n_trials_total, 50);  % Plot only first 50 trials
        for trial = 1:n_trials
            % Plot pre-stimulus spikes (negative times)
            if ~isempty(preAP_cs{trial})
                plot(ax_cs, preAP_cs{trial}, trial*ones(size(preAP_cs{trial})), 'k.', 'MarkerSize', 6);
            end
            % Plot post-stimulus spikes (positive times)
            if ~isempty(postAP_cs{trial})
                plot(ax_cs, postAP_cs{trial}, trial*ones(size(postAP_cs{trial})), 'k.', 'MarkerSize', 6);
            end
        end
        % Stimulus line
        plot(ax_cs, [0 0], [0 n_trials+1], 'r-', 'LineWidth', g.xlinewidth);
        hold(ax_cs, 'off');

        xlim(ax_cs, [-g.plotwin(1) g.plotwin(2)]);
        ylim(ax_cs, [0 n_trials+1]);

        % Y-label for CS (left)
        ylabel(ax_cs, 'Trial #', 'FontSize', g.fontSize2);

        % X-label only for bottom row
        if ex > 2
            xlabel(ax_cs, 'Time (s)', 'FontSize', g.fontSize2);
        else
            set(ax_cs, 'XTickLabel', []);
        end

        title(ax_cs, 'CS', 'FontSize', g.fontSize2);
        set(ax_cs, 'FontSize', g.fontSize2);
        box(ax_cs, 'off');
    end

    % US raster (right column for this neuron)
    ax_us = nexttile(t_nested, 2);
    preAP_us = preAP_norm_us{neuron_idx};
    postAP_us = postAP_norm_us{neuron_idx};

    if ~isempty(preAP_us) || ~isempty(postAP_us)
        hold(ax_us, 'on');
        n_trials_total = length(preAP_us);
        n_trials = min(n_trials_total, 50);  % Plot only first 50 trials
        for trial = 1:n_trials
            % Plot pre-stimulus spikes (negative times)
            if ~isempty(preAP_us{trial})
                plot(ax_us, preAP_us{trial}, trial*ones(size(preAP_us{trial})), 'k.', 'MarkerSize', 6);
            end
            % Plot post-stimulus spikes (positive times)
            if ~isempty(postAP_us{trial})
                plot(ax_us, postAP_us{trial}, trial*ones(size(postAP_us{trial})), 'k.', 'MarkerSize', 6);
            end
        end
        % Stimulus line
        plot(ax_us, [0 0], [0 n_trials+1], 'r-', 'LineWidth', g.xlinewidth);
        hold(ax_us, 'off');

        xlim(ax_us, [-g.plotwin(1) g.plotwin(2)]);
        ylim(ax_us, [0 n_trials+1]);

        % No Y-labels for US (right)
        set(ax_us, 'YTickLabel', []);

        % X-label only for bottom row
        if ex > 2
            xlabel(ax_us, 'Time (s)', 'FontSize', g.fontSize2);
        else
            set(ax_us, 'XTickLabel', []);
        end

        title(ax_us, 'US', 'FontSize', g.fontSize2);
        set(ax_us, 'FontSize', g.fontSize2);
        box(ax_us, 'off');
    end
end

fprintf('\nDone.\n');
