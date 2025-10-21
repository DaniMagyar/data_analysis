% BAfc_correlation_CS_US_CSandUS
% Calculate correlations between neuronal responses to CS, US, and CS+US stimuli
%
% This script:
% 1. Loads neuronal data from all recordings
% 2. Calculates Z-scored PSTHs for CS, US, and CS+US conditions
% 3. Averages Z-scored responses during first 1s after stimulus onset
% 4. Calculates correlations between stimulus pairs (CS vs US, CS vs CS+US, US vs CS+US)
% 5. Performs analysis separately for each brain region
% 6. Allows selection of cell types (PN, IN, or both)

clear all; close all

%% ===== PARAMETERS =====

% Cell type selection: 'PN', 'IN', or 'both'
g.cellType = 'PN';  % Change this to 'IN' or 'both' as needed

% Brain regions to analyze
g.brainRegions = {'LA', 'BA'};  % Add more regions if needed

% All recordings
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

% Stimulus conditions
ttl = {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'};
stim_names = {'CS', 'US', 'CS+US'};

% PSTH parameters
g.pre_time = 5;
g.post_time = 5;
g.bin_time = 0.01;
g.smoothvalue = 5;

% Response window for averaging (0 to 1s after stimulus onset)
g.response_window_start = 0;
g.response_window_end = 0.5;

% Visualization parameters
g.fontSize1 = 13;
g.fontSize2 = 12;
g.colors = BAfc_colors;

%% ===== LOAD DATA =====

fprintf('\n========== LOADING DATA ==========\n');
fprintf('Cell type: %s\n', g.cellType);
fprintf('Number of recordings: %d\n', length(recordings));

% Load cell metrics
g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
g.cell_metrics = BAfc_match_celltypes('cell_metrics', g.cell_metrics);

fprintf('Total neurons loaded: %d\n', length(g.cell_metrics.brainRegion));

%% ===== PROCESS EACH BRAIN REGION =====

results = struct();

for bR_idx = 1:length(g.brainRegions)
    bR = g.brainRegions{bR_idx};

    fprintf('\n========== PROCESSING BRAIN REGION: %s ==========\n', bR);

    %% Filter neurons by brain region and cell type

    % Brain region filter
    idx_region = strcmp(g.cell_metrics.brainRegion, bR);

    % Cell type filter
    switch g.cellType
        case 'PN'
            idx_cellType = strcmp(g.cell_metrics.putativeCellType, 'Pyramidal Cell') | ...
                          strcmp(g.cell_metrics.putativeCellType, 'PN');
        case 'IN'
            idx_cellType = strcmp(g.cell_metrics.putativeCellType, 'Narrow Interneuron') | ...
                          strcmp(g.cell_metrics.putativeCellType, 'Wide Interneuron') | ...
                          strcmp(g.cell_metrics.putativeCellType, 'IN');
        case 'both'
            idx_cellType = true(size(idx_region));
        otherwise
            error('Invalid cell type. Use ''PN'', ''IN'', or ''both''');
    end

    % Combined filter
    idx_neurons = idx_region & idx_cellType;
    n_neurons = sum(idx_neurons);

    fprintf('Neurons in %s (%s): %d\n', bR, g.cellType, n_neurons);

    if n_neurons < 2
        fprintf('WARNING: Not enough neurons for correlation analysis in %s. Skipping.\n', bR);
        continue;
    end

    %% Calculate Z-scored PSTHs and average responses

    % Storage for average responses
    avg_responses = zeros(n_neurons, length(ttl));

    % Storage for full z-scored PSTHs (for visualization)
    psth_zscore_all = cell(1, length(ttl));

    for stim_idx = 1:length(ttl)
        fprintf('Processing stimulus: %s (%s)...\n', stim_names{stim_idx}, ttl{stim_idx});

        % Get raw PSTH
        psth_raw = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{stim_idx}, ...
            'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);

        % Z-score each neuron (relative to entire trial)
        psth_z = zscore(psth_raw, 0, 2);

        % Smooth
        psth_z = smoothdata(psth_z, 2, 'sgolay', g.smoothvalue);

        % Filter by region and cell type
        psth_z_filtered = psth_z(idx_neurons, :);

        % Store full PSTH
        psth_zscore_all{stim_idx} = psth_z_filtered;

        % Define response window indices
        response_start_idx = round((g.pre_time + g.response_window_start) / g.bin_time) + 1;
        response_end_idx = round((g.pre_time + g.response_window_end) / g.bin_time);

        % Average Z-scored response in the response window (first 1s after onset)
        avg_responses(:, stim_idx) = mean(psth_z_filtered(:, response_start_idx:response_end_idx), 2);

        fprintf('  Mean response (Z-score): %.3f ± %.3f\n', ...
            mean(avg_responses(:, stim_idx)), std(avg_responses(:, stim_idx)));
    end

    %% Calculate correlations between stimulus pairs

    fprintf('\nCalculating correlations...\n');

    % Remove neurons with NaN values
    valid_idx = ~any(isnan(avg_responses), 2);
    avg_responses_clean = avg_responses(valid_idx, :);
    n_valid = sum(valid_idx);

    fprintf('Valid neurons for correlation: %d / %d\n', n_valid, n_neurons);

    if n_valid < 2
        fprintf('WARNING: Not enough valid neurons for correlation in %s. Skipping.\n', bR);
        continue;
    end

    % Calculate correlations
    % 1. CS vs US
    [r_CS_US, p_CS_US] = corrcoef(avg_responses_clean(:, 1), avg_responses_clean(:, 2));
    r_CS_US = r_CS_US(1, 2);
    p_CS_US = p_CS_US(1, 2);

    % 2. CS vs CS+US
    [r_CS_CSUS, p_CS_CSUS] = corrcoef(avg_responses_clean(:, 1), avg_responses_clean(:, 3));
    r_CS_CSUS = r_CS_CSUS(1, 2);
    p_CS_CSUS = p_CS_CSUS(1, 2);

    % 3. US vs CS+US
    [r_US_CSUS, p_US_CSUS] = corrcoef(avg_responses_clean(:, 2), avg_responses_clean(:, 3));
    r_US_CSUS = r_US_CSUS(1, 2);
    p_US_CSUS = p_US_CSUS(1, 2);

    % Display results
    fprintf('\nCorrelation Results for %s (%s):\n', bR, g.cellType);
    fprintf('%-15s: r = %6.3f, p = %.4f %s\n', 'CS vs US', r_CS_US, p_CS_US, ...
        get_sig_str(p_CS_US));
    fprintf('%-15s: r = %6.3f, p = %.4f %s\n', 'CS vs CS+US', r_CS_CSUS, p_CS_CSUS, ...
        get_sig_str(p_CS_CSUS));
    fprintf('%-15s: r = %6.3f, p = %.4f %s\n', 'US vs CS+US', r_US_CSUS, p_US_CSUS, ...
        get_sig_str(p_US_CSUS));

    % Store results
    results.(bR).n_neurons = n_neurons;
    results.(bR).n_valid = n_valid;
    results.(bR).avg_responses = avg_responses_clean;
    results.(bR).psth_zscore_all = psth_zscore_all;
    results.(bR).r_CS_US = r_CS_US;
    results.(bR).p_CS_US = p_CS_US;
    results.(bR).r_CS_CSUS = r_CS_CSUS;
    results.(bR).p_CS_CSUS = p_CS_CSUS;
    results.(bR).r_US_CSUS = r_US_CSUS;
    results.(bR).p_US_CSUS = p_US_CSUS;
end

%% ===== VISUALIZATION =====

fprintf('\n========== CREATING VISUALIZATIONS ==========\n');

n_regions = length(fieldnames(results));
if n_regions == 0
    fprintf('No results to visualize.\n');
    return;
end

% Figure 1: Scatter plots showing correlations
fig1 = figure('Position', [100, 100, 1800, 600]);
t1 = tiledlayout(fig1, n_regions, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

region_names = fieldnames(results);

for region_idx = 1:n_regions
    bR = region_names{region_idx};
    data = results.(bR);

    % Panel 1: CS vs US
    nexttile(t1, (region_idx-1)*3 + 1);
    scatter(data.avg_responses(:, 1), data.avg_responses(:, 2), ...
        60, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
    hold on;

    % Add regression line
    p = polyfit(data.avg_responses(:, 1), data.avg_responses(:, 2), 1);
    x_fit = linspace(min(data.avg_responses(:, 1)), max(data.avg_responses(:, 1)), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

    % Unity line
    lims = [min([data.avg_responses(:, 1); data.avg_responses(:, 2)]), ...
            max([data.avg_responses(:, 1); data.avg_responses(:, 2)])];
    plot(lims, lims, '--k', 'LineWidth', 1);

    xlabel('CS Response (Z-score)', 'FontSize', g.fontSize2);
    ylabel('US Response (Z-score)', 'FontSize', g.fontSize2);
    title(sprintf('%s: CS vs US\nr = %.3f, p = %.4f %s', ...
        bR, data.r_CS_US, data.p_CS_US, get_sig_str(data.p_CS_US)), ...
        'FontSize', g.fontSize1);
    axis square;
    grid on;
    set(gca, 'FontSize', g.fontSize2);

    % Panel 2: CS vs CS+US
    nexttile(t1, (region_idx-1)*3 + 2);
    scatter(data.avg_responses(:, 1), data.avg_responses(:, 3), ...
        60, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
    hold on;

    % Add regression line
    p = polyfit(data.avg_responses(:, 1), data.avg_responses(:, 3), 1);
    x_fit = linspace(min(data.avg_responses(:, 1)), max(data.avg_responses(:, 1)), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

    % Unity line
    lims = [min([data.avg_responses(:, 1); data.avg_responses(:, 3)]), ...
            max([data.avg_responses(:, 1); data.avg_responses(:, 3)])];
    plot(lims, lims, '--k', 'LineWidth', 1);

    xlabel('CS Response (Z-score)', 'FontSize', g.fontSize2);
    ylabel('CS+US Response (Z-score)', 'FontSize', g.fontSize2);
    title(sprintf('%s: CS vs CS+US\nr = %.3f, p = %.4f %s', ...
        bR, data.r_CS_CSUS, data.p_CS_CSUS, get_sig_str(data.p_CS_CSUS)), ...
        'FontSize', g.fontSize1);
    axis square;
    grid on;
    set(gca, 'FontSize', g.fontSize2);

    % Panel 3: US vs CS+US
    nexttile(t1, (region_idx-1)*3 + 3);
    scatter(data.avg_responses(:, 2), data.avg_responses(:, 3), ...
        60, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
    hold on;

    % Add regression line
    p = polyfit(data.avg_responses(:, 2), data.avg_responses(:, 3), 1);
    x_fit = linspace(min(data.avg_responses(:, 2)), max(data.avg_responses(:, 2)), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

    % Unity line
    lims = [min([data.avg_responses(:, 2); data.avg_responses(:, 3)]), ...
            max([data.avg_responses(:, 2); data.avg_responses(:, 3)])];
    plot(lims, lims, '--k', 'LineWidth', 1);

    xlabel('US Response (Z-score)', 'FontSize', g.fontSize2);
    ylabel('CS+US Response (Z-score)', 'FontSize', g.fontSize2);
    title(sprintf('%s: US vs CS+US\nr = %.3f, p = %.4f %s', ...
        bR, data.r_US_CSUS, data.p_US_CSUS, get_sig_str(data.p_US_CSUS)), ...
        'FontSize', g.fontSize1);
    axis square;
    grid on;
    set(gca, 'FontSize', g.fontSize2);
end

t1.Title.String = sprintf('Cross-Stimulus Response Correlations (%s neurons)', g.cellType);
t1.Title.FontSize = g.fontSize1 + 2;

% Figure 2: Bar plot comparing correlation coefficients
fig2 = figure('Position', [200, 200, 800, 600]);
t2 = tiledlayout(fig2, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile(t2);
hold on;

% Prepare data for bar plot
n_regions = length(region_names);
bar_data = zeros(n_regions, 3);
p_values = zeros(n_regions, 3);

for region_idx = 1:n_regions
    bR = region_names{region_idx};
    bar_data(region_idx, 1) = results.(bR).r_CS_US;
    bar_data(region_idx, 2) = results.(bR).r_CS_CSUS;
    bar_data(region_idx, 3) = results.(bR).r_US_CSUS;

    p_values(region_idx, 1) = results.(bR).p_CS_US;
    p_values(region_idx, 2) = results.(bR).p_CS_CSUS;
    p_values(region_idx, 3) = results.(bR).p_US_CSUS;
end

% Create grouped bar plot
x = 1:n_regions;
bar_width = 0.25;
colors_bar = [0.2 0.4 0.8; 0.8 0.4 0.2; 0.4 0.8 0.2];

for i = 1:3
    bar(x + (i-2)*bar_width, bar_data(:, i), bar_width, ...
        'FaceColor', colors_bar(i, :), 'EdgeColor', 'k', 'LineWidth', 1);
end

% Add significance markers
for region_idx = 1:n_regions
    for i = 1:3
        x_pos = region_idx + (i-2)*bar_width;
        y_pos = bar_data(region_idx, i);

        if y_pos < 0
            y_offset = -0.05;
            valign = 'top';
        else
            y_offset = 0.05;
            valign = 'bottom';
        end

        text(x_pos, y_pos + y_offset, get_sig_str(p_values(region_idx, i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', valign, ...
            'FontSize', 11, 'FontWeight', 'bold');
    end
end

% Formatting
yline(0, '-k', 'LineWidth', 1);
xlim([0.5, n_regions + 0.5]);
xticks(1:n_regions);
xticklabels(region_names);
ylabel('Correlation Coefficient (r)', 'FontSize', g.fontSize1);
xlabel('Brain Region', 'FontSize', g.fontSize1);
title(sprintf('Cross-Stimulus Correlation Comparison (%s neurons)', g.cellType), ...
    'FontSize', g.fontSize1 + 1);
legend({'CS vs US', 'CS vs CS+US', 'US vs CS+US'}, 'Location', 'best', ...
    'FontSize', g.fontSize2);
grid on;
set(gca, 'FontSize', g.fontSize2);

% Figure 3: Distribution histograms
fig3 = figure('Position', [300, 300, 1800, 600]);
t3 = tiledlayout(fig3, n_regions, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for region_idx = 1:n_regions
    bR = region_names{region_idx};
    data = results.(bR);

    % Histogram 1: CS responses
    nexttile(t3, (region_idx-1)*3 + 1);
    histogram(data.avg_responses(:, 1), 20, 'FaceColor', [0.3 0.3 0.8], ...
        'EdgeColor', 'k', 'FaceAlpha', 0.7);
    xline(0, '--k', 'LineWidth', 2);
    xlabel('CS Response (Z-score)', 'FontSize', g.fontSize2);
    ylabel('Count', 'FontSize', g.fontSize2);
    title(sprintf('%s: CS\nMean = %.2f', bR, mean(data.avg_responses(:, 1))), ...
        'FontSize', g.fontSize1);
    grid on;
    set(gca, 'FontSize', g.fontSize2);

    % Histogram 2: US responses
    nexttile(t3, (region_idx-1)*3 + 2);
    histogram(data.avg_responses(:, 2), 20, 'FaceColor', [0.8 0.3 0.3], ...
        'EdgeColor', 'k', 'FaceAlpha', 0.7);
    xline(0, '--k', 'LineWidth', 2);
    xlabel('US Response (Z-score)', 'FontSize', g.fontSize2);
    ylabel('Count', 'FontSize', g.fontSize2);
    title(sprintf('%s: US\nMean = %.2f', bR, mean(data.avg_responses(:, 2))), ...
        'FontSize', g.fontSize1);
    grid on;
    set(gca, 'FontSize', g.fontSize2);

    % Histogram 3: CS+US responses
    nexttile(t3, (region_idx-1)*3 + 3);
    histogram(data.avg_responses(:, 3), 20, 'FaceColor', [0.3 0.8 0.3], ...
        'EdgeColor', 'k', 'FaceAlpha', 0.7);
    xline(0, '--k', 'LineWidth', 2);
    xlabel('CS+US Response (Z-score)', 'FontSize', g.fontSize2);
    ylabel('Count', 'FontSize', g.fontSize2);
    title(sprintf('%s: CS+US\nMean = %.2f', bR, mean(data.avg_responses(:, 3))), ...
        'FontSize', g.fontSize1);
    grid on;
    set(gca, 'FontSize', g.fontSize2);
end

t3.Title.String = sprintf('Response Distribution by Stimulus Type (%s neurons)', g.cellType);
t3.Title.FontSize = g.fontSize1 + 2;

%% ===== SUMMARY TABLE =====

fprintf('\n========== SUMMARY TABLE ==========\n');
fprintf('Cell Type: %s\n\n', g.cellType);
fprintf('%-10s | %-8s | %-15s | %-15s | %-15s\n', ...
    'Region', 'N Cells', 'CS vs US', 'CS vs CS+US', 'US vs CS+US');
fprintf('%s\n', repmat('-', 1, 85));

for region_idx = 1:n_regions
    bR = region_names{region_idx};
    data = results.(bR);

    fprintf('%-10s | %-8d | r=%.3f p=%.4f | r=%.3f p=%.4f | r=%.3f p=%.4f\n', ...
        bR, data.n_valid, ...
        data.r_CS_US, data.p_CS_US, ...
        data.r_CS_CSUS, data.p_CS_CSUS, ...
        data.r_US_CSUS, data.p_US_CSUS);
end

fprintf('\n========== ANALYSIS COMPLETE ==========\n');

%% ===== HELPER FUNCTION =====

function sig_str = get_sig_str(p_val)
    % Returns significance string based on p-value
    if p_val < 0.001
        sig_str = '***';
    elseif p_val < 0.01
        sig_str = '**';
    elseif p_val < 0.05
        sig_str = '*';
    else
        sig_str = 'n.s.';
    end
end
