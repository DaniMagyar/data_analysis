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

%% ===== FIGURE: MONOSYNAPTIC RESPONSE HEATMAPS =====

% Define brain regions and cell types to plot
brain_regions = {'LA'};
cell_types = {'PN', 'IN'};
stimulus_names = {'CS (Sound)', 'US (Shock)', 'CS+US (Both)'};

% Heatmap parameters
g.bin_time = params.bin_time;
g.pre_time = 5;
g.post_time = 5;
g.plotwin = [0.1 0.1];  % -0.1 to 0.1s window
g.timeaxis_hmp = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.clim = [-5.5 13];
g.smoothvalue = 5;
g.colors = BAfc_colors;
g.fontSize1 = 13;
g.fontSize2 = 12;
g.xlinewidth = 2;

% Calculate PSTHs ONCE for all stimuli (much faster)
fprintf('\nCalculating PSTHs for all stimuli...\n');
psth_all = cell(1, length(ttl));
for stim_idx = 1:length(ttl)
    fprintf('  %s...\n', ttl{stim_idx});
    psth_spx = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{stim_idx}, ...
        'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    psth_spx = zscore(psth_spx, 0, 2);
    psth_spx = smoothdata(psth_spx, 2, 'sgolay', g.smoothvalue);
    psth_all{stim_idx} = psth_spx;
end

% First, identify all neurons responsive to ANY stimulus for each group
fprintf('\nIdentifying neurons responsive to ANY stimulus...\n');
group_neurons = struct();

group_idx = 1;
for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        br = brain_regions{br_idx};
        ct = cell_types{ct_idx};

        % Collect all responsive neurons across all stimuli
        all_responsive_neurons = [];

        for stim_idx = 1:length(ttl)
            ttl_fn = matlab.lang.makeValidName(ttl{stim_idx});

            % Get responsive neurons for this stimulus
            responsive_idx = responsive_results.(ttl_fn).responsive_idx;

            if ~isempty(responsive_idx)
                % Get all neuron indices and filter by brain region and cell type
                all_neuron_idx = responsive_results.(ttl_fn).all_neuron_idx(responsive_idx);
                all_brain_regions = g.cell_metrics.brainRegion(all_neuron_idx);
                all_cell_types = g.cell_metrics.putativeCellType(all_neuron_idx);

                % Find neurons matching this brain region and cell type
                group_mask = strcmp(all_brain_regions, br) & strcmp(all_cell_types, ct);
                group_neuron_idx = all_neuron_idx(group_mask);

                all_responsive_neurons = [all_responsive_neurons; group_neuron_idx(:)];
            end
        end

        % Get unique neurons
        all_responsive_neurons = unique(all_responsive_neurons);

        if isempty(all_responsive_neurons)
            fprintf('  %s %s: No responsive neurons\n', br, ct);
            group_neurons(group_idx).neurons = [];
            group_neurons(group_idx).br = br;
            group_neurons(group_idx).ct = ct;
            group_idx = group_idx + 1;
            continue;
        end

        % Categorize neurons by CS and US responsiveness
        ttl_fn_cs = matlab.lang.makeValidName(ttl{1});
        ttl_fn_us = matlab.lang.makeValidName(ttl{2});
        ttl_fn_csus = matlab.lang.makeValidName(ttl{3});

        cs_responsive_neurons = responsive_results.(ttl_fn_cs).all_neuron_idx(responsive_results.(ttl_fn_cs).responsive_idx);
        us_responsive_neurons = responsive_results.(ttl_fn_us).all_neuron_idx(responsive_results.(ttl_fn_us).responsive_idx);
        csus_responsive_neurons = responsive_results.(ttl_fn_csus).all_neuron_idx(responsive_results.(ttl_fn_csus).responsive_idx);

        % Get latencies for CS, US, and CS+US
        cs_latencies = containers.Map('KeyType', 'double', 'ValueType', 'double');
        us_latencies = containers.Map('KeyType', 'double', 'ValueType', 'double');
        csus_latencies = containers.Map('KeyType', 'double', 'ValueType', 'double');

        for i = 1:length(responsive_results.(ttl_fn_cs).responsive_idx)
            neuron_idx = cs_responsive_neurons(i);
            cs_latencies(neuron_idx) = responsive_results.(ttl_fn_cs).mean_latency(i);
        end

        for i = 1:length(responsive_results.(ttl_fn_us).responsive_idx)
            neuron_idx = us_responsive_neurons(i);
            us_latencies(neuron_idx) = responsive_results.(ttl_fn_us).mean_latency(i);
        end

        for i = 1:length(responsive_results.(ttl_fn_csus).responsive_idx)
            neuron_idx = csus_responsive_neurons(i);
            csus_latencies(neuron_idx) = responsive_results.(ttl_fn_csus).mean_latency(i);
        end

        % Categorize each neuron
        category_1 = [];  % CS only
        category_2 = [];  % US only
        category_3 = [];  % Both CS and US

        latency_1 = [];
        latency_2 = [];
        latency_3 = [];

        for i = 1:length(all_responsive_neurons)
            neuron_idx = all_responsive_neurons(i);
            is_cs = ismember(neuron_idx, cs_responsive_neurons);
            is_us = ismember(neuron_idx, us_responsive_neurons);

            if is_cs && ~is_us
                % Category 1: CS only
                category_1 = [category_1; neuron_idx];
                latency_1 = [latency_1; cs_latencies(neuron_idx)];
            elseif ~is_cs && is_us
                % Category 2: US only
                category_2 = [category_2; neuron_idx];
                latency_2 = [latency_2; us_latencies(neuron_idx)];
            elseif is_cs && is_us
                % Category 3: Both CS and US - sort by CS+US latency
                category_3 = [category_3; neuron_idx];
                % Check if neuron is responsive to CS+US, otherwise use CS latency as fallback
                if isKey(csus_latencies, neuron_idx)
                    latency_3 = [latency_3; csus_latencies(neuron_idx)];
                else
                    latency_3 = [latency_3; cs_latencies(neuron_idx)];
                end
            end
        end

        % Sort within each category by latency (ascending - earliest first)
        sorted_neurons = [];
        if ~isempty(category_1)
            [~, sort_idx] = sort(latency_1, 'ascend');
            sorted_neurons = [sorted_neurons; category_1(sort_idx)];
        end
        if ~isempty(category_2)
            [~, sort_idx] = sort(latency_2, 'ascend');
            sorted_neurons = [sorted_neurons; category_2(sort_idx)];
        end
        if ~isempty(category_3)
            [~, sort_idx] = sort(latency_3, 'ascend');
            sorted_neurons = [sorted_neurons; category_3(sort_idx)];
        end

        % Store sorted neuron indices
        group_neurons(group_idx).neurons = sorted_neurons;
        group_neurons(group_idx).br = br;
        group_neurons(group_idx).ct = ct;

        fprintf('  %s %s: %d neurons (CS only: %d, US only: %d, Both: %d)\n', ...
            br, ct, length(all_responsive_neurons), length(category_1), length(category_2), length(category_3));

        group_idx = group_idx + 1;
    end
end

% Create figure: 4 groups x 3 stimuli
fig = figure('Position', [50 50 1400 1000]);
t = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Now plot heatmaps with consistent neuron ordering
fprintf('\nGenerating heatmaps...\n');
for group_idx = 1:length(group_neurons)
    br = group_neurons(group_idx).br;
    ct = group_neurons(group_idx).ct;
    neuron_list = group_neurons(group_idx).neurons;

    if isempty(neuron_list)
        % No neurons for this group - create empty tiles
        for stim_idx = 1:length(ttl)
            nexttile(t, (group_idx-1)*3 + stim_idx);
            text(0.5, 0.5, sprintf('No %s %s neurons', br, ct), 'HorizontalAlignment', 'center');
            axis off;
        end
        continue;
    end

    fprintf('%s %s: %d neurons\n', br, ct, length(neuron_list));

    % Loop through each stimulus
    for stim_idx = 1:length(ttl)
        % Use pre-calculated PSTH
        psth_spx = psth_all{stim_idx};

        % Extract data for the group neurons (already sorted)
        psth_group = psth_spx(neuron_list, :);

        % Extract plot window (-0.1 to 0.1s)
        plot_bins = round((g.pre_time - g.plotwin(1))/g.bin_time + 1 : (g.pre_time + g.plotwin(2))/g.bin_time);
        matrix = psth_group(:, plot_bins);

        % Create heatmap
        ax = nexttile(t, (group_idx-1)*3 + stim_idx);
        imagesc(g.timeaxis_hmp, 1:size(matrix,1), matrix);
        clim(g.clim);
        colormap(gca, g.colors.Heatmap);
        xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha', 1);

        % Labels
        if stim_idx == 1
            ylabel(sprintf('%s %s (n=%d)', br, ct, size(matrix,1)), 'FontSize', g.fontSize2);
        else
            set(gca, 'YTickLabel', []);
        end

        if group_idx == 1
            title(stimulus_names{stim_idx}, 'FontSize', g.fontSize1);
        end

        if group_idx == 4
            xlabel('Time from Stimulus Onset (s)', 'FontSize', g.fontSize2);
        end

        set(gca, 'FontSize', g.fontSize2);
    end
end

% Add colorbar
cb = colorbar(ax, 'eastoutside', 'FontSize', g.fontSize2);
ylabel(cb, 'Z-score Firing Rate', 'FontSize', g.fontSize2);
cb.Layout.Tile = 'east';

% Overall title
t.Title.String = 'Monosynaptically Responsive Neurons (-0.1 to 0.1s)';
t.Title.FontSize = g.fontSize1 + 2;

fprintf('\nHeatmap figure complete!\n');

%% ===== QUANTITATIVE ANALYSES =====

fprintf('\n========================================\n');
fprintf('QUANTITATIVE ANALYSES\n');
fprintf('========================================\n');

%% 1. CS-US CONVERGENCE METRICS

fprintf('\n--- 1. CS-US CONVERGENCE METRICS ---\n');

ttl_fn_cs = matlab.lang.makeValidName(ttl{1});
ttl_fn_us = matlab.lang.makeValidName(ttl{2});
ttl_fn_csus = matlab.lang.makeValidName(ttl{3});

cs_responsive = responsive_results.(ttl_fn_cs).all_neuron_idx(responsive_results.(ttl_fn_cs).responsive_idx);
us_responsive = responsive_results.(ttl_fn_us).all_neuron_idx(responsive_results.(ttl_fn_us).responsive_idx);
csus_responsive = responsive_results.(ttl_fn_csus).all_neuron_idx(responsive_results.(ttl_fn_csus).responsive_idx);

convergence_stats = struct();

for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        br = brain_regions{br_idx};
        ct = cell_types{ct_idx};
        group_name = sprintf('%s_%s', br, ct);

        % Get all neurons in this group
        all_idx = find(strcmp(g.cell_metrics.brainRegion, br) & strcmp(g.cell_metrics.putativeCellType, ct));

        % Find responsive neurons
        cs_group = intersect(cs_responsive, all_idx);
        us_group = intersect(us_responsive, all_idx);
        csus_group = intersect(csus_responsive, all_idx);

        % Calculate convergence
        cs_only = setdiff(cs_group, us_group);
        us_only = setdiff(us_group, cs_group);
        both = intersect(cs_group, us_group);

        convergence_stats.(group_name).cs_only = cs_only;
        convergence_stats.(group_name).us_only = us_only;
        convergence_stats.(group_name).both = both;
        convergence_stats.(group_name).total = union(cs_group, us_group);

        if ~isempty(convergence_stats.(group_name).total)
            pct_cs_only = 100 * length(cs_only) / length(convergence_stats.(group_name).total);
            pct_us_only = 100 * length(us_only) / length(convergence_stats.(group_name).total);
            pct_both = 100 * length(both) / length(convergence_stats.(group_name).total);

            fprintf('\n%s %s (n=%d responsive):\n', br, ct, length(convergence_stats.(group_name).total));
            fprintf('  CS only:  %d (%.1f%%)\n', length(cs_only), pct_cs_only);
            fprintf('  US only:  %d (%.1f%%)\n', length(us_only), pct_us_only);
            fprintf('  Both:     %d (%.1f%%)\n', length(both), pct_both);
            fprintf('  Convergence index: %.3f\n', length(both) / length(convergence_stats.(group_name).total));
        end
    end
end

%% 2. INTEGRATION TYPE ANALYSIS (Supralinear/Linear/Sublinear)

fprintf('\n--- 2. INTEGRATION TYPE ANALYSIS ---\n');
fprintf('(For dual-responsive neurons: compare CS+US response to CS+US prediction)\n');

integration_stats = struct();

for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        br = brain_regions{br_idx};
        ct = cell_types{ct_idx};
        group_name = sprintf('%s_%s', br, ct);

        dual_neurons = convergence_stats.(group_name).both;

        if isempty(dual_neurons)
            continue;
        end

        % Get peak z-scores for CS, US, CS+US (0-50ms window)
        peak_window = round((g.pre_time)/g.bin_time + 1 : (g.pre_time + 0.05)/g.bin_time);

        cs_peaks = [];
        us_peaks = [];
        csus_peaks = [];

        for n = 1:length(dual_neurons)
            neuron_idx = dual_neurons(n);
            cs_peaks(n) = max(psth_all{1}(neuron_idx, peak_window));
            us_peaks(n) = max(psth_all{2}(neuron_idx, peak_window));
            csus_peaks(n) = max(psth_all{3}(neuron_idx, peak_window));
        end

        % Calculate integration index
        predicted_sum = cs_peaks + us_peaks;
        integration_index = (csus_peaks - predicted_sum) ./ predicted_sum;

        % Classify integration type (threshold = 20%)
        supralinear = sum(integration_index > 0.2);
        linear = sum(abs(integration_index) <= 0.2);
        sublinear = sum(integration_index < -0.2);

        integration_stats.(group_name).integration_index = integration_index;
        integration_stats.(group_name).supralinear = supralinear;
        integration_stats.(group_name).linear = linear;
        integration_stats.(group_name).sublinear = sublinear;

        fprintf('\n%s %s (n=%d dual-responsive):\n', br, ct, length(dual_neurons));
        fprintf('  Supralinear (>20%% enhancement): %d (%.1f%%)\n', supralinear, 100*supralinear/length(dual_neurons));
        fprintf('  Linear (±20%%):                  %d (%.1f%%)\n', linear, 100*linear/length(dual_neurons));
        fprintf('  Sublinear (<-20%% suppression):  %d (%.1f%%)\n', sublinear, 100*sublinear/length(dual_neurons));
        fprintf('  Mean integration index: %.3f ± %.3f\n', mean(integration_index), std(integration_index)/sqrt(length(integration_index)));
    end
end

%% 3. TEMPORAL PRECISION & RELIABILITY

fprintf('\n--- 3. TEMPORAL PRECISION & RELIABILITY ---\n');
fprintf('(Using trial-to-trial latency data from responsive_results)\n');

temporal_stats = struct();

for stim_idx = 1:length(ttl)
    ttl_fn = matlab.lang.makeValidName(ttl{stim_idx});

    fprintf('\n%s:\n', stimulus_names{stim_idx});

    for br_idx = 1:length(brain_regions)
        for ct_idx = 1:length(cell_types)
            br = brain_regions{br_idx};
            ct = cell_types{ct_idx};
            group_name = sprintf('%s_%s', br, ct);

            % Get responsive neurons for this group
            all_neuron_idx = responsive_results.(ttl_fn).all_neuron_idx(responsive_results.(ttl_fn).responsive_idx);
            all_brain_regions = g.cell_metrics.brainRegion(all_neuron_idx);
            all_cell_types = g.cell_metrics.putativeCellType(all_neuron_idx);

            group_mask = strcmp(all_brain_regions, br) & strcmp(all_cell_types, ct);

            if sum(group_mask) == 0
                continue;
            end

            % Get latency data
            mean_latencies = responsive_results.(ttl_fn).mean_latency(group_mask);

            % Note: std_latency not available in responsive_results
            % Using latency variability across population as proxy for temporal precision
            population_latency_std = std(mean_latencies);

            % Response probability
            response_prob = responsive_results.(ttl_fn).response_probability(group_mask);
            mean_prob = mean(response_prob);

            % Store stats
            if ~isfield(temporal_stats, group_name)
                temporal_stats.(group_name) = struct();
            end

            temporal_stats.(group_name).(ttl_fn).mean_latency = mean(mean_latencies);
            temporal_stats.(group_name).(ttl_fn).std_latency = std(mean_latencies);
            temporal_stats.(group_name).(ttl_fn).population_std = population_latency_std;
            temporal_stats.(group_name).(ttl_fn).response_prob = mean_prob;
            temporal_stats.(group_name).(ttl_fn).n = sum(group_mask);

            fprintf('  %s %s (n=%d):\n', br, ct, sum(group_mask));
            fprintf('    Latency: %.2f ± %.2f ms\n', 1000*mean(mean_latencies), 1000*std(mean_latencies));
            fprintf('    Population latency SD: %.2f ms\n', 1000*population_latency_std);
            fprintf('    Response probability: %.3f\n', mean_prob);
        end
    end
end

%% 4. LATENCY COMPARISONS

fprintf('\n--- 4. LATENCY COMPARISONS ---\n');

latency_stats = struct();

for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        br = brain_regions{br_idx};
        ct = cell_types{ct_idx};
        group_name = sprintf('%s_%s', br, ct);

        dual_neurons = convergence_stats.(group_name).both;

        if isempty(dual_neurons)
            continue;
        end

        % Get latencies for CS and US in dual-responsive neurons
        cs_latencies_dual = [];
        us_latencies_dual = [];

        for n = 1:length(dual_neurons)
            neuron_idx = dual_neurons(n);

            % Find this neuron in CS responsive results
            cs_resp_idx = find(responsive_results.(ttl_fn_cs).all_neuron_idx(responsive_results.(ttl_fn_cs).responsive_idx) == neuron_idx);
            if ~isempty(cs_resp_idx)
                cs_latencies_dual(n) = responsive_results.(ttl_fn_cs).mean_latency(cs_resp_idx);
            else
                cs_latencies_dual(n) = NaN;
            end

            % Find this neuron in US responsive results
            us_resp_idx = find(responsive_results.(ttl_fn_us).all_neuron_idx(responsive_results.(ttl_fn_us).responsive_idx) == neuron_idx);
            if ~isempty(us_resp_idx)
                us_latencies_dual(n) = responsive_results.(ttl_fn_us).mean_latency(us_resp_idx);
            else
                us_latencies_dual(n) = NaN;
            end
        end

        % Remove NaN values
        valid_idx = ~isnan(cs_latencies_dual) & ~isnan(us_latencies_dual);
        cs_latencies_dual = cs_latencies_dual(valid_idx);
        us_latencies_dual = us_latencies_dual(valid_idx);

        if isempty(cs_latencies_dual)
            continue;
        end

        % Calculate latency difference (CS - US)
        latency_diff = cs_latencies_dual - us_latencies_dual;

        % Statistical test
        [h, p] = ttest(cs_latencies_dual, us_latencies_dual);

        latency_stats.(group_name).cs_latencies = cs_latencies_dual;
        latency_stats.(group_name).us_latencies = us_latencies_dual;
        latency_stats.(group_name).latency_diff = latency_diff;
        latency_stats.(group_name).ttest_p = p;

        fprintf('\n%s %s (n=%d dual-responsive):\n', br, ct, length(cs_latencies_dual));
        fprintf('  CS latency: %.2f ± %.2f ms\n', 1000*mean(cs_latencies_dual), 1000*std(cs_latencies_dual)/sqrt(length(cs_latencies_dual)));
        fprintf('  US latency: %.2f ± %.2f ms\n', 1000*mean(us_latencies_dual), 1000*std(us_latencies_dual)/sqrt(length(us_latencies_dual)));
        fprintf('  Difference (CS-US): %.2f ± %.2f ms\n', 1000*mean(latency_diff), 1000*std(latency_diff)/sqrt(length(latency_diff)));
        if p < 0.001
            fprintf('  Paired t-test: p < 0.001 ***\n');
        elseif p < 0.01
            fprintf('  Paired t-test: p = %.3f **\n', p);
        elseif p < 0.05
            fprintf('  Paired t-test: p = %.3f *\n', p);
        else
            fprintf('  Paired t-test: p = %.3f (n.s.)\n', p);
        end

        % Count neurons with faster CS vs US response
        faster_cs = sum(cs_latencies_dual < us_latencies_dual);
        faster_us = sum(us_latencies_dual < cs_latencies_dual);
        fprintf('  Faster CS response: %d (%.1f%%)\n', faster_cs, 100*faster_cs/length(cs_latencies_dual));
        fprintf('  Faster US response: %d (%.1f%%)\n', faster_us, 100*faster_us/length(cs_latencies_dual));
    end
end

%% 5. RESPONSE MAGNITUDE ANALYSIS

fprintf('\n--- 5. RESPONSE MAGNITUDE ANALYSIS ---\n');

magnitude_stats = struct();
peak_window = round((g.pre_time)/g.bin_time + 1 : (g.pre_time + 0.05)/g.bin_time);

for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        br = brain_regions{br_idx};
        ct = cell_types{ct_idx};
        group_name = sprintf('%s_%s', br, ct);

        dual_neurons = convergence_stats.(group_name).both;

        if isempty(dual_neurons)
            continue;
        end

        % Get peak z-scores for CS, US, CS+US
        cs_peaks = [];
        us_peaks = [];
        csus_peaks = [];

        for n = 1:length(dual_neurons)
            neuron_idx = dual_neurons(n);
            cs_peaks(n) = max(psth_all{1}(neuron_idx, peak_window));
            us_peaks(n) = max(psth_all{2}(neuron_idx, peak_window));
            csus_peaks(n) = max(psth_all{3}(neuron_idx, peak_window));
        end

        % Statistical comparisons
        [~, p_cs_vs_us] = ttest(cs_peaks, us_peaks);
        [~, p_csus_vs_cs] = ttest(csus_peaks, cs_peaks);
        [~, p_csus_vs_us] = ttest(csus_peaks, us_peaks);

        % Summation analysis
        predicted_sum = cs_peaks + us_peaks;
        [~, p_summation] = ttest(csus_peaks, predicted_sum);

        magnitude_stats.(group_name).cs_peaks = cs_peaks;
        magnitude_stats.(group_name).us_peaks = us_peaks;
        magnitude_stats.(group_name).csus_peaks = csus_peaks;

        fprintf('\n%s %s (n=%d dual-responsive):\n', br, ct, length(dual_neurons));
        fprintf('  CS peak z-score:     %.2f ± %.2f\n', mean(cs_peaks), std(cs_peaks)/sqrt(length(cs_peaks)));
        fprintf('  US peak z-score:     %.2f ± %.2f\n', mean(us_peaks), std(us_peaks)/sqrt(length(us_peaks)));
        fprintf('  CS+US peak z-score:  %.2f ± %.2f\n', mean(csus_peaks), std(csus_peaks)/sqrt(length(csus_peaks)));
        fprintf('  CS vs US: p = %.4f %s\n', p_cs_vs_us, get_sig_str(p_cs_vs_us));
        fprintf('  CS+US vs CS: p = %.4f %s\n', p_csus_vs_cs, get_sig_str(p_csus_vs_cs));
        fprintf('  CS+US vs US: p = %.4f %s\n', p_csus_vs_us, get_sig_str(p_csus_vs_us));
        fprintf('  CS+US vs (CS+US) predicted sum: p = %.4f %s\n', p_summation, get_sig_str(p_summation));
    end
end

%% 6. SELECTIVITY INDEX

fprintf('\n--- 6. SELECTIVITY INDEX ---\n');
fprintf('(Positive = CS-selective, Negative = US-selective)\n');

selectivity_stats = struct();

for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        br = brain_regions{br_idx};
        ct = cell_types{ct_idx};
        group_name = sprintf('%s_%s', br, ct);

        dual_neurons = convergence_stats.(group_name).both;

        if isempty(dual_neurons)
            continue;
        end

        cs_peaks = magnitude_stats.(group_name).cs_peaks;
        us_peaks = magnitude_stats.(group_name).us_peaks;

        % Selectivity index: (CS - US) / (CS + US)
        selectivity_index = (cs_peaks - us_peaks) ./ (cs_peaks + us_peaks);

        % Count selective neurons (threshold = 0.2)
        cs_selective = sum(selectivity_index > 0.2);
        us_selective = sum(selectivity_index < -0.2);
        non_selective = sum(abs(selectivity_index) <= 0.2);

        selectivity_stats.(group_name).selectivity_index = selectivity_index;

        fprintf('\n%s %s (n=%d dual-responsive):\n', br, ct, length(dual_neurons));
        fprintf('  CS-selective (>0.2):  %d (%.1f%%)\n', cs_selective, 100*cs_selective/length(dual_neurons));
        fprintf('  Non-selective (±0.2): %d (%.1f%%)\n', non_selective, 100*non_selective/length(dual_neurons));
        fprintf('  US-selective (<-0.2): %d (%.1f%%)\n', us_selective, 100*us_selective/length(dual_neurons));
        fprintf('  Mean selectivity index: %.3f ± %.3f\n', mean(selectivity_index), std(selectivity_index)/sqrt(length(selectivity_index)));
    end
end

%% 7. POPULATION LATENCY DISTRIBUTIONS

fprintf('\n--- 7. POPULATION LATENCY DISTRIBUTIONS ---\n');

% Create figure for latency distributions
fig_latency = figure('Position', [100 100 1200 800]);
t_lat = tiledlayout(fig_latency, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for stim_idx = 1:length(ttl)
    ttl_fn = matlab.lang.makeValidName(ttl{stim_idx});

    % Collect all latencies across groups
    all_latencies = struct();
    for br_idx = 1:length(brain_regions)
        for ct_idx = 1:length(cell_types)
            br = brain_regions{br_idx};
            ct = cell_types{ct_idx};
            group_name = sprintf('%s_%s', br, ct);

            if isfield(temporal_stats, group_name) && isfield(temporal_stats.(group_name), ttl_fn)
                % Get latencies for this group
                all_neuron_idx = responsive_results.(ttl_fn).all_neuron_idx(responsive_results.(ttl_fn).responsive_idx);
                all_brain_regions = g.cell_metrics.brainRegion(all_neuron_idx);
                all_cell_types = g.cell_metrics.putativeCellType(all_neuron_idx);

                group_mask = strcmp(all_brain_regions, br) & strcmp(all_cell_types, ct);
                latencies = responsive_results.(ttl_fn).mean_latency(group_mask);

                all_latencies.(group_name) = latencies;  % Already in ms from responsive_results
            end
        end
    end

    % Plot histogram
    ax = nexttile(t_lat, stim_idx);
    hold on;

    colors_plot = struct('LA_PN', [0 0.4470 0.7410], 'LA_IN', [0.8500 0.3250 0.0980], ...
                        'BA_PN', [0.4940 0.1840 0.5560], 'BA_IN', [0.4660 0.6740 0.1880]);

    for br_idx = 1:length(brain_regions)
        for ct_idx = 1:length(cell_types)
            group_name = sprintf('%s_%s', brain_regions{br_idx}, cell_types{ct_idx});
            if isfield(all_latencies, group_name) && ~isempty(all_latencies.(group_name))
                histogram(all_latencies.(group_name), 0:1:50, 'FaceColor', colors_plot.(group_name), ...
                         'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', strrep(group_name, '_', ' '));
            end
        end
    end

    xlabel('Latency (ms)', 'FontSize', 12);
    ylabel('Count', 'FontSize', 12);
    title(stimulus_names{stim_idx}, 'FontSize', 13);
    if stim_idx == 1
        legend('Location', 'northeast', 'FontSize', 10);
    end
    xlim([0 50]);
    grid on;
    hold off;
end

% Plot latency comparison for dual-responsive neurons
for stim_idx = 1:2
    ax = nexttile(t_lat, 3 + stim_idx);
    hold on;

    if stim_idx == 1
        % CS vs US latency scatter
        for br_idx = 1:length(brain_regions)
            for ct_idx = 1:length(cell_types)
                group_name = sprintf('%s_%s', brain_regions{br_idx}, cell_types{ct_idx});
                if isfield(latency_stats, group_name)
                    cs_lat = latency_stats.(group_name).cs_latencies;  % Already in ms
                    us_lat = latency_stats.(group_name).us_latencies;  % Already in ms
                    scatter(cs_lat, us_lat, 40, colors_plot.(group_name), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
        plot([0 50], [0 50], 'k--', 'LineWidth', 1.5);
        xlabel('CS Latency (ms)', 'FontSize', 12);
        ylabel('US Latency (ms)', 'FontSize', 12);
        title('CS vs US Latency (Dual-responsive)', 'FontSize', 13);
        xlim([0 50]);
        ylim([0 50]);
        axis square;
        grid on;
    else
        % Latency difference by group - use scatter plot instead of boxplot
        group_labels = {};
        latency_diffs_all = [];
        group_positions = [];
        group_colors_all = [];

        pos = 1;
        for br_idx = 1:length(brain_regions)
            for ct_idx = 1:length(cell_types)
                group_name = sprintf('%s_%s', brain_regions{br_idx}, cell_types{ct_idx});
                if isfield(latency_stats, group_name)
                    group_labels{end+1} = strrep(group_name, '_', ' ');
                    diffs = latency_stats.(group_name).latency_diff;  % Already in ms

                    % Add jitter to x-axis for visibility
                    x_jitter = pos + 0.2*(rand(size(diffs))-0.5);

                    scatter(x_jitter, diffs, 40, colors_plot.(group_name), 'filled', 'MarkerFaceAlpha', 0.6);

                    % Add mean and SEM
                    errorbar(pos, mean(diffs), std(diffs)/sqrt(length(diffs)), ...
                            'k', 'LineWidth', 2, 'CapSize', 10);

                    pos = pos + 1;
                end
            end
        end
        plot([0 pos], [0 0], 'k--', 'LineWidth', 1.5);
        xlim([0.5 pos-0.5]);
        xticks(1:length(group_labels));
        xticklabels(group_labels);
        xtickangle(45);
        ylabel('Latency Difference CS-US (ms)', 'FontSize', 12);
        title('Latency Difference (Dual-responsive)', 'FontSize', 13);
        grid on;
    end
    hold off;
end

t_lat.Title.String = 'Population Latency Analysis';
t_lat.Title.FontSize = 15;

fprintf('Population latency distribution figure created.\n');

%% FIGURE: CONVERGENCE AND INTEGRATION ANALYSIS

fprintf('\n--- Creating convergence and integration figure ---\n');

fig_conv = figure('Position', [150 150 1400 900]);
t_conv = tiledlayout(fig_conv, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. Convergence proportions (stacked bar)
ax1 = nexttile(t_conv, 1);
group_names_plot = {};
cs_only_pct = [];
us_only_pct = [];
both_pct = [];

for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        group_name = sprintf('%s_%s', brain_regions{br_idx}, cell_types{ct_idx});
        if isfield(convergence_stats, group_name) && ~isempty(convergence_stats.(group_name).total)
            group_names_plot{end+1} = sprintf('%s %s', brain_regions{br_idx}, cell_types{ct_idx});
            total = length(convergence_stats.(group_name).total);
            cs_only_pct(end+1) = 100 * length(convergence_stats.(group_name).cs_only) / total;
            us_only_pct(end+1) = 100 * length(convergence_stats.(group_name).us_only) / total;
            both_pct(end+1) = 100 * length(convergence_stats.(group_name).both) / total;
        end
    end
end

if ~isempty(group_names_plot)
    b = bar(1:length(group_names_plot), [cs_only_pct; us_only_pct; both_pct]', 'stacked');
    b(1).FaceColor = [0.3 0.6 0.9];
    b(2).FaceColor = [0.9 0.5 0.3];
    b(3).FaceColor = [0.5 0.8 0.5];
    xticks(1:length(group_names_plot));
    xticklabels(group_names_plot);
    xtickangle(45);
    ylabel('Percentage (%)', 'FontSize', 11);
    title('CS-US Convergence', 'FontSize', 12);
    legend({'CS only', 'US only', 'Both'}, 'Location', 'northwest', 'FontSize', 10);
    ylim([0 100]);
    grid on;
end

% 2. Integration type (supralinear/linear/sublinear)
ax2 = nexttile(t_conv, 2);
group_names_int = {};
supralin_pct = [];
linear_pct = [];
sublin_pct = [];

for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        group_name = sprintf('%s_%s', brain_regions{br_idx}, cell_types{ct_idx});
        if isfield(integration_stats, group_name)
            group_names_int{end+1} = sprintf('%s %s', brain_regions{br_idx}, cell_types{ct_idx});
            n_total = length(integration_stats.(group_name).integration_index);
            supralin_pct(end+1) = 100 * integration_stats.(group_name).supralinear / n_total;
            linear_pct(end+1) = 100 * integration_stats.(group_name).linear / n_total;
            sublin_pct(end+1) = 100 * integration_stats.(group_name).sublinear / n_total;
        end
    end
end

if ~isempty(group_names_int)
    b = bar(1:length(group_names_int), [supralin_pct; linear_pct; sublin_pct]', 'stacked');
    b(1).FaceColor = [0.8 0.3 0.3];
    b(2).FaceColor = [0.7 0.7 0.7];
    b(3).FaceColor = [0.3 0.3 0.8];
    xticks(1:length(group_names_int));
    xticklabels(group_names_int);
    xtickangle(45);
    ylabel('Percentage (%)', 'FontSize', 11);
    title('Integration Type (Dual-responsive)', 'FontSize', 12);
    legend({'Supralinear', 'Linear', 'Sublinear'}, 'Location', 'northwest', 'FontSize', 10);
    ylim([0 100]);
    grid on;
end

% 3. Selectivity index distribution
ax3 = nexttile(t_conv, 3);
group_names_sel = {};
cs_sel_pct = [];
non_sel_pct = [];
us_sel_pct = [];

for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        group_name = sprintf('%s_%s', brain_regions{br_idx}, cell_types{ct_idx});
        if isfield(selectivity_stats, group_name)
            group_names_sel{end+1} = sprintf('%s %s', brain_regions{br_idx}, cell_types{ct_idx});
            sel_idx = selectivity_stats.(group_name).selectivity_index;
            n_total = length(sel_idx);
            cs_sel_pct(end+1) = 100 * sum(sel_idx > 0.2) / n_total;
            non_sel_pct(end+1) = 100 * sum(abs(sel_idx) <= 0.2) / n_total;
            us_sel_pct(end+1) = 100 * sum(sel_idx < -0.2) / n_total;
        end
    end
end

if ~isempty(group_names_sel)
    b = bar(1:length(group_names_sel), [cs_sel_pct; non_sel_pct; us_sel_pct]', 'stacked');
    b(1).FaceColor = [0.3 0.6 0.9];
    b(2).FaceColor = [0.7 0.7 0.7];
    b(3).FaceColor = [0.9 0.5 0.3];
    xticks(1:length(group_names_sel));
    xticklabels(group_names_sel);
    xtickangle(45);
    ylabel('Percentage (%)', 'FontSize', 11);
    title('Selectivity Index (Dual-responsive)', 'FontSize', 12);
    legend({'CS-selective', 'Non-selective', 'US-selective'}, 'Location', 'northwest', 'FontSize', 10);
    ylim([0 100]);
    grid on;
end

% 4. Response magnitude comparison (CS vs US vs CS+US)
ax4 = nexttile(t_conv, 4, [1 2]);
hold on;
pos = 1;
group_names_mag = {};
colors_mag = [0.3 0.6 0.9; 0.9 0.5 0.3; 0.5 0.8 0.5];

for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        group_name = sprintf('%s_%s', brain_regions{br_idx}, cell_types{ct_idx});
        if isfield(magnitude_stats, group_name)
            group_names_mag{end+1} = sprintf('%s %s', brain_regions{br_idx}, cell_types{ct_idx});

            cs_pk = magnitude_stats.(group_name).cs_peaks;
            us_pk = magnitude_stats.(group_name).us_peaks;
            csus_pk = magnitude_stats.(group_name).csus_peaks;

            % Plot individual points with jitter
            x_jitter = 0.15;
            scatter(pos-0.25 + x_jitter*(rand(size(cs_pk))-0.5), cs_pk, 20, colors_mag(1,:), 'filled', 'MarkerFaceAlpha', 0.4);
            scatter(pos + x_jitter*(rand(size(us_pk))-0.5), us_pk, 20, colors_mag(2,:), 'filled', 'MarkerFaceAlpha', 0.4);
            scatter(pos+0.25 + x_jitter*(rand(size(csus_pk))-0.5), csus_pk, 20, colors_mag(3,:), 'filled', 'MarkerFaceAlpha', 0.4);

            % Plot means with error bars
            errorbar(pos-0.25, mean(cs_pk), std(cs_pk)/sqrt(length(cs_pk)), 'Color', colors_mag(1,:), 'LineWidth', 2, 'CapSize', 8);
            errorbar(pos, mean(us_pk), std(us_pk)/sqrt(length(us_pk)), 'Color', colors_mag(2,:), 'LineWidth', 2, 'CapSize', 8);
            errorbar(pos+0.25, mean(csus_pk), std(csus_pk)/sqrt(length(csus_pk)), 'Color', colors_mag(3,:), 'LineWidth', 2, 'CapSize', 8);

            pos = pos + 1;
        end
    end
end
hold off;
xlim([0.5 pos-0.5]);
xticks(1:length(group_names_mag));
xticklabels(group_names_mag);
xtickangle(45);
ylabel('Peak Z-score', 'FontSize', 11);
title('Response Magnitude (Dual-responsive)', 'FontSize', 12);
legend({'CS', 'US', 'CS+US'}, 'Location', 'northwest', 'FontSize', 10);
grid on;

% 5. Integration index scatter
ax5 = nexttile(t_conv, 6);
hold on;
pos = 1;
for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        group_name = sprintf('%s_%s', brain_regions{br_idx}, cell_types{ct_idx});
        if isfield(integration_stats, group_name)
            int_idx = integration_stats.(group_name).integration_index;

            x_jitter = pos + 0.2*(rand(size(int_idx))-0.5);
            scatter(x_jitter, int_idx, 40, 'filled', 'MarkerFaceAlpha', 0.5);

            % Mean line
            plot([pos-0.3 pos+0.3], [mean(int_idx) mean(int_idx)], 'k-', 'LineWidth', 2);

            pos = pos + 1;
        end
    end
end
plot([0 pos], [0 0], 'k--', 'LineWidth', 1.5);
plot([0 pos], [0.2 0.2], 'r--', 'LineWidth', 1);
plot([0 pos], [-0.2 -0.2], 'b--', 'LineWidth', 1);
hold off;
xlim([0.5 pos-0.5]);
xticks(1:length(group_names_int));
xticklabels(group_names_int);
xtickangle(45);
ylabel('Integration Index', 'FontSize', 11);
title('CS+US Integration (Dual-responsive)', 'FontSize', 12);
grid on;

% 6. Selectivity index scatter
ax6 = nexttile(t_conv, 7, [1 2]);
hold on;
pos = 1;
for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        group_name = sprintf('%s_%s', brain_regions{br_idx}, cell_types{ct_idx});
        if isfield(selectivity_stats, group_name)
            sel_idx = selectivity_stats.(group_name).selectivity_index;

            x_jitter = pos + 0.2*(rand(size(sel_idx))-0.5);
            scatter(x_jitter, sel_idx, 40, 'filled', 'MarkerFaceAlpha', 0.5);

            % Mean line
            plot([pos-0.3 pos+0.3], [mean(sel_idx) mean(sel_idx)], 'k-', 'LineWidth', 2);

            pos = pos + 1;
        end
    end
end
plot([0 pos], [0 0], 'k--', 'LineWidth', 1.5);
plot([0 pos], [0.2 0.2], 'r--', 'LineWidth', 1);
plot([0 pos], [-0.2 -0.2], 'b--', 'LineWidth', 1);
hold off;
xlim([0.5 pos-0.5]);
xticks(1:length(group_names_sel));
xticklabels(group_names_sel);
xtickangle(45);
ylabel('Selectivity Index', 'FontSize', 11);
title('CS vs US Selectivity (Dual-responsive)', 'FontSize', 12);
grid on;

t_conv.Title.String = 'Convergence, Integration, and Selectivity Analysis';
t_conv.Title.FontSize = 15;

fprintf('Convergence and integration figure created.\n');

%% 8. CELL TYPE AND REGION COMPARISONS

fprintf('\n--- 8. CELL TYPE AND REGION COMPARISONS ---\n');

% Compare convergence rates
fprintf('\nConvergence Rate Comparisons:\n');

% LA vs BA
la_convergence = [];
ba_convergence = [];
for ct_idx = 1:length(cell_types)
    la_name = sprintf('LA_%s', cell_types{ct_idx});
    ba_name = sprintf('BA_%s', cell_types{ct_idx});

    if isfield(convergence_stats, la_name) && ~isempty(convergence_stats.(la_name).total)
        la_conv = length(convergence_stats.(la_name).both) / length(convergence_stats.(la_name).total);
        la_convergence = [la_convergence; la_conv];
    end

    if isfield(convergence_stats, ba_name) && ~isempty(convergence_stats.(ba_name).total)
        ba_conv = length(convergence_stats.(ba_name).both) / length(convergence_stats.(ba_name).total);
        ba_convergence = [ba_convergence; ba_conv];
    end
end

fprintf('  LA convergence rate: %.3f ± %.3f (n=%d groups)\n', mean(la_convergence), std(la_convergence), length(la_convergence));
fprintf('  BA convergence rate: %.3f ± %.3f (n=%d groups)\n', mean(ba_convergence), std(ba_convergence), length(ba_convergence));

% PN vs IN
pn_convergence = [];
in_convergence = [];
for br_idx = 1:length(brain_regions)
    pn_name = sprintf('%s_PN', brain_regions{br_idx});
    in_name = sprintf('%s_IN', brain_regions{br_idx});

    if isfield(convergence_stats, pn_name) && ~isempty(convergence_stats.(pn_name).total)
        pn_conv = length(convergence_stats.(pn_name).both) / length(convergence_stats.(pn_name).total);
        pn_convergence = [pn_convergence; pn_conv];
    end

    if isfield(convergence_stats, in_name) && ~isempty(convergence_stats.(in_name).total)
        in_conv = length(convergence_stats.(in_name).both) / length(convergence_stats.(in_name).total);
        in_convergence = [in_convergence; in_conv];
    end
end

fprintf('  PN convergence rate: %.3f ± %.3f (n=%d groups)\n', mean(pn_convergence), std(pn_convergence), length(pn_convergence));
fprintf('  IN convergence rate: %.3f ± %.3f (n=%d groups)\n', mean(in_convergence), std(in_convergence), length(in_convergence));

% Compare latencies (CS responses)
fprintf('\nLatency Comparisons (CS responses):\n');

if isfield(temporal_stats, 'LA_PN') && isfield(temporal_stats.LA_PN, ttl_fn_cs)
    fprintf('  LA PN: %.2f ± %.2f ms (n=%d)\n', 1000*temporal_stats.LA_PN.(ttl_fn_cs).mean_latency, ...
           1000*temporal_stats.LA_PN.(ttl_fn_cs).std_latency, temporal_stats.LA_PN.(ttl_fn_cs).n);
end

if isfield(temporal_stats, 'LA_IN') && isfield(temporal_stats.LA_IN, ttl_fn_cs)
    fprintf('  LA IN: %.2f ± %.2f ms (n=%d)\n', 1000*temporal_stats.LA_IN.(ttl_fn_cs).mean_latency, ...
           1000*temporal_stats.LA_IN.(ttl_fn_cs).std_latency, temporal_stats.LA_IN.(ttl_fn_cs).n);
end

if isfield(temporal_stats, 'BA_PN') && isfield(temporal_stats.BA_PN, ttl_fn_cs)
    fprintf('  BA PN: %.2f ± %.2f ms (n=%d)\n', 1000*temporal_stats.BA_PN.(ttl_fn_cs).mean_latency, ...
           1000*temporal_stats.BA_PN.(ttl_fn_cs).std_latency, temporal_stats.BA_PN.(ttl_fn_cs).n);
end

if isfield(temporal_stats, 'BA_IN') && isfield(temporal_stats.BA_IN, ttl_fn_cs)
    fprintf('  BA IN: %.2f ± %.2f ms (n=%d)\n', 1000*temporal_stats.BA_IN.(ttl_fn_cs).mean_latency, ...
           1000*temporal_stats.BA_IN.(ttl_fn_cs).std_latency, temporal_stats.BA_IN.(ttl_fn_cs).n);
end

% Compare integration indices
fprintf('\nIntegration Index Comparisons:\n');

all_integration = [];
all_groups = {};
for br_idx = 1:length(brain_regions)
    for ct_idx = 1:length(cell_types)
        group_name = sprintf('%s_%s', brain_regions{br_idx}, cell_types{ct_idx});
        if isfield(integration_stats, group_name)
            all_integration{end+1} = integration_stats.(group_name).integration_index;
            all_groups{end+1} = group_name;
            fprintf('  %s: %.3f ± %.3f\n', strrep(group_name, '_', ' '), ...
                   mean(integration_stats.(group_name).integration_index), ...
                   std(integration_stats.(group_name).integration_index));
        end
    end
end

fprintf('\n========================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('========================================\n');

% Helper function for significance stars
function sig_str = get_sig_str(p)
    if p < 0.001
        sig_str = '***';
    elseif p < 0.01
        sig_str = '**';
    elseif p < 0.05
        sig_str = '*';
    else
        sig_str = 'n.s.';
    end
end
