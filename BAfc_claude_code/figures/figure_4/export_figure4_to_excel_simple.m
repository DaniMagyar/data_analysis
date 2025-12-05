function export_figure4_to_excel_simple(results_all, kw_data_storage, kw_results, contingency_table, ...
    chi2_obs, p_perm, cramers_v, brain_regions, cluster_names, cell_metrics, g, output_filename)

    if exist(output_filename, 'file')
        delete(output_filename);
    end

    stim_names = {'CS', 'US', 'CS+US'};

    %% PANELS B & D: ΔFR Bar Charts (main figure)
    sheet_data = {};
    sheet_data{1, 1} = 'PANELS B & D: ΔFR Bar Charts (Monosynaptic 12-25ms window)';
    sheet_data{2, 1} = 'ΔFR calculated as: (nSpikes/trial) / window_duration';
    sheet_data{3, 1} = '';
    row = 4;

    for br = 1:2
        region_name = get_region_name(brain_regions{br});
        sheet_data{row, 1} = sprintf('=== %s ===', region_name);
        row = row + 1;

        for c = 1:3
            if isempty(kw_data_storage{br, c, 1})
                continue;
            end

            sheet_data{row, 1} = cluster_names{c};
            sheet_data{row, 2} = sprintf('(n = %d neurons)', length(kw_data_storage{br, c, 1}));
            row = row + 1;

            % Header
            sheet_data{row, 1} = 'Stimulus';
            sheet_data{row, 2} = 'Mean (Hz)';
            sheet_data{row, 3} = 'SEM (Hz)';
            sheet_data{row, 4} = 'Median (Hz)';
            sheet_data{row, 5} = 'SD (Hz)';
            row = row + 1;

            % Data for each stimulus
            for stim = 1:3
                data_fr = kw_data_storage{br, c, stim};
                sheet_data{row, 1} = stim_names{stim};
                sheet_data{row, 2} = mean(data_fr);
                sheet_data{row, 3} = std(data_fr)/sqrt(length(data_fr));
                sheet_data{row, 4} = median(data_fr);
                sheet_data{row, 5} = std(data_fr);
                row = row + 1;
            end

            % Statistical tests
            sheet_data{row, 1} = '';
            row = row + 1;

            if ~isempty(kw_results) && numel(kw_results) >= (br + (c-1)*2)
                p_friedman = kw_results(br, c).p_friedman;
                p_values = kw_results(br, c).p_values;

                sheet_data{row, 1} = 'Friedman test p-value:';
                sheet_data{row, 2} = p_friedman;
                sheet_data{row, 3} = format_significance(p_friedman);
                row = row + 1;

                if p_friedman < 0.05
                    sheet_data{row, 1} = 'Post-hoc (Wilcoxon signed-rank):';
                    row = row + 1;

                    sheet_data{row, 1} = '  CS vs US:';
                    sheet_data{row, 2} = p_values(1);
                    sheet_data{row, 3} = format_significance(p_values(1));
                    row = row + 1;

                    sheet_data{row, 1} = '  CS vs CS+US:';
                    sheet_data{row, 2} = p_values(2);
                    sheet_data{row, 3} = format_significance(p_values(2));
                    row = row + 1;

                    sheet_data{row, 1} = '  US vs CS+US:';
                    sheet_data{row, 2} = p_values(3);
                    sheet_data{row, 3} = format_significance(p_values(3));
                    row = row + 1;
                else
                    sheet_data{row, 1} = '  (Post-hoc not performed - Friedman test n.s.)';
                    row = row + 1;
                end
            end
            row = row + 1;
        end
        row = row + 1;
    end

    writecell(sheet_data, output_filename, 'Sheet', 'PanelsBD_DeltaFR');

    %% PANEL E: Pie Charts
    sheet_data = {};
    sheet_data{1, 1} = 'PANEL E: Pie Charts - Cluster Proportions (Monosynaptic)';
    sheet_data{2, 1} = '';
    sheet_data{3, 1} = 'Region';
    sheet_data{3, 2} = 'CS-sel (n)';
    sheet_data{3, 3} = 'CS-sel (%)';
    sheet_data{3, 4} = 'US-sel (n)';
    sheet_data{3, 5} = 'US-sel (%)';
    sheet_data{3, 6} = 'Multi (n)';
    sheet_data{3, 7} = 'Multi (%)';
    sheet_data{3, 8} = 'Total';

    row = 4;
    for br = 1:2
        region_name = get_region_name(brain_regions{br});
        total = sum(contingency_table(br, :));
        sheet_data{row, 1} = region_name;
        sheet_data{row, 2} = contingency_table(br, 1);
        sheet_data{row, 3} = 100 * contingency_table(br, 1) / total;
        sheet_data{row, 4} = contingency_table(br, 2);
        sheet_data{row, 5} = 100 * contingency_table(br, 2) / total;
        sheet_data{row, 6} = contingency_table(br, 3);
        sheet_data{row, 7} = 100 * contingency_table(br, 3) / total;
        sheet_data{row, 8} = total;
        row = row + 1;
    end

    % Chi-square test results
    sheet_data{row, 1} = '';
    row = row + 1;
    sheet_data{row, 1} = '=== Chi-square Test (LA vs AStria) ===';
    row = row + 1;
    sheet_data{row, 1} = 'Chi-square statistic:';
    sheet_data{row, 2} = chi2_obs;
    row = row + 1;
    sheet_data{row, 1} = 'Permutation p-value:';
    sheet_data{row, 2} = p_perm;
    sheet_data{row, 3} = format_significance(p_perm);
    row = row + 1;
    sheet_data{row, 1} = 'Cramer''s V (effect size):';
    sheet_data{row, 2} = cramers_v;
    row = row + 1;

    writecell(sheet_data, output_filename, 'Sheet', 'PanelE_PieCharts');

    %% PANEL F: Across-region comparison bars
    sheet_data = {};
    sheet_data{1, 1} = 'PANEL F: LA vs AStria Comparison (Monosynaptic 12-25ms)';
    sheet_data{2, 1} = '';
    row = 3;

    for stim = 1:3
        sheet_data{row, 1} = sprintf('=== %s trials ===', stim_names{stim});
        row = row + 1;

        sheet_data{row, 1} = 'Region';
        sheet_data{row, 2} = 'Mean (Hz)';
        sheet_data{row, 3} = 'SEM (Hz)';
        sheet_data{row, 4} = 'Median (Hz)';
        sheet_data{row, 5} = 'n neurons';
        row = row + 1;

        % Collect data for both regions
        LA_data = [];
        Astria_data = [];
        for br = 1:2
            if isempty(results_all{br})
                continue;
            end
            delta_fr = calculate_delta_fr_for_stim(results_all{br}, stim, g);
            if ~isempty(delta_fr)
                region_name = get_region_name(brain_regions{br});
                sheet_data{row, 1} = region_name;
                sheet_data{row, 2} = mean(delta_fr);
                sheet_data{row, 3} = std(delta_fr)/sqrt(length(delta_fr));
                sheet_data{row, 4} = median(delta_fr);
                sheet_data{row, 5} = length(delta_fr);
                row = row + 1;
                if br == 1
                    LA_data = delta_fr;
                else
                    Astria_data = delta_fr;
                end
            end
        end

        % Statistical test
        sheet_data{row, 1} = '';
        row = row + 1;
        if ~isempty(LA_data) && ~isempty(Astria_data)
            [p_val, ~] = ranksum(LA_data, Astria_data);
            sheet_data{row, 1} = 'Wilcoxon rank-sum (LA vs AStria):';
            sheet_data{row, 2} = p_val;
            sheet_data{row, 3} = format_significance(p_val);
            row = row + 1;
        end
        row = row + 1;
    end

    writecell(sheet_data, output_filename, 'Sheet', 'PanelF_RegionComparison');

    %% PANEL G: Onset latency comparison
    sheet_data = {};
    sheet_data{1, 1} = 'PANEL G: LA vs AStria Latency Comparison (Monosynaptic 0-25ms)';
    sheet_data{2, 1} = '';
    row = 3;

    for stim = 1:3
        sheet_data{row, 1} = sprintf('=== %s trials ===', stim_names{stim});
        row = row + 1;

        sheet_data{row, 1} = 'Region';
        sheet_data{row, 2} = 'Mean (ms)';
        sheet_data{row, 3} = 'SEM (ms)';
        sheet_data{row, 4} = 'Median (ms)';
        sheet_data{row, 5} = 'SD (ms)';
        sheet_data{row, 6} = 'Range (ms)';
        sheet_data{row, 7} = 'n neurons';
        row = row + 1;

        % Collect data for both regions
        LA_data = [];
        Astria_data = [];
        for br = 1:2
            if isempty(results_all{br})
                continue;
            end
            onset_lat = extract_onset_latency(results_all{br}, stim);
            if ~isempty(onset_lat)
                region_name = get_region_name(brain_regions{br});
                sheet_data{row, 1} = region_name;
                sheet_data{row, 2} = mean(onset_lat);
                sheet_data{row, 3} = std(onset_lat)/sqrt(length(onset_lat));
                sheet_data{row, 4} = median(onset_lat);
                sheet_data{row, 5} = std(onset_lat);
                sheet_data{row, 6} = sprintf('[%.2f, %.2f]', min(onset_lat), max(onset_lat));
                sheet_data{row, 7} = length(onset_lat);
                row = row + 1;
                if br == 1
                    LA_data = onset_lat;
                else
                    Astria_data = onset_lat;
                end
            end
        end

        % Statistical test
        sheet_data{row, 1} = '';
        row = row + 1;
        if ~isempty(LA_data) && ~isempty(Astria_data)
            [p_val, ~] = ranksum(LA_data, Astria_data);
            sheet_data{row, 1} = 'Wilcoxon rank-sum (LA vs AStria):';
            sheet_data{row, 2} = p_val;
            sheet_data{row, 3} = format_significance(p_val);
            row = row + 1;
        else
            sheet_data{row, 1} = 'Insufficient data for comparison';
            row = row + 1;
        end
        row = row + 1;
    end

    writecell(sheet_data, output_filename, 'Sheet', 'PanelG_LatencyComparison');

    %% RAW DATA: ΔFR and Onset Latency values
    sheet_data = {};
    sheet_data{1, 1} = 'RAW DATA: Individual Neuron ΔFR and Onset Latency Values';
    sheet_data{2, 1} = 'ΔFR calculated as: (nSpikes/trial) / window_duration';
    sheet_data{3, 1} = 'These are the individual data points used to calculate means and SEMs in Panels B, D & G';
    sheet_data{4, 1} = '';
    row = 5;

    for br = 1:2
        region_name = get_region_name(brain_regions{br});
        sheet_data{row, 1} = sprintf('========== %s ==========', region_name);
        row = row + 1;
        row = row + 1;

        % Get global indices for this region (from cell_metrics)
        idx_neurons = strcmp(cell_metrics.brainRegion, brain_regions{br});
        global_indices = find(idx_neurons);
        res = results_all{br};

        for c = 1:3
            if isempty(kw_data_storage{br, c, 1})
                continue;
            end

            sheet_data{row, 1} = sprintf('%s (n = %d neurons)', cluster_names{c}, length(kw_data_storage{br, c, 1}));
            row = row + 1;

            % Get neurons in this cluster (these are indices into results_all{br})
            clust_idx = find(res.Clusters_all == c);
            n_neurons = length(clust_idx);

            % Calculate onset latencies for each neuron
            CS_onset_lat = nan(n_neurons, 1);
            US_onset_lat = nan(n_neurons, 1);
            Both_onset_lat = nan(n_neurons, 1);

            for n = 1:n_neurons
                local_idx = clust_idx(n);
                CS_onset_lat(n) = res.CS_onset_lat_all(local_idx) * 1000;
                US_onset_lat(n) = res.US_onset_lat_all(local_idx) * 1000;
                Both_onset_lat(n) = res.Both_onset_lat_all(local_idx) * 1000;
            end

            % Header
            sheet_data{row, 1} = 'Local #';
            sheet_data{row, 2} = 'Global Index';
            sheet_data{row, 3} = 'Animal ID';
            sheet_data{row, 4} = 'ΔFR CS (Hz)';
            sheet_data{row, 5} = 'ΔFR US (Hz)';
            sheet_data{row, 6} = 'ΔFR CS+US (Hz)';
            sheet_data{row, 7} = '';
            sheet_data{row, 8} = 'Onset Lat CS (ms)';
            sheet_data{row, 9} = 'Onset Lat US (ms)';
            sheet_data{row, 10} = 'Onset Lat CS+US (ms)';
            row = row + 1;

            for n = 1:n_neurons
                local_idx = clust_idx(n);
                global_idx = global_indices(local_idx);
                animal_id = cell_metrics.animal{global_idx};

                sheet_data{row, 1} = n;
                sheet_data{row, 2} = global_idx;
                sheet_data{row, 3} = animal_id;
                sheet_data{row, 4} = kw_data_storage{br, c, 1}(n);
                sheet_data{row, 5} = kw_data_storage{br, c, 2}(n);
                sheet_data{row, 6} = kw_data_storage{br, c, 3}(n);
                sheet_data{row, 7} = '';
                sheet_data{row, 8} = CS_onset_lat(n);
                sheet_data{row, 9} = US_onset_lat(n);
                sheet_data{row, 10} = Both_onset_lat(n);
                row = row + 1;
            end

            % Add summary stats
            sheet_data{row, 1} = '';
            row = row + 1;
            sheet_data{row, 1} = 'Mean:';
            sheet_data{row, 2} = '';
            sheet_data{row, 3} = '';
            sheet_data{row, 4} = mean(kw_data_storage{br, c, 1});
            sheet_data{row, 5} = mean(kw_data_storage{br, c, 2});
            sheet_data{row, 6} = mean(kw_data_storage{br, c, 3});
            sheet_data{row, 7} = '';
            sheet_data{row, 8} = mean(CS_onset_lat, 'omitnan');
            sheet_data{row, 9} = mean(US_onset_lat, 'omitnan');
            sheet_data{row, 10} = mean(Both_onset_lat, 'omitnan');
            row = row + 1;

            sheet_data{row, 1} = 'SEM:';
            sheet_data{row, 2} = '';
            sheet_data{row, 3} = '';
            sheet_data{row, 4} = std(kw_data_storage{br, c, 1})/sqrt(n_neurons);
            sheet_data{row, 5} = std(kw_data_storage{br, c, 2})/sqrt(n_neurons);
            sheet_data{row, 6} = std(kw_data_storage{br, c, 3})/sqrt(n_neurons);
            sheet_data{row, 7} = '';
            CS_valid = CS_onset_lat(~isnan(CS_onset_lat));
            US_valid = US_onset_lat(~isnan(US_onset_lat));
            Both_valid = Both_onset_lat(~isnan(Both_onset_lat));
            sheet_data{row, 8} = std(CS_valid)/sqrt(length(CS_valid));
            sheet_data{row, 9} = std(US_valid)/sqrt(length(US_valid));
            sheet_data{row, 10} = std(Both_valid)/sqrt(length(Both_valid));
            row = row + 1;
            row = row + 1;
        end
        row = row + 1;
    end

    writecell(sheet_data, output_filename, 'Sheet', 'RawData_DeltaFR_OnsetLat');
end

%% Helper functions
function region_name = get_region_name(region)
    if strcmp(region, 'Astria')
        region_name = 'AStria';
    else
        region_name = region;
    end
end

function delta_fr = calculate_delta_fr_for_stim(res, stim, g)
    % Select responsive neurons based on stimulus (same as Panel F in figure)
    % CS: CS-sel + Multi (clusters 1, 3)
    % US: US-sel + Multi (clusters 2, 3)
    % CS+US: All responsive (clusters 1, 2, 3)
    if stim == 1
        responsive_idx = find(ismember(res.Clusters_all, [1 3]));
        postAP_norm = res.postAP_norm_CS;
    elseif stim == 2
        responsive_idx = find(ismember(res.Clusters_all, [2 3]));
        postAP_norm = res.postAP_norm_US;
    else
        responsive_idx = find(ismember(res.Clusters_all, [1 2 3]));
        postAP_norm = res.postAP_norm_Both;
    end

    if isempty(responsive_idx)
        delta_fr = [];
        return;
    end

    % Calculate ΔFR from nSpikes per trial in monosyn window
    window_duration = g.monosyn_window - g.artifact_exclusion;
    delta_fr = zeros(length(responsive_idx), 1);

    for n = 1:length(responsive_idx)
        idx_n = responsive_idx(n);
        global_idx = res.neuron_indices_all(idx_n);

        spike_count = 0;
        trial_count = 0;
        if ~isempty(postAP_norm{global_idx})
            for trial = 1:length(postAP_norm{global_idx})
                trial_spikes = postAP_norm{global_idx}{trial};
                spike_count = spike_count + sum(trial_spikes >= g.artifact_exclusion & trial_spikes <= g.monosyn_window);
            end
            trial_count = length(postAP_norm{global_idx});
        end
        if trial_count > 0
            delta_fr(n) = (spike_count / trial_count) / window_duration;
        end
    end
end

function onset_lat = extract_onset_latency(res, stim)
    % Select appropriate latencies based on stimulus and cluster
    if stim == 1  % CS - from CS-selective and Multisensory
        responsive_idx = find(ismember(res.Clusters_all, [1 3]));
        onset_lat = res.CS_onset_lat_all(responsive_idx);
    elseif stim == 2  % US - from US-selective and Multisensory
        responsive_idx = find(ismember(res.Clusters_all, [2 3]));
        onset_lat = res.US_onset_lat_all(responsive_idx);
    else  % CS+US - from all responsive neurons
        responsive_idx = find(ismember(res.Clusters_all, [1 2 3]));
        onset_lat = res.Both_onset_lat_all(responsive_idx);
    end

    % Remove NaN values and convert to ms
    onset_lat = onset_lat(~isnan(onset_lat)) * 1000;
end

function sig_str = format_significance(p_val)
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
