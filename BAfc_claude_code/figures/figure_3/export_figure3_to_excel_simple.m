function export_figure3_to_excel_simple(results_all, kw_data_storage, kw_results, contingency_table, ...
    chi2_obs, p_perm, cramers_v, brain_regions, cluster_names, g, output_filename)

    if exist(output_filename, 'file')
        delete(output_filename);
    end

    stim_names = {'CS', 'US', 'CS+US'};

    %% PANEL C & F: Peak FR and Response Length Bar Charts (main figure)
    sheet_data = {};
    sheet_data{1, 1} = 'PANELS C & F: Peak FR and Response Length Bar Charts';
    sheet_data{2, 1} = '';
    row = 3;

    for br = 1:2
        if strcmp(brain_regions{br}, 'Astria')
            region_name = 'AStria';
        else
            region_name = brain_regions{br};
        end

        sheet_data{row, 1} = sprintf('=== %s ===', region_name);
        row = row + 1;

        for c = 1:3
            if isempty(kw_data_storage{br, c, 1})
                continue;
            end

            res = results_all{br};
            clust_idx = find(res.Clusters == c);

            sheet_data{row, 1} = cluster_names{c};
            sheet_data{row, 2} = sprintf('(n = %d neurons)', length(kw_data_storage{br, c, 1}));
            row = row + 1;

            % Calculate response length for each stimulus
            CS_metric = zeros(length(clust_idx), 1);
            US_metric = zeros(length(clust_idx), 1);
            Both_metric = zeros(length(clust_idx), 1);

            for n = 1:length(clust_idx)
                idx_n = clust_idx(n);

                % CS response length
                if ~isnan(res.CS_onset_lat(idx_n)) && ~isnan(res.CS_offset_lat(idx_n))
                    CS_metric(n) = (res.CS_offset_lat(idx_n) - res.CS_onset_lat(idx_n)) * 1000;  % Convert to ms
                else
                    CS_metric(n) = 0;  % No response detected
                end

                % US response length
                if ~isnan(res.US_onset_lat(idx_n)) && ~isnan(res.US_offset_lat(idx_n))
                    US_metric(n) = (res.US_offset_lat(idx_n) - res.US_onset_lat(idx_n)) * 1000;
                else
                    US_metric(n) = 0;
                end

                % CS+US response length
                if ~isnan(res.Both_onset_lat(idx_n)) && ~isnan(res.Both_offset_lat(idx_n))
                    Both_metric(n) = (res.Both_offset_lat(idx_n) - res.Both_onset_lat(idx_n)) * 1000;
                else
                    Both_metric(n) = 0;
                end
            end

            % Perform Friedman test for response length
            response_length_data = [CS_metric, US_metric, Both_metric];
            [p_friedman_rl, ~, stats_rl] = friedman(response_length_data, 1, 'off');

            % Get post-hoc for response length
            if p_friedman_rl < 0.05
                [p_cs_us_rl, ~] = signrank(CS_metric, US_metric);
                [p_cs_both_rl, ~] = signrank(CS_metric, Both_metric);
                [p_us_both_rl, ~] = signrank(US_metric, Both_metric);
            end

            % Header with both Delta FR and Response Length
            sheet_data{row, 1} = 'Stimulus';
            sheet_data{row, 2} = 'ΔPeak FR Mean (Hz)';
            sheet_data{row, 3} = 'ΔPeak FR SEM (Hz)';
            sheet_data{row, 4} = 'ΔPeak FR Median (Hz)';
            sheet_data{row, 5} = 'ΔPeak FR SD (Hz)';
            sheet_data{row, 6} = '';  % Empty column separator
            sheet_data{row, 7} = 'Resp Length Mean (ms)';
            sheet_data{row, 8} = 'Resp Length SEM (ms)';
            sheet_data{row, 9} = 'Resp Length Median (ms)';
            sheet_data{row, 10} = 'Resp Length SD (ms)';
            row = row + 1;

            % Data for each stimulus (side by side)
            stim_metrics = {kw_data_storage{br, c, 1}, kw_data_storage{br, c, 2}, kw_data_storage{br, c, 3}};
            resp_metrics = {CS_metric, US_metric, Both_metric};

            for stim = 1:3
                data_fr = stim_metrics{stim};
                data_rl = resp_metrics{stim};

                sheet_data{row, 1} = stim_names{stim};
                sheet_data{row, 2} = mean(data_fr);
                sheet_data{row, 3} = std(data_fr)/sqrt(length(data_fr));
                sheet_data{row, 4} = median(data_fr);
                sheet_data{row, 5} = std(data_fr);
                sheet_data{row, 6} = '';  % Empty column separator
                sheet_data{row, 7} = mean(data_rl);
                sheet_data{row, 8} = std(data_rl)/sqrt(length(data_rl));
                sheet_data{row, 9} = median(data_rl);
                sheet_data{row, 10} = std(data_rl);
                row = row + 1;
            end

            % Add statistical test results side by side
            sheet_data{row, 1} = '';
            row = row + 1;

            % Statistical test headers
            sheet_data{row, 1} = 'Statistical Tests:';
            sheet_data{row, 2} = 'ΔPeak FR';
            sheet_data{row, 3} = '';
            sheet_data{row, 4} = '';
            sheet_data{row, 5} = '';
            sheet_data{row, 6} = '';  % Empty column separator
            sheet_data{row, 7} = 'Response Length';
            row = row + 1;

            % Friedman test results
            if ~isempty(kw_results) && numel(kw_results) >= (br + (c-1)*2)
                p_friedman = kw_results(br, c).p_friedman;
                p_values = kw_results(br, c).p_values;

                sheet_data{row, 1} = 'Friedman test p-value:';
                sheet_data{row, 2} = p_friedman;
                sheet_data{row, 3} = format_significance(p_friedman);
                sheet_data{row, 4} = '';
                sheet_data{row, 5} = '';
                sheet_data{row, 6} = '';  % Empty column separator
                sheet_data{row, 7} = p_friedman_rl;
                sheet_data{row, 8} = format_significance(p_friedman_rl);
                row = row + 1;

                % Post-hoc header
                if p_friedman < 0.05 || p_friedman_rl < 0.05
                    sheet_data{row, 1} = 'Post-hoc (Wilcoxon signed-rank):';
                    row = row + 1;

                    % CS vs US
                    sheet_data{row, 1} = '  CS vs US:';
                    if p_friedman < 0.05
                        sheet_data{row, 2} = p_values(1);
                        sheet_data{row, 3} = format_significance(p_values(1));
                    else
                        sheet_data{row, 2} = 'n/a';
                        sheet_data{row, 3} = '';
                    end
                    sheet_data{row, 4} = '';
                    sheet_data{row, 5} = '';
                    sheet_data{row, 6} = '';  % Empty column separator
                    if p_friedman_rl < 0.05
                        sheet_data{row, 7} = p_cs_us_rl;
                        sheet_data{row, 8} = format_significance(p_cs_us_rl);
                    else
                        sheet_data{row, 7} = 'n/a';
                        sheet_data{row, 8} = '';
                    end
                    row = row + 1;

                    % CS vs CS+US
                    sheet_data{row, 1} = '  CS vs CS+US:';
                    if p_friedman < 0.05
                        sheet_data{row, 2} = p_values(2);
                        sheet_data{row, 3} = format_significance(p_values(2));
                    else
                        sheet_data{row, 2} = 'n/a';
                        sheet_data{row, 3} = '';
                    end
                    sheet_data{row, 4} = '';
                    sheet_data{row, 5} = '';
                    sheet_data{row, 6} = '';  % Empty column separator
                    if p_friedman_rl < 0.05
                        sheet_data{row, 7} = p_cs_both_rl;
                        sheet_data{row, 8} = format_significance(p_cs_both_rl);
                    else
                        sheet_data{row, 7} = 'n/a';
                        sheet_data{row, 8} = '';
                    end
                    row = row + 1;

                    % US vs CS+US
                    sheet_data{row, 1} = '  US vs CS+US:';
                    if p_friedman < 0.05
                        sheet_data{row, 2} = p_values(3);
                        sheet_data{row, 3} = format_significance(p_values(3));
                    else
                        sheet_data{row, 2} = 'n/a';
                        sheet_data{row, 3} = '';
                    end
                    sheet_data{row, 4} = '';
                    sheet_data{row, 5} = '';
                    sheet_data{row, 6} = '';  % Empty column separator
                    if p_friedman_rl < 0.05
                        sheet_data{row, 7} = p_us_both_rl;
                        sheet_data{row, 8} = format_significance(p_us_both_rl);
                    else
                        sheet_data{row, 7} = 'n/a';
                        sheet_data{row, 8} = '';
                    end
                    row = row + 1;
                else
                    sheet_data{row, 1} = '  (Post-hoc not performed - Both Friedman tests n.s.)';
                    row = row + 1;
                end
            end
            row = row + 1;
        end
        row = row + 1;
    end

    writecell(sheet_data, output_filename, 'Sheet', 'PanelsCF_PeakFR_RespLen');

    %% PANEL G: Pie Charts
    sheet_data = {};
    sheet_data{1, 1} = 'PANEL G: Pie Charts - Cluster Proportions';
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
        if strcmp(brain_regions{br}, 'Astria')
            region_name = 'AStria';
        else
            region_name = brain_regions{br};
        end

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

    % Add chi-square test results
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

    writecell(sheet_data, output_filename, 'Sheet', 'PanelG_PieCharts');

    %% PANEL H: Across-region comparison bars
    sheet_data = {};
    sheet_data{1, 1} = 'PANEL H: LA vs AStria Comparison';
    sheet_data{2, 1} = '';
    row = 3;

    for stim = 1:3
        sheet_data{row, 1} = sprintf('=== %s trials ===', stim_names{stim});
        row = row + 1;

        % Header
        sheet_data{row, 1} = 'Region';
        sheet_data{row, 2} = 'Mean (Hz)';
        sheet_data{row, 3} = 'SEM (Hz)';
        sheet_data{row, 4} = 'Median (Hz)';
        sheet_data{row, 5} = 'n neurons';
        row = row + 1;

        % Collect data
        for br = 1:2
            if isempty(results_all{br})
                continue;
            end

            res = results_all{br};

            % Select appropriate peak FR data based on stimulus
            cluster_mask = get_cluster_mask_for_stimulus(res, stim);
            if stim == 1
                peak_fr_data = res.CS_peak_Hz;
            elseif stim == 2
                peak_fr_data = res.US_peak_Hz;
            else
                peak_fr_data = res.Both_peak_Hz;
            end

            responsive_idx = find(cluster_mask);

            if ~isempty(responsive_idx)
                peak_fr = peak_fr_data(responsive_idx);

                if strcmp(brain_regions{br}, 'Astria')
                    region_name = 'AStria';
                else
                    region_name = brain_regions{br};
                end

                sheet_data{row, 1} = region_name;
                sheet_data{row, 2} = mean(peak_fr);
                sheet_data{row, 3} = std(peak_fr)/sqrt(length(peak_fr));
                sheet_data{row, 4} = median(peak_fr);
                sheet_data{row, 5} = length(peak_fr);
                row = row + 1;
            end
        end

        % Add statistical test (Wilcoxon rank-sum between regions)
        sheet_data{row, 1} = '';
        row = row + 1;
        sheet_data{row, 1} = 'Statistical Test:';
        row = row + 1;

        % Collect LA and AStria data for this stimulus
        LA_data = [];
        Astria_data = [];
        for br = 1:2
            if isempty(results_all{br})
                continue;
            end
            res = results_all{br};
            cluster_mask = get_cluster_mask_for_stimulus(res, stim);
            if stim == 1
                peak_fr_data = res.CS_peak_Hz;
            elseif stim == 2
                peak_fr_data = res.US_peak_Hz;
            else
                peak_fr_data = res.Both_peak_Hz;
            end
            responsive_idx = find(cluster_mask);
            if ~isempty(responsive_idx)
                peak_fr = peak_fr_data(responsive_idx);
                if br == 1
                    LA_data = peak_fr;
                else
                    Astria_data = peak_fr;
                end
            end
        end

        if ~isempty(LA_data) && ~isempty(Astria_data)
            [p_val, ~] = ranksum(LA_data, Astria_data);
            sheet_data{row, 1} = 'Wilcoxon rank-sum (LA vs AStria):';
            sheet_data{row, 2} = p_val;
            sheet_data{row, 3} = format_significance(p_val);
            row = row + 1;
        end

        row = row + 1;
    end

    writecell(sheet_data, output_filename, 'Sheet', 'PanelH_RegionComparison');

    %% RAW DATA: Individual neuron values for bar charts
    sheet_data = {};
    sheet_data{1, 1} = 'RAW DATA: Individual Neuron Peak FR and Response Length Values';
    sheet_data{2, 1} = 'These are the individual data points used to calculate means and SEMs in Panels C & F';
    sheet_data{3, 1} = '';
    row = 4;

    for br = 1:2
        if strcmp(brain_regions{br}, 'Astria')
            region_name = 'AStria';
        else
            region_name = brain_regions{br};
        end

        sheet_data{row, 1} = sprintf('========== %s ==========', region_name);
        row = row + 1;
        row = row + 1;

        % Get global indices for this region
        idx_neurons = strcmp(g.cell_metrics.brainRegion, brain_regions{br});
        global_indices = find(idx_neurons);

        res = results_all{br};

        for c = 1:3
            if isempty(kw_data_storage{br, c, 1})
                continue;
            end

            sheet_data{row, 1} = sprintf('%s (n = %d neurons)', cluster_names{c}, length(kw_data_storage{br, c, 1}));
            row = row + 1;

            % Get neurons in this cluster
            clust_idx = find(res.Clusters == c);
            n_neurons = length(clust_idx);

            % Calculate response lengths
            CS_resp_len = zeros(n_neurons, 1);
            US_resp_len = zeros(n_neurons, 1);
            Both_resp_len = zeros(n_neurons, 1);

            for n = 1:n_neurons
                local_idx = clust_idx(n);

                % CS response length
                if ~isnan(res.CS_onset_lat(local_idx)) && ~isnan(res.CS_offset_lat(local_idx))
                    CS_resp_len(n) = (res.CS_offset_lat(local_idx) - res.CS_onset_lat(local_idx)) * 1000;
                else
                    CS_resp_len(n) = 0;
                end

                % US response length
                if ~isnan(res.US_onset_lat(local_idx)) && ~isnan(res.US_offset_lat(local_idx))
                    US_resp_len(n) = (res.US_offset_lat(local_idx) - res.US_onset_lat(local_idx)) * 1000;
                else
                    US_resp_len(n) = 0;
                end

                % CS+US response length
                if ~isnan(res.Both_onset_lat(local_idx)) && ~isnan(res.Both_offset_lat(local_idx))
                    Both_resp_len(n) = (res.Both_offset_lat(local_idx) - res.Both_onset_lat(local_idx)) * 1000;
                else
                    Both_resp_len(n) = 0;
                end
            end

            % Header
            sheet_data{row, 1} = 'Local #';
            sheet_data{row, 2} = 'Global Index';
            sheet_data{row, 3} = 'Animal ID';
            sheet_data{row, 4} = 'ΔPeak FR CS (Hz)';
            sheet_data{row, 5} = 'ΔPeak FR US (Hz)';
            sheet_data{row, 6} = 'ΔPeak FR CS+US (Hz)';
            sheet_data{row, 7} = '';  % Empty separator
            sheet_data{row, 8} = 'Resp Length CS (ms)';
            sheet_data{row, 9} = 'Resp Length US (ms)';
            sheet_data{row, 10} = 'Resp Length CS+US (ms)';
            row = row + 1;

            for n = 1:n_neurons
                local_idx = clust_idx(n);
                global_idx = global_indices(local_idx);
                animal_id = g.cell_metrics.animal{global_idx};

                sheet_data{row, 1} = n;
                sheet_data{row, 2} = global_idx;
                sheet_data{row, 3} = animal_id;
                sheet_data{row, 4} = kw_data_storage{br, c, 1}(n);  % CS ΔFR
                sheet_data{row, 5} = kw_data_storage{br, c, 2}(n);  % US ΔFR
                sheet_data{row, 6} = kw_data_storage{br, c, 3}(n);  % CS+US ΔFR
                sheet_data{row, 7} = '';  % Empty separator
                sheet_data{row, 8} = CS_resp_len(n);
                sheet_data{row, 9} = US_resp_len(n);
                sheet_data{row, 10} = Both_resp_len(n);
                row = row + 1;
            end

            % Add summary stats
            sheet_data{row, 1} = '';
            row = row + 1;
            sheet_data{row, 1} = 'Mean:';
            sheet_data{row, 4} = mean(kw_data_storage{br, c, 1});
            sheet_data{row, 5} = mean(kw_data_storage{br, c, 2});
            sheet_data{row, 6} = mean(kw_data_storage{br, c, 3});
            sheet_data{row, 7} = '';  % Empty separator
            sheet_data{row, 8} = mean(CS_resp_len);
            sheet_data{row, 9} = mean(US_resp_len);
            sheet_data{row, 10} = mean(Both_resp_len);
            row = row + 1;
            sheet_data{row, 1} = 'SEM:';
            sheet_data{row, 4} = std(kw_data_storage{br, c, 1})/sqrt(n_neurons);
            sheet_data{row, 5} = std(kw_data_storage{br, c, 2})/sqrt(n_neurons);
            sheet_data{row, 6} = std(kw_data_storage{br, c, 3})/sqrt(n_neurons);
            sheet_data{row, 7} = '';  % Empty separator
            sheet_data{row, 8} = std(CS_resp_len)/sqrt(n_neurons);
            sheet_data{row, 9} = std(US_resp_len)/sqrt(n_neurons);
            sheet_data{row, 10} = std(Both_resp_len)/sqrt(n_neurons);
            row = row + 1;

            row = row + 2;
        end
        row = row + 1;
    end

    writecell(sheet_data, output_filename, 'Sheet', 'RawData_PeakFR_RespLen');
end

function cluster_mask = get_cluster_mask_for_stimulus(res, stim)
    % Return cluster mask for Panel H region comparison
    % stim: 1=CS, 2=US, 3=CS+US
    if stim == 1
        % CS: CS-sel (cluster 1) + Multi (cluster 3)
        cluster_mask = (res.Clusters == 1) | (res.Clusters == 3);
    elseif stim == 2
        % US: US-sel (cluster 2) + Multi (cluster 3)
        cluster_mask = (res.Clusters == 2) | (res.Clusters == 3);
    else
        % CS+US: CS-sel (cluster 1) + US-sel (cluster 2) + Multi (cluster 3)
        cluster_mask = (res.Clusters == 1) | (res.Clusters == 2) | (res.Clusters == 3);
    end
end

function sig_str = format_significance(p_val)
    % Format p-value as significance stars
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
