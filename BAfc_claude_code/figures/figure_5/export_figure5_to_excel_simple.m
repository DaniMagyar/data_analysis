function export_figure5_to_excel_simple(comparison_results, cell_metrics, results_all, region_names, stim_names, g)
    output_filename = 'figure_5_data.xlsx';

    if exist(output_filename, 'file')
        delete(output_filename);
    end

    %% Sheet 1: Panel E - Pie Charts (Proportions)
    sheet_data = {};
    sheet_data{1,1} = 'PANEL E: PIE CHARTS - PROPORTIONS OF ENHANCED NEURONS';
    sheet_data{2,1} = sprintf('Test window: %.0f-%.0f ms', g.artifact_end*1000, g.monosyn_window*1000);
    sheet_data{5,1} = 'Region';
    sheet_data{5,2} = 'Stimulus';
    sheet_data{5,3} = 'Total Monosynaptic';
    sheet_data{5,4} = 'Enhanced';
    sheet_data{5,5} = 'Enhanced (%)';
    sheet_data{5,6} = 'Non-Enhanced';
    sheet_data{5,7} = 'Non-Enhanced (%)';

    row_idx = 6;
    for r = 1:2
        for s = 1:2
            result_field = sprintf('%s_%s', region_names{r}, stim_names{s});
            result = comparison_results.(result_field);

            n_enhanced = result.n_increased;
            n_non_enhanced = result.n_decreased + result.n_unchanged;

            sheet_data{row_idx, 1} = region_names{r};
            sheet_data{row_idx, 2} = stim_names{s};
            sheet_data{row_idx, 3} = result.n_total;
            sheet_data{row_idx, 4} = n_enhanced;
            if result.n_total > 0
                sheet_data{row_idx, 5} = 100 * n_enhanced / result.n_total;
            else
                sheet_data{row_idx, 5} = 0;
            end
            sheet_data{row_idx, 6} = n_non_enhanced;
            if result.n_total > 0
                sheet_data{row_idx, 7} = 100 * n_non_enhanced / result.n_total;
            else
                sheet_data{row_idx, 7} = 0;
            end
            row_idx = row_idx + 1;
        end
    end

    writecell(sheet_data, output_filename, 'Sheet', 'PanelE_PieCharts');

    %% Sheet 2: Panel F - Spaghetti Plots (FR by Classification)
    sheet_data = {};
    sheet_data{1,1} = 'PANEL F: SPAGHETTI PLOTS - FR BY CLASSIFICATION';
    sheet_data{2,1} = 'FR calculated as (nSpikes/trial) / window_duration (no baseline subtraction)';
    sheet_data{3,1} = sprintf('Test window: %.0f-%.0f ms', g.artifact_end*1000, g.monosyn_window*1000);
    sheet_data{5,1} = 'Region';
    sheet_data{5,2} = 'Stimulus';
    sheet_data{5,3} = 'Classification';
    sheet_data{5,4} = 'N Neurons';
    sheet_data{5,5} = '';
    sheet_data{5,6} = 'No Light Mean (Hz)';
    sheet_data{5,7} = 'No Light SEM (Hz)';
    sheet_data{5,8} = 'No Light Median (Hz)';
    sheet_data{5,9} = 'No Light SD (Hz)';
    sheet_data{5,10} = '';
    sheet_data{5,11} = 'With Light Mean (Hz)';
    sheet_data{5,12} = 'With Light SEM (Hz)';
    sheet_data{5,13} = 'With Light Median (Hz)';
    sheet_data{5,14} = 'With Light SD (Hz)';
    sheet_data{5,15} = '';
    sheet_data{5,16} = 'Paired t-test p';
    sheet_data{5,17} = 'Significance';

    row_idx = 6;
    for r = 1:2
        for s = 1:2
            result_field = sprintf('%s_%s', region_names{r}, stim_names{s});
            result = comparison_results.(result_field);

            res = results_all{r, s};
            postAP_norm_nolight = res.postAP_norm_nolight;
            postAP_norm_light = res.postAP_norm_light;

            % Process all three categories
            categories = {'increased', 'decreased', 'unchanged'};
            category_names = {'Increased', 'Decreased', 'Unchanged'};

            for cat = 1:3
                idx_list = result.([categories{cat} '_idx']);
                n_neurons = length(idx_list);

                sheet_data{row_idx, 1} = region_names{r};
                sheet_data{row_idx, 2} = stim_names{s};
                sheet_data{row_idx, 3} = category_names{cat};
                sheet_data{row_idx, 4} = n_neurons;

                if n_neurons > 0
                    [fr_nolight, fr_light] = calculate_fr(idx_list, postAP_norm_nolight, postAP_norm_light, g);

                    sheet_data{row_idx, 6} = mean(fr_nolight);
                    sheet_data{row_idx, 7} = std(fr_nolight) / sqrt(n_neurons);
                    sheet_data{row_idx, 8} = median(fr_nolight);
                    sheet_data{row_idx, 9} = std(fr_nolight);
                    sheet_data{row_idx, 11} = mean(fr_light);
                    sheet_data{row_idx, 12} = std(fr_light) / sqrt(n_neurons);
                    sheet_data{row_idx, 13} = median(fr_light);
                    sheet_data{row_idx, 14} = std(fr_light);

                    if n_neurons > 1
                        [~, p_pop, ~, ~] = ttest(fr_nolight, fr_light);
                        sheet_data{row_idx, 16} = p_pop;
                        sheet_data{row_idx, 17} = format_significance(p_pop);
                    else
                        sheet_data{row_idx, 16} = 'n/a';
                        sheet_data{row_idx, 17} = 'n/a';
                    end
                end

                row_idx = row_idx + 1;
            end

            % Add blank row between region/stimulus combinations
            row_idx = row_idx + 1;
        end
    end

    writecell(sheet_data, output_filename, 'Sheet', 'PanelF_SpaghetttiPlots');

    %% Sheet 3: Raw Data - Individual Neurons
    sheet_data = {};
    sheet_data{1,1} = 'RAW DATA: INDIVIDUAL NEURON CLASSIFICATIONS';
    sheet_data{2,1} = sprintf('Test window: %.0f-%.0f ms', g.artifact_end*1000, g.monosyn_window*1000);
    sheet_data{4,1} = 'Global Index';
    sheet_data{4,2} = 'Cell ID';
    sheet_data{4,3} = 'Animal ID';
    sheet_data{4,4} = 'Region';
    sheet_data{4,5} = 'Cell Type';
    sheet_data{4,6} = 'Stimulus';
    sheet_data{4,7} = 'Classification';
    sheet_data{4,8} = 'p-value';
    sheet_data{4,9} = '';
    sheet_data{4,10} = 'FR No Light (Hz)';
    sheet_data{4,11} = 'FR Light (Hz)';
    sheet_data{4,12} = 'FR Diff (Hz)';

    row_idx = 5;
    for r = 1:2
        for s = 1:2
            result_field = sprintf('%s_%s', region_names{r}, stim_names{s});
            result = comparison_results.(result_field);

            if result.n_total == 0
                continue;
            end

            res = results_all{r, s};
            postAP_norm_nolight = res.postAP_norm_nolight;
            postAP_norm_light = res.postAP_norm_light;

            % Process all three categories
            categories = {'increased_idx', 'decreased_idx', 'unchanged_idx'};
            category_names = {'Increased', 'Decreased', 'Unchanged'};

            for cat = 1:3
                idx_list = result.(categories{cat});

                for i = 1:length(idx_list)
                    idx = idx_list(i);
                    neuron_pos = find(result.all_neuron_indices == idx);

                    sheet_data{row_idx, 1} = idx;
                    sheet_data{row_idx, 2} = cell_metrics.cellID(idx);
                    sheet_data{row_idx, 3} = cell_metrics.animal{idx};
                    sheet_data{row_idx, 4} = region_names{r};
                    sheet_data{row_idx, 5} = cell_metrics.putativeCellType{idx};
                    sheet_data{row_idx, 6} = stim_names{s};
                    sheet_data{row_idx, 7} = category_names{cat};

                    if ~isempty(neuron_pos)
                        sheet_data{row_idx, 8} = result.neuron_pvalues{neuron_pos};

                        test_window = [g.artifact_end, g.monosyn_window];
                        window_duration = test_window(2) - test_window(1);
                        [fr_nolight, fr_light] = calculate_fr_single(idx, postAP_norm_nolight, postAP_norm_light, test_window, window_duration);

                        sheet_data{row_idx, 10} = fr_nolight;
                        sheet_data{row_idx, 11} = fr_light;
                        sheet_data{row_idx, 12} = fr_light - fr_nolight;
                    end
                    row_idx = row_idx + 1;
                end
            end
        end
    end

    writecell(sheet_data, output_filename, 'Sheet', 'RawData_Classifications');

    %% Sheet 4: Chi-square tests - LA vs AStria comparisons
    sheet_data = {};
    sheet_data{1,1} = 'CHI-SQUARE TESTS: LA VS ASTRIA';
    sheet_data{2,1} = 'Comparison of classification distributions between brain regions';

    LA_CS = comparison_results.LA_CS;
    LA_US = comparison_results.LA_US;
    Astria_CS = comparison_results.Astria_CS;
    Astria_US = comparison_results.Astria_US;

    % Test 1: LA CS vs AStria CS
    sheet_data{4,1} = 'TEST 1: LA CS vs AStria CS';
    sheet_data{5,1} = '';
    sheet_data{6,1} = 'Contingency Table:';
    sheet_data{7,1} = '';
    sheet_data{7,2} = 'Increased';
    sheet_data{7,3} = 'Decreased';
    sheet_data{7,4} = 'Unchanged';
    sheet_data{7,5} = 'Total';

    sheet_data{8,1} = 'LA CS';
    sheet_data{8,2} = LA_CS.n_increased;
    sheet_data{8,3} = LA_CS.n_decreased;
    sheet_data{8,4} = LA_CS.n_unchanged;
    sheet_data{8,5} = LA_CS.n_total;

    sheet_data{9,1} = 'AStria CS';
    sheet_data{9,2} = Astria_CS.n_increased;
    sheet_data{9,3} = Astria_CS.n_decreased;
    sheet_data{9,4} = Astria_CS.n_unchanged;
    sheet_data{9,5} = Astria_CS.n_total;

    contingency_table_CS = [LA_CS.n_increased, LA_CS.n_decreased, LA_CS.n_unchanged;
                            Astria_CS.n_increased, Astria_CS.n_decreased, Astria_CS.n_unchanged];
    [chi2_CS, p_CS] = calculate_chi_square(contingency_table_CS);
    df_CS = (size(contingency_table_CS,1)-1) * (size(contingency_table_CS,2)-1);

    sheet_data{11,1} = 'Chi-square statistic:';
    sheet_data{11,2} = chi2_CS;
    sheet_data{12,1} = 'Degrees of freedom:';
    sheet_data{12,2} = df_CS;
    sheet_data{13,1} = 'p-value:';
    sheet_data{13,2} = p_CS;
    sheet_data{13,3} = format_significance(p_CS);

    % Test 2: LA US vs AStria US
    sheet_data{16,1} = 'TEST 2: LA US vs AStria US';
    sheet_data{17,1} = '';
    sheet_data{18,1} = 'Contingency Table:';
    sheet_data{19,1} = '';
    sheet_data{19,2} = 'Increased';
    sheet_data{19,3} = 'Decreased';
    sheet_data{19,4} = 'Unchanged';
    sheet_data{19,5} = 'Total';

    sheet_data{20,1} = 'LA US';
    sheet_data{20,2} = LA_US.n_increased;
    sheet_data{20,3} = LA_US.n_decreased;
    sheet_data{20,4} = LA_US.n_unchanged;
    sheet_data{20,5} = LA_US.n_total;

    sheet_data{21,1} = 'AStria US';
    sheet_data{21,2} = Astria_US.n_increased;
    sheet_data{21,3} = Astria_US.n_decreased;
    sheet_data{21,4} = Astria_US.n_unchanged;
    sheet_data{21,5} = Astria_US.n_total;

    contingency_table_US = [LA_US.n_increased, LA_US.n_decreased, LA_US.n_unchanged;
                            Astria_US.n_increased, Astria_US.n_decreased, Astria_US.n_unchanged];
    [chi2_US, p_US] = calculate_chi_square(contingency_table_US);
    df_US = (size(contingency_table_US,1)-1) * (size(contingency_table_US,2)-1);

    sheet_data{23,1} = 'Chi-square statistic:';
    sheet_data{23,2} = chi2_US;
    sheet_data{24,1} = 'Degrees of freedom:';
    sheet_data{24,2} = df_US;
    sheet_data{25,1} = 'p-value:';
    sheet_data{25,2} = p_US;
    sheet_data{25,3} = format_significance(p_US);

    writecell(sheet_data, output_filename, 'Sheet', 'ChiSquare_RegionComparisons');
end

%% Helper functions
function [fr_nolight, fr_light] = calculate_fr(neuron_idx_list, postAP_norm_nolight, postAP_norm_light, g)
    % Calculate FR for multiple neurons
    n_neurons = length(neuron_idx_list);
    fr_nolight = zeros(n_neurons, 1);
    fr_light = zeros(n_neurons, 1);

    test_window = [g.artifact_end, g.monosyn_window];
    window_duration = test_window(2) - test_window(1);

    for i = 1:n_neurons
        idx = neuron_idx_list(i);
        [fr_nolight(i), fr_light(i)] = calculate_fr_single(idx, postAP_norm_nolight, postAP_norm_light, test_window, window_duration);
    end
end

function [fr_nolight, fr_light] = calculate_fr_single(idx, postAP_norm_nolight, postAP_norm_light, test_window, window_duration)
    spike_count_nolight = 0;
    trial_count_nolight = 0;
    if ~isempty(postAP_norm_nolight{idx})
        for trial = 1:length(postAP_norm_nolight{idx})
            trial_spikes = postAP_norm_nolight{idx}{trial};
            spike_count_nolight = spike_count_nolight + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
        end
        trial_count_nolight = length(postAP_norm_nolight{idx});
    end
    fr_nolight = 0;
    if trial_count_nolight > 0
        fr_nolight = (spike_count_nolight / trial_count_nolight) / window_duration;
    end

    spike_count_light = 0;
    trial_count_light = 0;
    if ~isempty(postAP_norm_light{idx})
        for trial = 1:length(postAP_norm_light{idx})
            trial_spikes = postAP_norm_light{idx}{trial};
            spike_count_light = spike_count_light + sum(trial_spikes >= test_window(1) & trial_spikes <= test_window(2));
        end
        trial_count_light = length(postAP_norm_light{idx});
    end
    fr_light = 0;
    if trial_count_light > 0
        fr_light = (spike_count_light / trial_count_light) / window_duration;
    end
end

function [chi2, p] = calculate_chi_square(contingency_table)
    row_totals = sum(contingency_table, 2);
    col_totals = sum(contingency_table, 1);
    grand_total = sum(row_totals);

    expected = (row_totals * col_totals) / grand_total;
    chi2 = sum(sum((contingency_table - expected).^2 ./ expected));
    df = (size(contingency_table,1)-1) * (size(contingency_table,2)-1);
    p = 1 - chi2cdf(chi2, df);
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
