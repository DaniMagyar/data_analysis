cell_metrics = BAfc_load_neurons('NP_BAfc_triptest');
% Find excited neurons
[spx_sounds, idx_sounds, zResp_sound] =   BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_sound_only',    'exc', 5, 5, [0.051 0.5], 0.001, 'artefactLength',     0, 'psth_out', [0.051 0.5]);
[spx_shocks, idx_shocks, zResp_shock] =   BAfc_find_response('NP_BAfc_triptest', 'shock','TTL_triptest_shocks_only',       'exc', 5, 5, [0.051 0.5], 0.001, 'artefactLength', 0.012, 'psth_out', [0.051 0.5]);
[spx_both, idx_both, zResp_both] =        BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_both',          'exc', 5, 5, [0.051 0.5], 0.001, 'artefactLength', 0.012, 'psth_out', [0.051 0.5]);

% Find location of neurons
idx_LA = find(contains(cell_metrics.brainRegion, 'LA'));
idx_BA = find(contains(cell_metrics.brainRegion, 'BA'));



%% Zscore plot
% Combine the data
nexttile([2,1])
data_zResp = {zResp_sound(intersect(idx_sounds,idx_LA)),...
        zResp_sound(intersect(idx_sounds,idx_BA)),...
        zResp_shock(intersect(idx_shocks,idx_LA)),...
        zResp_shock(intersect(idx_shocks,idx_BA)),...
        zResp_both(intersect(idx_both,idx_LA)),...      
        zResp_both(intersect(idx_both,idx_BA))};
% Determine the maximum length among the data arrays
maxLength_zresp = max(cellfun(@length, data_zResp));
% Pad shorter arrays with NaNs
paddedData_zResp = nan(maxLength_zresp, length(data_zResp)); % Initialize with NaNs
for i = 1:length(data_zResp)
    len = length(data_zResp{i});
    paddedData_zResp(1:len, i) = data_zResp{i}; % Fill with actual data
end
h2 = boxplot(paddedData_zResp, 'Colors',[0.8500 0.3250 0.0980; 0 0.4470 0.7410;0.8500 0.3250 0.0980; 0 0.4470 0.7410;0.8500 0.3250 0.0980; 0 0.4470 0.7410], 'Symbol', '', 'Widths', 0.6);
set(findobj(h2, 'type', 'line'), 'LineWidth', 2); % Adjust the value (e.g., 2) as desired
hold on;
hLegend = [patch(nan, nan, [0.8500 0.3250 0.0980]), patch(nan, nan, [0 0.4470 0.7410])];
legend(hLegend, 'LA', 'BA', 'Location', 'northeast');
set(gca, 'xtick', [1.5, 3.5, 5.5], 'xticklabel', {'Tone', 'Shock', 'Tone + Shock'});
set(h2,{'linew'},{2})
%ylim([0 100])
ylabel('Z-score change (median)');
title('Response distribution');
hold off;

%% Spike number plot
nexttile([2,1])
allspx_sounds = sum(spx_sounds,2);
allspx_shocks = sum(spx_shocks,2);
allspx_both = sum(spx_both,2);
data_spx = {allspx_sounds(intersect(idx_sounds,idx_LA)),...
        allspx_sounds(intersect(idx_sounds,idx_BA)),...
        allspx_shocks(intersect(idx_shocks,idx_LA)),...
        allspx_shocks(intersect(idx_shocks,idx_BA)),...
        allspx_both(intersect(idx_both,idx_LA)),...      
        allspx_both(intersect(idx_both,idx_BA))};
% Determine the maximum length among the data arrays
maxLength_spx = max(cellfun(@length, data_spx));
% Pad shorter arrays with NaNs
paddedData_spx = nan(maxLength_spx, length(data_spx)); % Initialize with NaNs
for i = 1:length(data_spx)
    len = length(data_spx{i});
    paddedData_spx(1:len, i) = data_spx{i}; % Fill with actual data
end
h3 = boxplot(paddedData_spx, 'Colors',[0.8500 0.3250 0.0980; 0 0.4470 0.7410;0.8500 0.3250 0.0980; 0 0.4470 0.7410;0.8500 0.3250 0.0980; 0 0.4470 0.7410], 'Symbol', '', 'Widths', 0.6);
set(findobj(h3, 'type', 'line'), 'LineWidth', 2); % Adjust the value (e.g., 2) as desired
hold on;
hLegend = [patch(nan, nan, [0.8500 0.3250 0.0980]), patch(nan, nan, [0 0.4470 0.7410])];
legend(hLegend, 'LA', 'BA', 'Location', 'northeast');
set(gca, 'xtick', [1.5, 3.5, 5.5], 'xticklabel', {'Tone', 'Shock', 'Tone + Shock'});
set(h3,{'linew'},{2})
%ylim([0 300])
ylabel('Number of spikes');
title('Spikes in 500 ms');
hold off;

cell_metrics.labels(1:numel(cell_metrics.cellID)) = {'0'};
%cell_metrics.labels(intersect(idx_sounds_exc,idx_LA)) = {'LA_sounds_exc'};
cell_metrics.labels(intersect(idx_shocks,idx_LA)) = {'LA_shocks_exc'};
%cell_metrics.labels(intersect(idx_both_exc,idx_LA)) = {'LA_both_exc'};
%cell_metrics.labels(intersect(idx_sounds_exc,idx_BA)) = {'BA_sound_exc'};
cell_metrics.labels(intersect(idx_shocks,idx_BA)) = {'BA_shocks_exc'};
%cell_metrics.labels(intersect(idx_both_exc,idx_BA)) = {'BA_both_exc'};




CellExplorer('metrics',cell_metrics);
