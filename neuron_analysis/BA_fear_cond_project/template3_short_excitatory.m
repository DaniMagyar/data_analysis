cell_metrics = BAfc_load_neurons('NP_BAfc_triptest');
% Find excited neurons
[spx_sounds, idx_sounds, zResp_sounds] =   BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_sound_only',    'exc', 0.5, 0.5, [0.013 0.05], 0.001, 'artefactLength',     0, 'psth_out', [0.013 0.05]);
[spx_shocks, idx_shocks, zResp_shocks] =   BAfc_find_response('NP_BAfc_triptest', 'shock','TTL_triptest_shocks_only',       'exc', 0.5, 0.5, [0.013 0.05], 0.001, 'artefactLength', 0.012, 'psth_out', [0.013 0.05]);
[spx_both, idx_both, zResp_both] =        BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_both',          'exc', 0.5, 0.5, [0.013 0.05], 0.001, 'artefactLength', 0.012, 'psth_out', [0.013 0.05]);

% Calculate excitatory response latencies (ms)
[lat_sounds] =  BAfc_calc_latency('NP_BAfc_triptest', 'TTL_triptest_sound_only',     0.5, 0.05, 0.001);
[lat_shocks] =  BAfc_calc_latency('NP_BAfc_triptest', 'TTL_triptest_shocks_only',    0.5, 0.05, 0.001);
[lat_both] =    BAfc_calc_latency('NP_BAfc_triptest', 'TTL_triptest_both',           0.5, 0.05, 0.001);
% Find location of neurons
idx_LA = find(contains(cell_metrics.brainRegion, 'LA'));
idx_BA = find(contains(cell_metrics.brainRegion, 'BA'));

plot_50 = tiledlayout(2,7);
nexttile([2,1])
%% Latency plot
% Combine the data
data_lat = {lat_sounds(intersect(idx_sounds,idx_LA)),...
        lat_sounds(intersect(idx_sounds,idx_BA)),...
        lat_shocks(intersect(idx_shocks,idx_LA)),...
        lat_shocks(intersect(idx_shocks,idx_BA)),...
        lat_both(intersect(idx_both,idx_LA)),...      
        lat_both(intersect(idx_both,idx_BA))};
% Determine the maximum length among the data arrays
maxLength_lat = max(cellfun(@length, data_lat));
% Pad shorter arrays with NaNs
paddedData_lat = nan(maxLength_lat, length(data_lat)); % Initialize with NaNs
for i = 1:length(data_lat)
    len = length(data_lat{i});
    paddedData_lat(1:len, i) = data_lat{i}; % Fill with actual data
end
% Create boxplot
paddedData_lat(find(paddedData_lat<5)) = NaN;
h = boxplot(paddedData_lat, 'Colors',[0.8500 0.3250 0.0980; 0 0.4470 0.7410;0.8500 0.3250 0.0980; 0 0.4470 0.7410;0.8500 0.3250 0.0980; 0 0.4470 0.7410], 'Symbol', '', 'Widths', 0.6);
set(findobj(h, 'type', 'line'), 'LineWidth', 2); % Adjust the value (e.g., 2) as desired
hold on;
hLegend = [patch(nan, nan, [0.8500 0.3250 0.0980]), patch(nan, nan, [0 0.4470 0.7410])];
legend(hLegend, 'LA', 'BA', 'Location', 'northeast');
set(gca, 'xtick', [1.5, 3.5, 5.5], 'xticklabel', {'Tone', 'Shock', 'Tone + Shock'});
set(h,{'linew'},{2})
ylim([0 40])
ylabel('Latency (ms)');
title('Latency Distributions');
hold off;


%% Zscore plot
% Combine the data
nexttile([2,1])
data_zResp = {zResp_sounds(intersect(idx_sounds,idx_LA)),...
        zResp_sounds(intersect(idx_sounds,idx_BA)),...
        zResp_shocks(intersect(idx_shocks,idx_LA)),...
        zResp_shocks(intersect(idx_shocks,idx_BA)),...
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
title('Spikes in 50 ms');
hold off;

%% Zscores with lines 
nexttile
hold on;
idx_all = intersect(unique([idx_sounds; idx_shocks; idx_both]), idx_LA);
for i = idx_all
    % Plot individual data points
    plot([1, 2, 3], [zResp_sounds(i), zResp_shocks(i), zResp_both(i)], 'k-'); % Line connecting paired points
    plot(1, zResp_sounds(i), 'ro'); % Data point from vector1
    plot(2, zResp_shocks(i), 'bo'); % Data point from vector2
    plot(3, zResp_both(i), 'go'); % Data point from vector2
end
hold off;
ylim([0 100])
title('Lineplots of Z-score change LA');

nexttile
hold on;
idx_all = intersect(unique([idx_sounds; idx_shocks; idx_both]), idx_BA);
for i = idx_all
    % Plot individual data points
    plot([1, 2, 3], [zResp_sounds(i), zResp_shocks(i), zResp_both(i)], 'k-'); % Line connecting paired points
    plot(1, zResp_sounds(i), 'ro'); % Data point from vector1
    plot(2, zResp_shocks(i), 'bo'); % Data point from vector2
    plot(3, zResp_both(i), 'go'); % Data point from vector2
end
hold off;
ylim([0 100])
title('Lineplots of Z-score change in BA');


%% Allspikes with lines 
nexttile
hold on;
idx_all = intersect(unique([idx_sounds; idx_shocks; idx_both]), idx_LA);
for i = idx_all
    % Plot individual data points
    plot([1, 2, 3], [allspx_sounds(i), allspx_shocks(i), allspx_both(i)], 'k-'); % Line connecting paired points
    plot(1, allspx_sounds(i), 'ro'); % Data point from vector1
    plot(2, allspx_shocks(i), 'bo'); % Data point from vector2
    plot(3, allspx_both(i), 'go'); % Data point from vector2
end
hold off;
ylim([0 300])
title('Lineplots of spike number change LA');

nexttile
hold on;
idx_all = intersect(unique([idx_sounds; idx_shocks; idx_both]), idx_BA);
for i = idx_all
    % Plot individual data points
    plot([1, 2, 3], [allspx_sounds(i), allspx_shocks(i), allspx_both(i)], 'k-'); % Line connecting paired points
    plot(1, allspx_sounds(i), 'ro'); % Data point from vector1
    plot(2, allspx_shocks(i), 'bo'); % Data point from vector2
    plot(3, allspx_both(i), 'go'); % Data point from vector2
end
hold off;
ylim([0 300])
title('Lineplots of spike number change in BA');


%% Piechart 
idx_only_LA_sounds = intersect(setdiff(idx_sounds, idx_shocks), idx_LA);
idx_only_LA_shocks = intersect(setdiff(idx_shocks, idx_sounds), idx_LA);
idx_only_LA_both =   intersect(setdiff(idx_both, [idx_sounds; idx_shocks]), idx_LA);
idx_all_LA = intersect(intersect(idx_sounds, intersect(idx_shocks, idx_both)), idx_LA);

idx_only_BA_sounds = intersect(setdiff(idx_sounds, idx_shocks), idx_BA);
idx_only_BA_shocks = intersect(setdiff(idx_shocks, idx_sounds), idx_BA);
idx_only_BA_both =   intersect(setdiff(idx_both, [idx_sounds; idx_shocks]), idx_BA);
idx_all_BA = intersect(intersect(idx_sounds, intersect(idx_shocks, idx_both)), idx_BA);



nexttile
labels_LA = {['Tone (' num2str(numel(idx_only_LA_sounds)) ')'],...
            ['Shock (' num2str(numel(idx_only_LA_shocks)) ')'],... 
            ['Tone + Shock (' num2str(numel(idx_only_LA_both)) ')'],...
            ['All (' num2str(numel(idx_all_LA)) ')'],...
            ['Non-responseive (' num2str(numel(idx_LA)-numel(unique([idx_only_LA_sounds; idx_only_LA_shocks; idx_only_LA_both; idx_all_LA]))) ')']};

pie([numel(idx_only_LA_sounds) numel(idx_only_LA_shocks) numel(idx_only_LA_both) numel(idx_all_LA) numel(idx_LA)-numel(unique([idx_only_LA_sounds; idx_only_LA_shocks; idx_only_LA_both]))], labels_LA)
title('LA');
nexttile
labels_BA = {['Tone (' num2str(numel(idx_only_BA_sounds)) ')'],...
            ['Shock (' num2str(numel(idx_only_BA_shocks)) ')'],... 
            ['Tone + Shock (' num2str(numel(idx_only_BA_both)) ')'],...
            ['All (' num2str(numel(idx_all_BA)) ')'],...
            ['Non-responseive (' num2str(numel(idx_BA)-numel(unique([idx_only_BA_sounds; idx_only_BA_shocks; idx_only_BA_both; idx_all_BA]))) ')']};

pie([numel(idx_only_BA_sounds) numel(idx_only_BA_shocks) numel(idx_only_BA_both) numel(idx_all_BA) numel(idx_BA)-numel(unique([idx_only_BA_sounds; idx_only_BA_shocks; idx_only_BA_both]))], labels_BA)
title('BA');






% 
% 
% cell_metrics.labels(1:numel(cell_metrics.cellID)) = {'0'};
% cell_metrics.labels(idx_only_LA_sounds_50) = {'idx_only_LA_sounds_50'};
% cell_metrics.labels(idx_only_LA_shocks_50) = {'idx_only_LA_shocks_50'};
% cell_metrics.labels(idx_only_LA_both_50) = {'idx_only_LA_both_50'};
% cell_metrics.labels(idx_all_LA) = {'idx_all_LA'};
% cell_metrics.labels(idx_only_BA_sounds_50) = {'idx_only_BA_sounds_50'};
% cell_metrics.labels(idx_only_BA_shocks_50) = {'idx_only_BA_shocks_50'};
% cell_metrics.labels(idx_only_BA_both_50) = {'idx_only_BA_both_50'};
% cell_metrics.labels(idx_all_BA) = {'idx_all_BA'};
% 
% cell_metrics.labels(1:numel(cell_metrics.cellID)) = {'0'};
% %cell_metrics.labels(intersect(idx_sounds_exc_50,idx_LA)) = {'LA_sounds_exc'};
% cell_metrics.labels(intersect(idx_shocks_exc_50,idx_LA)) = {'LA_shocks_exc'};
% %cell_metrics.labels(intersect(idx_both_exc_50,idx_LA)) = {'LA_both_exc'};
% %cell_metrics.labels(intersect(idx_sounds_exc_50,idx_BA)) = {'BA_sound_exc'};
% cell_metrics.labels(intersect(idx_shocks_exc_50,idx_BA)) = {'BA_shocks_exc'};
% %cell_metrics.labels(intersect(idx_both_exc_50,idx_BA)) = {'BA_both_exc'};
% 
% 
% 
% 
% CellExplorer('metrics',cell_metrics);
