% Gergo_figure_juxta
clear all

cd('C:\Users\dmagyar\Desktop\Gergo\Juxta_results')
load('cell_metrics_juxta.mat')
g.mainFolder = 'C:\Users\dmagyar\Desktop\Gergo';
% 
% g.celldata = readtable('cell_list.xlsx');
% g.networkfolder = 'Z:\HajosLab\Dani\Magyar Dániel\Analysis2\Matlab_files';
% for ii = 1:size(g.celldata,1)
%     load([g.networkfolder '\' g.celldata.filename{ii}])
%     vars = who;                         % get all variable names
%     ttl = vars(contains(vars,'Ch31')); % find variables with 'Ch13'
%     ttl = eval(ttl{1});
%     g.cell_metrics.general.shockTTL{ii} = ttl.times;
%     spikes = vars(contains(vars,'Ch13')); % find variables with 'Ch13'
%     spikes = eval(spikes{1});
%     g.cell_metrics.spikes.times{ii} = spikes.times;
%     disp(ii)
%     clearvars -except g ii 
% end
g.pre_time = 0.1;
g.post_time = 0.1;
g.bin_time = 0.001;
g.smoothvalue = 7;
g.timeaxis = -g.pre_time:g.bin_time:g.post_time;
g.colors = BAfc_colors;
g.test_time = 0.1;
g.exctestvalue = 3; % if the mean change is negative, but there is a peak where zscore is above this, than counted as excited
g.inhtestvalue = 1;
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
g.fontSize2 = 12;
g.fontSize1 = 15;
g.xlinewidth = 2;

g.idx_LA = strcmp(g.cell_metrics.brainRegion, 'LA')';
g.idx_BA = strcmp(g.cell_metrics.brainRegion, 'BA')';


%% Initiate figure
fig = figure('Position', [400, 100, 800, 1200]);
t = tiledlayout(fig,11,8,'TileSpacing', 'tight', 'Padding', 'none');

%% (1,1) - Schematic figure
ax1 = nexttile(t,1,[3 8]);
[img, cmap] = imread([g.mainFolder '\juxta_image.png']);
if ~isempty(cmap)
    img = ind2rgb(img, cmap);  % Convert to RGB
end
imshow(img, 'Parent', ax1);
title(ax1, 'Experimental setup', 'FontSize', g.fontSize1)
clearvars -except t g % clear variables

cellID_1 = 24; 
cellID_2 = 6; 

spike_times{1} = g.cell_metrics.spikes.times{cellID_1}; 
spike_times{2} = g.cell_metrics.spikes.times{cellID_2}; 

stimulus_times{1} = g.cell_metrics.general.shockTTL{cellID_1};
stimulus_times{2} = g.cell_metrics.general.shockTTL{cellID_2};

% Plot raster
time_bins = -g.pre_time:g.bin_time:g.post_time;    
spike_counts = zeros(1, length(time_bins) - 1);   
rasterloc = [25 29];
rastertitles = {'Example LA neuron (Cell 1)','Example BA neuron (Cell 2)'};
for ii = 1:2
    ax = nexttile(t,rasterloc(ii),[2 4]);
    hold on;    
    for trial = 1:length(stimulus_times{ii})
        aligned_spikes = spike_times{ii} - stimulus_times{ii}(trial);
        valid_spikes = aligned_spikes(aligned_spikes >= -g.pre_time & aligned_spikes <= g.post_time);
        scatter(valid_spikes, trial * ones(size(valid_spikes)), 10, 'k', 'filled');
        spike_counts = spike_counts + histcounts(valid_spikes, time_bins);
    end
    xline(0, 'r--', 'LineWidth', 2)
    xlim([-g.pre_time g.post_time])
    % xlabel('Time (s)')
    xticks([-g.pre_time 0 g.post_time]);
    title(rastertitles{ii})
    ylabel('Trials')
    set(gca, 'FontSize', g.fontSize2)
end

% PSTH

psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'shockTTL', 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
resp_spx = sum(psth_spx(:, round(g.pre_time/g.bin_time+1):round((g.pre_time+g.test_time)/g.bin_time)),2);
g.idx_LA(resp_spx<3) = 0;
g.idx_BA(resp_spx<3) = 0;
psth_spx_zscore = zscore(psth_spx,0,2);
psth_spx_zscore = smoothdata(psth_spx_zscore,2,'sgolay', g.smoothvalue);

psth_spx_zscore_LA = psth_spx_zscore(g.idx_LA,:);
psth_spx_zscore_BA = psth_spx_zscore(g.idx_BA,:);

[results, clusters] = Gergo_psth_sorter(g,psth_spx_zscore);

[~, order_LA] = sortrows([results.respDir(g.idx_LA) results.onsetIdx(g.idx_LA)*g.bin_time*1000], [1 2], 'ascend');
[~, order_BA] = sortrows([results.respDir(g.idx_BA) results.onsetIdx(g.idx_BA)*g.bin_time*1000], [1 2], 'ascend');

n_LA_excited = sum(clusters(g.idx_LA) == 1) + sum(clusters(g.idx_LA) == 2);
n_BA_excited = sum(clusters(g.idx_BA) == 1) + sum(clusters(g.idx_BA) == 2);

psth_spx_zscore_LA = psth_spx_zscore_LA(order_LA(1:n_LA_excited),:);
psth_spx_zscore_BA = psth_spx_zscore_BA(order_BA(1:n_BA_excited),:);

g.clim = [min(psth_spx_zscore(:)) max(psth_spx_zscore(:))*2];

% LA heatmap
nexttile(t,41,[2 4]);
imagesc(g.timeaxis,1:size(psth_spx_zscore_LA,1),psth_spx_zscore_LA);
clim(g.clim)
colormap(g.colors.Heatmap);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
xlabel('Time (s)')
xticks([-g.pre_time 0 g.post_time]);
yticks([1 size(psth_spx_zscore_LA,1)])
set(gca, 'FontSize', g.fontSize2);
title('LA neurons', 'FontSize', g.fontSize1)

ylabel('Neuron #')

% BA heatmapt
nexttile(t,45,[2 4]);
imagesc(g.timeaxis,1:size(psth_spx_zscore_BA,1),psth_spx_zscore_BA);
clim(g.clim)
colormap(g.colors.Heatmap);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
xlabel('Time (s)')
xticks([-g.pre_time 0 g.post_time]);
yticks([1 size(psth_spx_zscore_BA,1)])
set(gca, 'FontSize', g.fontSize2);
title('BA neurons', 'FontSize', g.fontSize1)

cb = colorbar('eastoutside', 'FontSize', g.fontSize2);



% clearvars -except g t
%% Bar graph, bin_time changed sometimes

bartimeaxis = -g.pre_time:g.bin_time:g.post_time;
[psth_spx_bar, ~, ~, preAP_norm, postAP_norm] =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'shockTTL', 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
%psth_spx_bar = smoothdata(psth_spx_bar,2,'sgolay', g.smoothvalue);
psth_spx_bar_zscore = zscore(smoothdata(psth_spx_bar,2,'sgolay', g.smoothvalue),0,2);
psth_spx_LA = psth_spx_bar(g.idx_LA,:);
psth_spx_BA = psth_spx_bar(g.idx_BA,:);
psth_spx_LA_zscore = psth_spx_bar_zscore(g.idx_LA,:);
psth_spx_BA_zscore = psth_spx_bar_zscore(g.idx_BA,:);

postAP_spx_LA = postAP_norm(g.idx_LA);
postAP_spx_LA = postAP_spx_LA(order_LA(1:n_LA_excited));
postAP_spx_BA = postAP_norm(g.idx_BA);
postAP_spx_BA = postAP_spx_BA(order_BA(1:n_BA_excited));
postAP_spx_LA = sort(cell2mat(postAP_spx_LA'));
postAP_spx_BA = sort(cell2mat(postAP_spx_BA'));

% LA firing rate histogram
nexttile(t,57,[2 4]);
psth_spx_excited_LA = psth_spx_LA(order_LA(1:n_LA_excited),:);
psth_spx_excited_LA_zscore = psth_spx_LA_zscore(order_LA(1:n_LA_excited),:);
time = bartimeaxis(2:end);
allspikes_LA = sum(psth_spx_excited_LA,1);       % Sum across trials
allspikes_LA_zscore = mean(psth_spx_excited_LA_zscore,1);       % Mean across trials
bar(time, allspikes_LA)
xlabel('Time (s)')
xticks([-g.pre_time 0 g.post_time]);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
ylabel('Z-score')
set(gca, 'FontSize', g.fontSize2);
%title('PFC projecting group', 'FontSize', g.fontSize1);
box off;


% BA firing rate histogram
nexttile(t,61,[2 4]);
psth_spx_excited_BA = psth_spx_BA(order_BA(1:n_BA_excited),:);
psth_spx_excited_BA_zscore = psth_spx_BA_zscore(order_BA(1:n_BA_excited),:);
time = bartimeaxis(2:end);
allspikes_BA = sum(psth_spx_excited_BA,1);       % Sum across trials
allspikes_BA_zscore = mean(psth_spx_excited_BA_zscore,1);       % Mean across trials
bar(time, allspikes_BA)
xlabel('Time (s)')
xticks([-g.pre_time 0 g.post_time]);
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
set(gca, 'FontSize', g.fontSize2);
%title('PFC projecting group', 'FontSize', g.fontSize1);
box off;

LA_neurons_sorted_idx = find(g.idx_LA);
LA_neurons_sorted_idx = LA_neurons_sorted_idx(order_LA(1:n_LA_excited));

BA_neurons_sorted_idx = find(g.idx_BA);
BA_neurons_sorted_idx = BA_neurons_sorted_idx(order_BA(1:n_BA_excited));

%% calculate latencies
g.cell_metrics.t_onset = results.onsetIdx*g.bin_time;

ax = nexttile(t,73,[2 2]);
jitter_amount = 0.1;
jitter_LA = (rand(sum(g.idx_LA), 1) - 0.5) * jitter_amount;
jitter_BA = (rand(sum(g.idx_BA), 1) - 0.5) * jitter_amount;

% Get data and filter out values < 0 or > 0.5
data_LA_all = g.cell_metrics.t_onset(g.idx_LA);
data_BA_all = g.cell_metrics.t_onset(g.idx_BA);

% Filter data
valid_LA = data_LA_all >= 0.005 & data_LA_all <= 0.05;
valid_BA = data_BA_all >= 0.005 & data_BA_all <= 0.05;

data_LA = data_LA_all(valid_LA) * 1000;  % Convert to milliseconds
data_BA = data_BA_all(valid_BA) * 1000;  % Convert to milliseconds

% Update jitter to match filtered data
jitter_LA = jitter_LA(valid_LA);
jitter_BA = jitter_BA(valid_BA);

% Calculate statistics
mean_LA = mean(data_LA);
sd_LA = std(data_LA);
mean_BA = mean(data_BA);
sd_BA = std(data_BA);

% Create scatter plots
scatter(1 + jitter_LA, data_LA)
hold on
scatter(2 + jitter_BA, data_BA)

% Add mean as black dots and SD as whiskers (offset to avoid overlap)
plot(1.1, mean_LA, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
plot([1.1 1.1], [mean_LA - sd_LA, mean_LA + sd_LA], 'k-', 'LineWidth', 2)
% Add whiskers (horizontal caps) at the ends of SD lines
plot([1.08 1.12], [mean_LA - sd_LA, mean_LA - sd_LA], 'k-', 'LineWidth', 2)
plot([1.08 1.12], [mean_LA + sd_LA, mean_LA + sd_LA], 'k-', 'LineWidth', 2)

plot(2.1, mean_BA, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
plot([2.1 2.1], [mean_BA - sd_BA, mean_BA + sd_BA], 'k-', 'LineWidth', 2)
% Add whiskers (horizontal caps) at the ends of SD lines
plot([2.08 2.12], [mean_BA - sd_BA, mean_BA - sd_BA], 'k-', 'LineWidth', 2)
plot([2.08 2.12], [mean_BA + sd_BA, mean_BA + sd_BA], 'k-', 'LineWidth', 2)

% Perform Wilcoxon rank-sum test (non-paired)
[p_value, h] = ranksum(data_LA, data_BA);

% Add significance indicator close to mean values
y_line = (mean_LA + mean_BA) / 2;  % Position text at average of the two means

% Determine significance level and symbol
if p_value < 0.001
    sig_symbol = '***';
elseif p_value < 0.01
    sig_symbol = '**';
elseif p_value < 0.05
    sig_symbol = '*';
else
    sig_symbol = 'ns';
end

text(1.5, y_line, sprintf('%s (p=%.3f)', sig_symbol, p_value), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12)

% Add labels and title
title('Resp. latency')
xlabel('')
ylabel('Time (ms)')
set(gca, 'XTick', [1 2], 'XTickLabel', {'LA', 'BA'})
set(gca, 'FontSize', 12)
ylim([0 60])
xlim([0.8 2.2])

clearvars -except t g postAP_spx_LA postAP_spx_BA valid_LA valid_BA
%% response magnitude

ax = nexttile(t,75,[2 2]);
jitter_amount = 0.1;
jitter_LA = (rand(sum(g.idx_LA), 1) - 0.5) * jitter_amount;
jitter_BA = (rand(sum(g.idx_BA), 1) - 0.5) * jitter_amount;
psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'shockTTL', 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
response_spikes = sum(psth_spx(:, round(g.pre_time/g.bin_time+23):round(g.pre_time/g.bin_time+34)),2);
response_spikes = response_spikes./ cellfun(@numel,g.cell_metrics.general.shockTTL);
data_LA_all = response_spikes(g.idx_LA);
data_BA_all = response_spikes(g.idx_BA);
data_LA = data_LA_all(valid_LA);  % Convert to milliseconds
data_BA = data_BA_all(valid_BA);  % Convert to milliseconds

% Update jitter to match filtered data
jitter_LA = jitter_LA(valid_LA);
jitter_BA = jitter_BA(valid_BA);

% Calculate statistics
mean_LA = mean(data_LA);
sd_LA = std(data_LA);
mean_BA = mean(data_BA);
sd_BA = std(data_BA);

% Create scatter plots
scatter(1 + jitter_LA, data_LA)
hold on
scatter(2 + jitter_BA, data_BA)

% Add mean as black dots and SD as whiskers (offset to avoid overlap)
plot(1.1, mean_LA, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
plot([1.1 1.1], [mean_LA - sd_LA, mean_LA + sd_LA], 'k-', 'LineWidth', 2)
% Add whiskers (horizontal caps) at the ends of SD lines
plot([1.08 1.12], [mean_LA - sd_LA, mean_LA - sd_LA], 'k-', 'LineWidth', 2)
plot([1.08 1.12], [mean_LA + sd_LA, mean_LA + sd_LA], 'k-', 'LineWidth', 2)

plot(2.1, mean_BA, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
plot([2.1 2.1], [mean_BA - sd_BA, mean_BA + sd_BA], 'k-', 'LineWidth', 2)
% Add whiskers (horizontal caps) at the ends of SD lines
plot([2.08 2.12], [mean_BA - sd_BA, mean_BA - sd_BA], 'k-', 'LineWidth', 2)
plot([2.08 2.12], [mean_BA + sd_BA, mean_BA + sd_BA], 'k-', 'LineWidth', 2)

% Perform Wilcoxon rank-sum test (non-paired)
[p_value, h] = ranksum(data_LA, data_BA);

% Add significance indicator close to mean values
y_line = (mean_LA + mean_BA) / 2;  % Position text at average of the two means

% Determine significance level and symbol
if p_value < 0.001
    sig_symbol = '***';
elseif p_value < 0.01
    sig_symbol = '**';
elseif p_value < 0.05
    sig_symbol = '*';
else
    sig_symbol = 'ns';
end

text(1.5, y_line, sprintf('%s (p=%.3f)', sig_symbol, p_value), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12)

% Add labels and title
title('Resp. magnitude')
xlabel('')
ylabel('n spikes/trial')
set(gca, 'XTick', [1 2], 'XTickLabel', {'LA', 'BA'})
set(gca, 'FontSize', 12)
ylim([-1 5])
xlim([0.8 2.2])


clearvars -except t g postAP_spx_LA postAP_spx_BA valid_LA valid_BA

%% caluclate baseline firing rate
excludesec = 1;
for ii = 1:size(g.cell_metrics.spikes.times,2)
    spikes = g.cell_metrics.spikes.times{ii};
    stimulus = g.cell_metrics.general.shockTTL{ii};
    for jj = 1:size(stimulus)
        spikes(spikes>stimulus(jj) & spikes<stimulus(jj)+excludesec) = [];
    end
    spikes = sort(spikes);
    n_spikes_baseline = numel(spikes);
    t_active = spikes(end) - spikes(1);
    baseline_fr(ii) = n_spikes_baseline/t_active;
end

ax = nexttile(t,77,[2 2]);
jitter_amount = 0.1;
jitter_LA = (rand(sum(g.idx_LA), 1) - 0.5) * jitter_amount;
jitter_BA = (rand(sum(g.idx_BA), 1) - 0.5) * jitter_amount;
data_LA_all = baseline_fr(g.idx_LA);
data_BA_all = baseline_fr(g.idx_BA);
data_LA = data_LA_all(valid_LA);  % Convert to milliseconds
data_BA = data_BA_all(valid_BA);  % Convert to milliseconds

% Update jitter to match filtered data
jitter_LA = jitter_LA(valid_LA);
jitter_BA = jitter_BA(valid_BA);

data_LA(data_LA==Inf) = 0;
data_BA(data_BA==Inf) = 0;
% Calculate statistics
mean_LA = mean(data_LA);
sd_LA = std(data_LA);
mean_BA = mean(data_BA);
sd_BA = std(data_BA);

% Create scatter plots
scatter(1 + jitter_LA, data_LA)
hold on
scatter(2 + jitter_BA, data_BA)

% Add mean as black dots and SD as whiskers (offset to avoid overlap)
plot(1.1, mean_LA, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
plot([1.1 1.1], [mean_LA - sd_LA, mean_LA + sd_LA], 'k-', 'LineWidth', 2)
% Add whiskers (horizontal caps) at the ends of SD lines
plot([1.08 1.12], [mean_LA - sd_LA, mean_LA - sd_LA], 'k-', 'LineWidth', 2)
plot([1.08 1.12], [mean_LA + sd_LA, mean_LA + sd_LA], 'k-', 'LineWidth', 2)

plot(2.1, mean_BA, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
plot([2.1 2.1], [mean_BA - sd_BA, mean_BA + sd_BA], 'k-', 'LineWidth', 2)
% Add whiskers (horizontal caps) at the ends of SD lines
plot([2.08 2.12], [mean_BA - sd_BA, mean_BA - sd_BA], 'k-', 'LineWidth', 2)
plot([2.08 2.12], [mean_BA + sd_BA, mean_BA + sd_BA], 'k-', 'LineWidth', 2)

% Perform Wilcoxon rank-sum test (non-paired)
[p_value, h] = ranksum(data_LA, data_BA);

% Add significance indicator close to mean values
y_line = (mean_LA + mean_BA) / 2;  % Position text at average of the two means

% Determine significance level and symbol
if p_value < 0.001
    sig_symbol = '***';
elseif p_value < 0.01
    sig_symbol = '**';
elseif p_value < 0.05
    sig_symbol = '*';
else
    sig_symbol = 'ns';
end

text(1.5, y_line, sprintf('%s (p=%.3f)', sig_symbol, p_value), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12)

% Add labels and title
title('Baseline fr.')
xlabel('')
ylabel('Hz')
set(gca, 'XTick', [1 2], 'XTickLabel', {'LA', 'BA'})
set(gca, 'FontSize', 12)
ylim([-5 20])
xlim([0.8 2.2])





%% cumulative diagram
ax = nexttile(t,79,[2 2]);
% allspikes_LA_ks= allspikes_LA(123:134);
% allspikes_BA_ks= allspikes_BA(123:134);

allspikes_LA_ks= postAP_spx_LA(postAP_spx_LA>0.023& postAP_spx_LA<0.034);
allspikes_BA_ks= postAP_spx_BA(postAP_spx_BA>0.023& postAP_spx_BA<0.034);

% Plot cumulative frequency diagrams with significance testing
[f1, x1] = ecdf(allspikes_LA_ks);
[f2, x2] = ecdf(allspikes_BA_ks);

plot(x1, f1, 'b-', 'LineWidth', 2, 'DisplayName', 'LA');
hold on;
plot(x2, f2, 'r-', 'LineWidth', 2, 'DisplayName', 'BA');

xlabel('Spike latency');
ylabel('Cumulative Prob.');
title('Cumul. Freq. Comp.');
legend('show', 'Location','northwest');
grid on;
xlim([0 0.05])
% Statistical significance testing
% Kolmogorov-Smirnov test for two samples
[h, p, ks_stat] = kstest2(allspikes_LA_ks, allspikes_BA_ks);

text(0.03, 0.2, sprintf('%s (p=%.3f)', sig_symbol, p), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12)

set(gca, 'FontSize', 12)

