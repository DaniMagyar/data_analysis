function [frame] = BAfc_plot_spikes_and_LFP(varargin)

prs = inputParser;
addParameter(prs,'cell_metrics',[],@isstruct)
addParameter(prs,'ttl',@ischar)
addParameter(prs,'animal',[],@iscell) % selected animals from cell_metrics
addParameter(prs,'twin',[-0.1 0.3],@isnumeric)
addParameter(prs,'visible','on',@ischar)
addParameter(prs,'response','exc',@ischar) % plot excited or inhiboted neurons' response
parse(prs,varargin{:})
g = prs.Results;


cell_metrics = g.cell_metrics;

LFP_LA = [];
LFP_BA = [];
if ~isempty(g.animal)
    fnames = g.animal;
else
    fnames = fieldnames(cell_metrics.LFP);
end
for ii = 1:numel(fnames)
    indicesLA = find(strcmp(cell_metrics.LFP.(fnames{ii}).(g.ttl).channelBR, 'LA'));
    temp(1:numel(indicesLA),:) =  cell_metrics.LFP.(fnames{ii}).(g.ttl).data(indicesLA,:);
    LFP_LA = [LFP_LA;temp]; clear temp;
    indicesBA = find(strcmp(cell_metrics.LFP.(fnames{ii}).(g.ttl).channelBR, 'BA'));
    temp(1:numel(indicesBA),:) =  cell_metrics.LFP.(fnames{ii}).(g.ttl).data(indicesBA,:);
    LFP_BA = [LFP_BA;temp]; clear temp;
end
LFPS{1} = mean(LFP_LA,1);
LFPS{2} = mean(LFP_BA,1);
idx_LA = find(contains(cell_metrics.brainRegion, 'LA'));
idx_BA = find(contains(cell_metrics.brainRegion, 'BA'));
% Find excited neurons
[psth_spx, res_1, ~] =          BAfc_find_response('cell_metrics', cell_metrics, 'ttl', g.ttl, 'pre_time', 0.5, 'post_time', 0.5, 'test_win', [0.013 0.027], 'artefactLength',     0, 'psth_out', g.twin);
[~, res_2, ~] =                 BAfc_find_response('cell_metrics', cell_metrics, 'ttl', g.ttl, 'pre_time', 0.5, 'post_time', 0.5, 'test_win', [0.028 0.039], 'artefactLength',     0);
[~, res_3, ~] =                 BAfc_find_response('cell_metrics', cell_metrics, 'ttl', g.ttl, 'pre_time', 0.5, 'post_time', 0.5, 'test_win', [0.040 0.092], 'artefactLength',     0);
idx_all_exc = unique([find(strcmp(res_1,'exc')); find(strcmp(res_2,'exc')); find(strcmp(res_3,'exc'))]);
idx_all_inh = unique([find(strcmp(res_1,'inh')); find(strcmp(res_2,'inh')); find(strcmp(res_3,'inh'))]);
if strcmp(g.response,'exc')
    idx_all = idx_all_exc;
elseif strcmp(g.response,'inh')
    idx_all = idx_all_exc;
elseif strcmp(g.response,'all')
    idx_all = unique([idx_all_exc; idx_all_inh]);
end

if ~isempty(g.animal)
    idx_included = zeros(size(cell_metrics.cellID))';
    for ii = 1:size(g.animal,2)
        idx_included(strcmp(cell_metrics.animal, g.animal{ii})) = 1;
    end
    idx_all = intersect(idx_all, find(idx_included == 1));
    idx_LA = intersect(idx_LA, find(idx_included ==1));
    idx_BA = intersect(idx_BA, find(idx_included ==1));
end

% LA sounds
idx_LA_resp = intersect(idx_all, idx_LA);
% BA sounds
idx_BA_resp = intersect(idx_all, idx_BA);

for ii = 1:size(psth_spx,1)
    psth_spx(ii,:) = smoothdata(psth_spx(ii,:), 'movmean', 10);
end

onPlot = {{idx_LA, idx_BA}, {idx_LA_resp, idx_BA_resp}};
onPlot_names = {{'LA all', 'BA all'},{'LA responsive', 'BA responsive'}};
fig = figure('Units', 'pixels', 'Position', [100, 100, 1200, 900], 'Visible', g.visible);
t = tiledlayout(2,2);
title(t, strjoin(string(g.animal)))
for oP = 1:numel(onPlot)
    for ii = 1:numel(onPlot)
        nexttile
        selected_rows = psth_spx(onPlot{oP}{ii}, :);
        % Compute row-wise sum
        row_sums = sum(selected_rows, 1);
        row_sums = row_sums/mean(row_sums(1:abs(g.twin(1)/0.001)));
        bin_size = 1;
        row_vector = row_sums;
        % Ensure the row length is even (truncate if necessary)
        num_bins = floor(size(row_vector, 2) / bin_size) * bin_size;
        row_vector = row_vector(1:num_bins); % Trim to ensure even binning
        % Reshape and sum adjacent elements
        binned_row(ii,:) = sum(reshape(row_vector, bin_size, []), 1);
        time_axis = linspace(g.twin(1), g.twin(2), size(binned_row(ii,:),2))*1000;
        % Plot histogram (bar plot)
        yyaxis left
        bar(time_axis, binned_row(ii,:), 'k'); % 'k' for black bars
        xlabel('Time (ms)');
        %ylim([0 max(binned_row(1,:)+0.5)])
        ylabel('Normalized Spike Count', 'Color', 'k');
        title([onPlot_names{oP}{ii} ', n = ' num2str(size(selected_rows,1))]);
        grid on;
        hold on
        yyaxis right
        plot(time_axis, LFPS{ii})
        %ylim([-500 2500])
        hold off       
    end
end
drawnow;
frame = getframe(fig);
