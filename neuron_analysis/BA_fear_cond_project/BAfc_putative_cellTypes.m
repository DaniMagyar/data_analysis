function [cell_metrics] = BAfc_putative_cellTypes(varargin)

% Hardcoded script for BAfc project's waweform plot.

prs = inputParser;
addParameter(prs,'cell_metrics',[],@isstruct)
addParameter(prs,'plot',false,@islogical) % plot results
addParameter(prs,'width_critical',0.475,@isnumeric) % border between narrow and wide waveforms
addParameter(prs,'fr_critical',6,@isnumerci) % border between fast and regular firing rate
parse(prs,varargin{:})
g = prs.Results;

cell_metrics = g.cell_metrics;
colors = BAfc_colors;

%cell_metrics = BAfc_load_neurons('recordings', g.recordings);
timeaxis = cell_metrics.waveforms.time{1};
narrow_idx = find(cell_metrics.troughToPeak<=g.width_critical);
wide_idx = find(cell_metrics.troughToPeak>g.width_critical);
idx_fs = find(cell_metrics.firingRate>=g.fr_critical);
idx_rs = find(cell_metrics.firingRate<g.fr_critical);
idx_LA = find(contains(cell_metrics.brainRegion, 'LA'));
idx_BA = find(contains(cell_metrics.brainRegion, 'BA'));

idx_IN = intersect(intersect(narrow_idx, idx_fs),[idx_LA idx_BA]);
idx_PN = intersect(intersect(wide_idx, idx_rs),[idx_LA idx_BA]);
idx_unknown = setdiff([idx_LA idx_BA], [idx_IN idx_PN]);

cell_metrics.putativeCellType(idx_IN) = {'IN'};
cell_metrics.putativeCellType(idx_PN) = {'PN'};
cell_metrics.putativeCellType(idx_unknown) = {'unknown'};

if g.plot
    fig = figure('Position', [900, 400, 800, 800]);
    tiledlayout(fig,3,2,'TileSpacing', 'compact', 'Padding', 'compact');
    %% plot example neurons
    % narrow example
    nexttile % waveform
    cellID_1 = intersect(find(strcmp(cell_metrics.animal, 'MD250')), find(cell_metrics.cluID == 321));
    plot(timeaxis(13:48),cell_metrics.waveforms.filt{cellID_1}(13:48), 'Color', colors.IN_primary, 'LineWidth', 2);  
    xlabel('Time (ms)', 'FontSize', 12);
    ylim([-400 120])
    ylabel('Voltage \muV','interpreter','Tex', 'FontSize', 12)
    title('Example IN waveform', 'FontSize', 12)
    
    % wide example
    nexttile % waveform
    cellID_2 = intersect(find(strcmp(cell_metrics.animal, 'MD289')), find(cell_metrics.cluID == 363));
    plot(timeaxis(13:48),cell_metrics.waveforms.filt{cellID_2}(13:48), 'Color', colors.PN_primary, 'LineWidth', 2);  
    xlabel('Time (ms)', 'FontSize', 12);
    ylim([-400 120])
    title('Example PN waveform', 'FontSize', 12)

    spike_timestamps_n = cell_metrics.spikes.times{cellID_1}; % example neuron
    ISI_n = diff(spike_timestamps_n); 
    spike_timestamps_w = cell_metrics.spikes.times{cellID_2}; % example neuron
    ISI_w = diff(spike_timestamps_w);
    num_bins = 50; % Increase the number of bins
    bin_edges = logspace(log10(min([ISI_n;ISI_w])), log10(max([ISI_n;ISI_w])), num_bins);    

    nexttile % narrow ISI   
    histogram(ISI_n, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', colors.IN_primary);  
    xlim([0.001 100])
    set(gca, 'XScale', 'log'); % Set x-axis to log scale
    xlabel('Time (ms)', 'FontSize', 12);
    ylim([0 0.1])    
    ylabel('Probability', 'FontSize', 12);
    title('Example IN ISI', 'FontSize', 12)

    nexttile % wide ISI
    histogram(ISI_w, 'BinEdges', bin_edges, 'Normalization', 'probability', 'FaceColor', colors.PN_primary);
    xlim([0.001 100])
    set(gca, 'XScale', 'log'); % Set x-axis to log scale
    xlabel('Time (ms)', 'FontSize', 12);
    ylim([0 0.1])
    title('Example PN ISI', 'FontSize', 12)

    %% Plot all trough-to-peak values against FR
    nexttile()
    hold on
    plot(cell_metrics.troughToPeak(idx_IN),cell_metrics.firingRate(idx_IN), 'o', 'Color',colors.IN_primary)
    plot(cell_metrics.troughToPeak(idx_PN),cell_metrics.firingRate(idx_PN), 'o', 'Color',colors.PN_primary)
    plot(cell_metrics.troughToPeak(idx_unknown),cell_metrics.firingRate(idx_unknown), 'o', 'Color',[.7 .7 .7])
    set(gca, 'YScale', 'log')  % Makes the Y-axis logarithmic
    yline(6, '--k', 'LineWidth', 1);
    xline(0.425, '--k', 'LineWidth', 1);   
    ylabel('Firing rate (Hz)', 'FontSize', 12);
    xlabel('Trough to peak time (ms)', 'FontSize', 12);
    title('Spike features accross unit type', 'FontSize', 12);

    %% plot narrow and wide spike waveforms separately
    nexttile()
    waveforms = zeros(size(cell_metrics.waveforms.filt,2),48);
    for ii = 1:size(cell_metrics.waveforms.filt,2)
        if size(cell_metrics.waveforms.time{ii},2) == 48
            waveforms(ii,:) = zscore(cell_metrics.waveforms.filt{ii});
        elseif size(cell_metrics.waveforms.time{ii},2) == 96
            waveforms(ii,:) = zscore(cell_metrics.waveforms.filt{ii}(30:77)); % 30:77 is the best, others shift towrads left or right
        end
    end
    % flipping reverse spikes
    for ii = 1:size(waveforms,1)
        if median(waveforms(ii,23:27)) > 1 && strcmp(cell_metrics.putativeCellType{ii}, 'Narrow Interneuron')
            waveforms(ii,:) = -waveforms(ii,:);
        end
    end
    % finding outliers based on Euclidean distance
    mean_wide = mean(waveforms(wide_idx,:),1);
    mean_narrow = mean(waveforms(narrow_idx,:),1);
    
    distances = zeros(size(waveforms,1), 1);
    for ii = wide_idx
        distances(ii) = sqrt(sum((waveforms(ii, :) - mean_wide).^2));
    end
    for ii = narrow_idx
        distances(ii) = sqrt(sum((waveforms(ii, :) - mean_narrow).^2));
    end    
    waveforms(distances>7,:) = NaN;    
    % Compute mean and standard deviation
    mean_IN = mean(waveforms(idx_IN, 13:48), 1, 'omitnan');
    std_IN = std(waveforms(idx_IN, 13:48), 0, 1, 'omitnan');
    mean_PN = mean(waveforms(idx_PN, 13:48), 1, 'omitnan');
    std_PN = std(waveforms(idx_PN, 13:48), 0, 1, 'omitnan');
    % Plot mean waveform with shaded standard deviation for IN
    hold on;
    fill([timeaxis(13:48), fliplr(timeaxis(13:48))], ...
         [mean_IN + std_IN, fliplr(mean_IN - std_IN)], ...
         colors.IN_secondary, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shaded region for SD
    plot(timeaxis(13:48), mean_IN, 'Color', colors.IN_primary, 'LineWidth', 1.5); % Mean waveform
    hold off;
    % Plot mean waveform with shaded standard deviation for IN
    %figure;
    hold on;
    fill([timeaxis(13:48), fliplr(timeaxis(13:48))], ...
         [mean_PN + std_PN, fliplr(mean_PN - std_PN)], ...
         colors.PN_secondary, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shaded region for SD
    plot(timeaxis(13:48), mean_PN, 'Color', colors.PN_primary, 'LineWidth', 1.5); % Mean waveform
    hold off;     
    xlabel('Time (ms)', 'FontSize', 12);
    ylabel('Z-score', 'FontSize', 12)
    title('Normalized waveforms accross unit type', 'FontSize', 12);
end   