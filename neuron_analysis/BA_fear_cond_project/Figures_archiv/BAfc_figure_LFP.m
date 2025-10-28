
% % Plotting region specific responses on heatmaps, proper clustering
clear all
recordings = {...
    'MD307_kilosort'};
g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', {'triptest_sound_only','triptest_sound_only_light', 'triptest_shocks_only', 'triptest_both'});
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);
g.plotwin = [0.2 0.4];
g.cell_metrics = BAfc_load_LFP('cell_metrics',g.cell_metrics, 'ttl', {'triptest_sound_only_light'},'twin', [-g.plotwin(1) g.plotwin(2)]);


clearvars -except g t hmpdata
g.bin_time = 0.005; % 0.01 volt sokat 5 os smoothal
g.smoothvalue = 5;
g.plotwin = [0.2 0.5];
g.fontSize2 = 12;
%% (3:4,1:4) - Population spikes with LFP
ttl = {'triptest_sound_only_light', 'triptest_sound_only_light'};
br = {'LA', 'LA'};

fig = figure('Position', [400, 100, 1400, 800]);
t = tiledlayout(fig,1,2,'TileSpacing', 'compact', 'Padding', 'none');


for ii = 1:2   
    LFP = [];
    fnames = fieldnames(g.cell_metrics.LFP);
    for fn = 1:numel(fnames)       
        region_indices = find(strcmp(g.cell_metrics.LFP.(fnames{fn}).(ttl{ii}).channelBR, br{ii}));
        temp(1:numel(region_indices),:) =  g.cell_metrics.LFP.(fnames{fn}).(ttl{ii}).data(region_indices,:);
        LFP = [LFP;temp]; clear temp;
    end
    LFP_mean = mean(LFP,1);
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{ii}, 'pre_time', g.plotwin(1), 'post_time', g.plotwin(2), 'bin_time', g.bin_time);
    time_axis = linspace(-g.plotwin(1),g.plotwin(2),size(psth_spx,2));   
    ax = nexttile;
    if ii == 1
        idx =  strcmp(g.cell_metrics.brainRegion,br{ii}) & strcmp(g.cell_metrics.putativeCellType,'PN');
    elseif ii == 2
        idx =  strcmp(g.cell_metrics.brainRegion,br{ii}) & strcmp(g.cell_metrics.putativeCellType,'IN');
    end
    spx = sum(psth_spx(idx,:),1);
    %spx = smoothdata(sum(psth_spx(idx,:),1), 'sgolay', g.smoothvalue);
    bar(time_axis,spx, 'k'); % 'k' for black bars
    %ylim([0 50])
    xlabel('Time (s)', 'FontSize', g.fontSize2)
    if ii == 1
        ylabel('Number of spikes', 'FontSize', g.fontSize2)  
        title('LA PN')
    elseif ii == 2
        title ('LA IN')
    end
    hold on
    yyaxis right
    plot(linspace(time_axis(1),time_axis(end),size(LFP_mean,2)), LFP_mean)
    hold off       
end

clearvars -except g t hmpdata