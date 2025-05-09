% function BAfc_figure_2_v3
% % Plotting region specific responses on heatmaps, proper clustering
clear all
recordings = {...
    'MD292_002_kilosort',...
    'MD293_kilosort',...
    'MD294_kilosort',...
    'MD295_kilosort',...
    'MD296_kilosort',...
    'MD297_kilosort'};
g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'});

clearvars -except g
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 12;
g.fontSize2 = 12;
g.fontsize3 = 10;
g.pre_time = 5;
g.post_time = 5;
g.test_time = 0.5; % 1 sec is good, beacuse captures better than 0.5
g.testvalue =3; % if the mean change is negative, but there is a peak where zscore is above this, than counted as excited
g.bin_time = 0.005; % 0.01 volt sokat 5 os smoothal
g.smoothvalue = 5;
g.plotwin = [0.05 0.5];
g.spxwin = [0.02 0.1];
g.spxbin = 0.002;
g.timeaxis = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.clim = [-5.5 13];
roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;

%% Initiate figure
fig = figure('Position', [400, 100, 1400, 1200]);
t = tiledlayout(fig,5,4,'TileSpacing', 'compact', 'Padding', 'none');
%% (1,1:4) - Heatmaps
mploc = [1 2 3 4];
ttl = {'triptest_sound_only', 'triptest_sound_only',  'triptest_shocks_only', 'triptest_shocks_only'};
br = {'LA', 'BA', 'LA', 'BA'};
hmptitles = {'LA population, CS reponse', 'BA population, CS response', 'LA population, US reponses', 'BA population, US response'};
for hmp = 1:4 
    ax = nexttile(t,mploc(hmp),[2 1]);
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    psth_spx = zscore(psth_spx,0,2);   
    psth_spx = smoothdata(psth_spx,2,'sgolay', g.smoothvalue);
    idx_PN = strcmp(g.cell_metrics.brainRegion,br{hmp}) & strcmp(g.cell_metrics.putativeCellType,'PN');
    idx_IN = strcmp(g.cell_metrics.brainRegion,br{hmp}) & strcmp(g.cell_metrics.putativeCellType,'IN');   

    for ii = 1:size(psth_spx,1)
        datasegment = psth_spx(ii,roi);
        if all(datasegment == datasegment(1)) ||  max(datasegment)<g.testvalue
            results.onsetIdx(ii,1) = -1;
            results.offsetIdx(ii,1) = length(datasegment); 
        else
            results.onsetIdx(ii,1) = find(datasegment>=g.testvalue,1,'first'); % begining of the response
            results.offsetIdx(ii,1) = find(datasegment>=g.testvalue,1,'last'); % end of response
            responseMean = mean(datasegment(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
            k = 0;
            while k == 0
                if responseMean < g.testvalue/2
                    datasegment(results.offsetIdx(ii,1):end) = NaN;
                    results.offsetIdx(ii,1) = find(datasegment>=g.testvalue,1,'last'); % end of response
                    responseMean = mean(datasegment(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
                else 
                    k = 1;
                end
            end
        end
    end
    % heatmaps
    psth_spx_PN = psth_spx(idx_PN,:);
    psth_spx_IN = psth_spx(idx_IN,:);
    edges = 0:100:g.test_time*1000;
    onsetTime_d = discretize(results.onsetIdx*g.bin_time*1000, edges); % round to 10 ms values.
    offsetTime_d = discretize(results.offsetIdx*g.bin_time*1000, edges); % round to 10 ms values.
    durationTime_d = discretize((results.offsetIdx-results.onsetIdx)*g.bin_time*1000, edges); % round to 10 ms values.

    clear clusters_PN psth_spx
    clusters(onsetTime_d == 1 & durationTime_d == 1,1) = 1;
    clusters(onsetTime_d == 1 & durationTime_d > 1,1) = 2;
    clusters(onsetTime_d > 1 ,1) = 2;
    clusters(isnan(onsetTime_d)) = 3;
    
    [~, order_PN] = sortrows([clusters(idx_PN) results.onsetIdx(idx_PN)*g.bin_time*1000], [1 2], 'ascend');

    % if hmp == 3
    %     order_PN = hmpdata.order_1.PN;
    % end
    % 

    psth_spx = [psth_spx_PN(order_PN,:)];
    matrix = psth_spx(:,(g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    imagesc(g.timeaxis,1:size(matrix,1),matrix);
    title(hmptitles{hmp}, 'FontSize', g.fontSize2)
    clim(g.clim)
    disp([-max(max(psth_spx,[],1))/2.5 max(max(psth_spx,[],1))])
    
    colormap(g.colors.Heatmap); 
    yline(sum(idx_PN)+0.5, 'k', 'LineWidth', 1);
    xlabel('Time (s)', 'FontSize', g.fontSize2)
    if hmp ==1 
        ylabel('Cell number', 'FontSize', g.fontSize2)
    end
    if hmp == 4
        cb = colorbar('eastoutside', 'FontSize',15);
    end

    % plotting x symbols at onset and offset
    onset_PN = results.onsetIdx(idx_PN);
    offset_PN = results.offsetIdx(idx_PN);
    onsets_sorted = [onset_PN(order_PN)];
    offsets_sorted = [offset_PN(order_PN)];
    hold on
    for ii = 1:size(matrix,1)
        plot(onsets_sorted(ii)*g.bin_time, ii, 'x', 'Color', 'r', 'MarkerSize', 4)
        plot(offsets_sorted(ii)*g.bin_time, ii, 'x', 'Color', 'b', 'MarkerSize', 4)
    end

    hmpdata.(['hmp_' num2str(hmp)]) = clusters;
    hmpdata.(['order_' num2str(hmp)]).PN = order_PN;
end
clearvars -except g t hmpdata

%% (3:4,1:4) - Population spikes with LFP
ttl = {'triptest_sound_only', 'triptest_sound_only',  'triptest_shocks_only', 'triptest_shocks_only'};
br = {'LA', 'BA', 'LA', 'BA'};

for ii = 1:4   
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{ii}, 'pre_time', g.plotwin(1), 'post_time', g.plotwin(2), 'bin_time', g.bin_time);
    time_axis = linspace(-g.plotwin(1),g.plotwin(2),size(psth_spx,2));   
    ax = nexttile;
    idx =  strcmp(g.cell_metrics.brainRegion,br{ii}) & strcmp(g.cell_metrics.putativeCellType,'PN');
    spx = sum(psth_spx(idx,:),1);
    spx = smoothdata(sum(psth_spx(idx,:),1), 'sgolay', g.smoothvalue);
    bar(time_axis,spx, 'k'); % 'k' for black bars
    ylim([0 150])
    xlabel('Time (s)', 'FontSize', g.fontSize2)
    if ii == 1
        ylabel('Number of spikes', 'FontSize', g.fontSize2)  
    end
end

clearvars -except g t hmpdata

%% (3:4,1:4) - cluster spikes
ttl = {'triptest_sound_only', 'triptest_sound_only',  'triptest_shocks_only', 'triptest_shocks_only'};
br = {'LA', 'BA', 'LA', 'BA'};
hmptitles = {'LA population, CS reponse', 'LA population, US response', 'BA population, CS reponses', 'BA population, US response'};
layoutTile = [13 14 15 16];

for ii = 1:4
    ts = tiledlayout(t,2,1);
    ts.Layout.Tile = layoutTile(ii);
    ts.Layout.TileSpan = [2 1];

    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{ii}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    time_axis = linspace(-g.plotwin(1),g.plotwin(2),(g.plotwin(1)+g.plotwin(2))/g.bin_time);
    idx_short = hmpdata.(['hmp_' num2str(ii)]) == 1;
    idx_long = hmpdata.(['hmp_' num2str(ii)]) == 2;

    idx_PN = strcmp(g.cell_metrics.brainRegion,br{ii}) & strcmp(g.cell_metrics.putativeCellType,'PN');
    idx_IN = strcmp(g.cell_metrics.brainRegion,br{ii}) & strcmp(g.cell_metrics.putativeCellType,'IN');   

    ax1 = nexttile(ts,[1 1]);
    idx =  idx_PN' & idx_short;
    spx = sum(psth_spx(idx,:),1);
    spx = smoothdata(sum(psth_spx(idx,:),1), 'sgolay', g.smoothvalue);

    spx = spx(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    bar(time_axis,spx, 'k'); % 'k' for black bars
    ylim([0 50])
    xlabel('Time (s)', 'FontSize', g.fontSize2)
    if ii == 1
        ylabel('Number of spikes', 'FontSize', g.fontSize2)
    end
    
    ax2 = nexttile(ts,[1 1]);
    idx =  idx_PN' & idx_long;
    spx = sum(psth_spx(idx,:),1);
    spx = smoothdata(sum(psth_spx(idx,:),1), 'sgolay', g.smoothvalue);

    spx = spx(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    bar(time_axis,spx, 'k'); % 'k' for black bars
    ylim([0 90])
    xlabel('Time (s)', 'FontSize', g.fontSize2)
    if ii == 1
        ylabel('Number of spikes', 'FontSize', g.fontSize2)
    end
end