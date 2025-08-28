% function BAfc_figure_2_v3
% % Plotting region specific responses on heatmaps, proper clustering
clear all
% recordings = {...
%     'MD292_002_kilosort',...
%     'MD293_kilosort',...
%     'MD294_kilosort',...
%     'MD295_kilosort',...
%     'MD296_kilosort',...
%     'MD297_kilosort',...
%     'MD298_kilosort',...
%     'MD299_kilosort',...
%     'MD300_kilosort',...
%     'MD304_kilosort'};

recordings = {...
    'MD292_002_kilosort',...
    'MD293_kilosort',...
    'MD294_kilosort',...
    'MD295_kilosort',...
    'MD296_kilosort',...
    'MD297_kilosort',...
    'MD298_kilosort',...
    'MD299_kilosort',...
    'MD300_kilosort',...
    'MD304_kilosort',...
    'MD307_kilosort',...
    'MD309_kilosort',...
    'MD311_kilosort',...
    'MD312_kilosort',...
    'MD313_kilosort',...
    'MD314_kilosort'};

g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'});

clearvars -except g
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 13;
g.fontSize2 = 12;
g.xlinewidth = 2;
g.pre_time = 5;
g.post_time = 5;
g.test_time = 0.5; % 1 sec is good, beacuse captures better than 0.5
g.exctestvalue =3; % if the mean change is negative, but there is a peak where zscore is above this, than counted as excited
g.inhtestvalue = 1;
g.bin_time = 0.005; % 0.01 volt sokat 5 os smoothal
g.smoothvalue = 5;
g.plotwin = [0.2 0.5];
g.spxwin = [0.02 0.1];
g.spxbin = 0.002;
g.timeaxis = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.clim = [-5.5 13];
roi = round((g.pre_time+0.01)/g.bin_time):(g.pre_time+g.test_time)/g.bin_time;

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
        if all(datasegment == datasegment(1))
            results.onsetIdx(ii,1) = -1;
            results.offsetIdx(ii,1) = length(datasegment); 
            if datasegment(1)>=g.exctestvalue
                results.respDir(ii,1) = -1;
            elseif datasegment(1)<=-g.inhtestvalue
                 results.respDir(ii,1) = 1;
            else
                 results.respDir(ii,1) = 0;
                 results.onsetIdx(ii,1) = NaN;
                 results.offsetIdx(ii,1) = NaN;
            end
        elseif max(datasegment)>=g.exctestvalue
            results.respDir(ii,1) = -1;
            results.onsetIdx(ii,1) = find(datasegment>=g.exctestvalue,1,'first'); % begining of the response
            results.offsetIdx(ii,1) = find(datasegment>=g.exctestvalue,1,'last'); % end of response
            responseMean = mean(datasegment(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
            k = 0;
            while k == 0
                if responseMean < g.exctestvalue/2
                    datasegment(results.offsetIdx(ii,1):end) = NaN;
                    results.offsetIdx(ii,1) = find(datasegment>=g.exctestvalue,1,'last'); % end of response
                    responseMean = mean(datasegment(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
                else 
                    k = 1;
                end
            end
        elseif max(datasegment)<g.exctestvalue
            datasegment_inverse = -smoothdata(datasegment);
            if max(datasegment_inverse)>=g.inhtestvalue
                results.respDir(ii,1) = 1;
                results.onsetIdx(ii,1) = find(datasegment_inverse>=g.inhtestvalue,1,'first'); % begining of the response
                results.offsetIdx(ii,1) = find(datasegment_inverse>=g.inhtestvalue,1,'last'); % end of response
                responseMean = mean(datasegment_inverse(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
                k = 0;
                while k == 0
                    if responseMean < g.inhtestvalue/2
                        datasegment_inverse(results.offsetIdx(ii,1):end) = NaN;
                        results.offsetIdx(ii,1) = find(datasegment_inverse>=g.inhtestvalue,1,'last'); % end of response
                        responseMean = mean(datasegment_inverse(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
                    else 
                        k = 1;
                    end
                end
            else
                results.respDir(ii,1) = 0;
                results.onsetIdx(ii,1) = NaN;
                results.offsetIdx(ii,1) = NaN;
            end
        end
    end
    % heatmaps
    psth_spx_PN = psth_spx(idx_PN,:);
    psth_spx_IN = psth_spx(idx_IN,:);
    edges = 0:50:g.test_time*1000;
    onsetTime_d = discretize(results.onsetIdx*g.bin_time*1000, edges); % round to 10 ms values.
    offsetTime_d = discretize(results.offsetIdx*g.bin_time*1000, edges); % round to 10 ms values.
    durationTime_d = discretize((results.offsetIdx-results.onsetIdx)*g.bin_time*1000, edges); % round to 10 ms values.

    clear clusters psth_spx
    clusters(results.respDir == -1 & onsetTime_d == 1 & durationTime_d == 1,1) = 1;
    clusters(results.respDir == -1 & onsetTime_d == 1 & durationTime_d > 1,1) = 2;
    clusters(results.respDir == -1 & onsetTime_d > 1 ,1) = 2;
    clusters(results.respDir == 0) = 3;
    clusters(results.respDir == 1 & onsetTime_d == 1 & durationTime_d == 1,1) = 4;
    clusters(results.respDir == 1 & onsetTime_d == 1 & durationTime_d > 1,1) = 4;
    clusters(results.respDir == 1 & onsetTime_d > 1 ,1) = 4;
    clusters(results.respDir == 1 & isnan(onsetTime_d)) = 5;   

    [~, order_PN] = sortrows([results.respDir(idx_PN) clusters(idx_PN) results.onsetIdx(idx_PN)*g.bin_time*1000 results.offsetIdx(idx_PN)*g.bin_time*1000], [1 2 3 4], 'ascend');
    psth_spx = [psth_spx_PN(order_PN,:)];
    matrix = psth_spx(:,(g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    imagesc(g.timeaxis,1:size(matrix,1),matrix);
    hold on
    cludat = clusters(idx_PN);
    n_clu = find(diff(cludat(order_PN))~=0);
    yline(n_clu+0.5, 'Color', 'k', 'LineWidth', 1);
    hold off
    clim(g.clim)
    disp([-max(max(psth_spx,[],1))/2.5 max(max(psth_spx,[],1))])
    colormap(g.colors.Heatmap); 
    xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
    xticks([-0.2 0 0.2 0.4]);
    if hmp ==1 
        ylabel('Cell number')
    end
    if hmp == 4
        cb = colorbar('eastoutside', 'FontSize', g.fontSize2);
    end
    set(gca, 'FontSize', g.fontSize2);
    title(hmptitles{hmp}, 'FontSize', g.fontSize1)
    % plotting x symbols at onset and offset
    onset_PN = results.onsetIdx(idx_PN);
    offset_PN = results.offsetIdx(idx_PN);
    onsets_sorted = [onset_PN(order_PN)];
    offsets_sorted = [offset_PN(order_PN)];
    hold on
    % for ii = 1:size(matrix,1)
    %     plot(onsets_sorted(ii)*g.bin_time, ii, 'x', 'Color', 'r', 'MarkerSize', 4)
    %     plot(offsets_sorted(ii)*g.bin_time, ii, 'x', 'Color', 'b', 'MarkerSize', 4)
    % end

    hmpdata.(['hmp_' num2str(hmp)]) = clusters;
    hmpdata.(['order_' num2str(hmp)]).PN = order_PN;
end
clearvars -except g t hmpdata clusters


%% (3:5,1:4) - cluster spikes
ttl = {'triptest_sound_only', 'triptest_sound_only',  'triptest_shocks_only', 'triptest_shocks_only'};
br = {'LA', 'BA', 'LA', 'BA'};
hmptitles = {'LA population, CS reponse', 'LA population, US response', 'BA population, CS reponses', 'BA population, US response'};
layoutTile = [9 10 11 12];

for ii = 1:4
    ts = tiledlayout(t,3,1);
    ts.Layout.Tile = layoutTile(ii);
    ts.Layout.TileSpan = [3 1];
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{ii}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    time_axis = linspace(-g.plotwin(1),g.plotwin(2),(g.plotwin(1)+g.plotwin(2))/g.bin_time);
    idx_short = hmpdata.(['hmp_' num2str(ii)]) == 1;
    idx_long = hmpdata.(['hmp_' num2str(ii)]) == 2;
    idx_inh = hmpdata.(['hmp_' num2str(ii)]) == 4;
    idx_PN = strcmp(g.cell_metrics.brainRegion,br{ii}) & strcmp(g.cell_metrics.putativeCellType,'PN');
    idx_IN = strcmp(g.cell_metrics.brainRegion,br{ii}) & strcmp(g.cell_metrics.putativeCellType,'IN');   

    ax1 = nexttile(ts,[1 1]);
    idx =  idx_PN' & idx_short;
    spx = sum(psth_spx(idx,:),1);
    spx = smoothdata(sum(psth_spx(idx,:),1), 'sgolay', g.smoothvalue);
    spx = spx(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    b1 = bar(time_axis,spx, 'k'); % 'k' for black bars
    b1.FaceColor =  g.colors.barcolor;
    b1.EdgeColor = g.colors.barcolor; 
    ylim([0 150])
    yticks([0 75 150])
    xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
    xticks([-0.2 0 0.2 0.4]);
    if ii == 1
        ylabel('Number of spikes')
    end
    set(gca, 'FontSize', g.fontSize2);
    
    ax2 = nexttile(ts,[1 1]);
    idx =  idx_PN' & idx_long;
    spx = sum(psth_spx(idx,:),1);
    spx = smoothdata(sum(psth_spx(idx,:),1), 'sgolay', g.smoothvalue);
    spx = spx(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    b2 = bar(time_axis,spx, 'k'); % 'k' for black bars
    b2.FaceColor =  g.colors.barcolor;
    b2.EdgeColor = g.colors.barcolor; 
    ylim([0 500])
    yticks([0 250 500])
    xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
    xticks([-0.2 0 0.2 0.4]);
    if ii == 1
        ylabel('Number of spikes')
    end
    set(gca, 'FontSize', g.fontSize2);

    ax3 = nexttile(ts,[1 1]);
    idx =  idx_PN' & idx_inh;
    spx = sum(psth_spx(idx,:),1);
    spx = smoothdata(sum(psth_spx(idx,:),1), 'sgolay', g.smoothvalue);
    spx = spx(:, (g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    b3 = bar(time_axis,spx, 'k'); % 'k' for black bars
    b3.FaceColor =  g.colors.barcolor;
    b3.EdgeColor = g.colors.barcolor; 
    ylim([0 120])
    yticks([0 60 120])
    xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
    xticks([-0.2 0 0.2 0.4]);
    if ii == 1
        ylabel('Number of spikes')
    end
    xlabel('Time (s)')
    set(gca, 'FontSize', g.fontSize2);
end