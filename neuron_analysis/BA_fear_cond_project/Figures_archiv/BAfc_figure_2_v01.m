% function BAfc_figure_2
% % Plotting region specific responses on heatmaps
% % 
% recordings = {...
%     'MD292_002_kilosort',...
%     'MD293_kilosort',...
%     'MD294_kilosort',...
%     'MD295_kilosort',...
%     'MD296_kilosort',...
%     'MD297_kilosort'};
% g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'});

clearvars -except g
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 12;
g.fontSize2 = 12;
g.fontsize3 = 10;
g.colors = BAfc_colors;
g.pre_time = 4;
g.post_time = 4;
g.test_time = 1; % 1 sec is good, beacuse captures better than 0.5
g.exc_iflarger = 6; % if the mean change is negative, but there is a peak where zscore is above this, than counted as excited
g.bin_time = 0.01; % 0.01 volt sokat 5 os smoothal
g.smoothvalue = 5;
g.plotwin = [0.5 1];
g.spxwin = [0.02 0.1];
g.spxbin = 0.001;
g.timeaxis = -g.plotwin(1):g.bin_time:g.plotwin(2);
g.clim = [-4.1 10.2];

%% Initiate figure
fig = figure('Position', [400, 100, 1800, 1200]);
t = tiledlayout(fig,3,6,'TileSpacing', 'compact', 'Padding', 'none');
%% (1,[1 2 4 5]) - Heatmaps
mploc = [1 2 4 5];
ttl = {'triptest_sound_only', 'triptest_sound_only',  'triptest_shocks_only', 'triptest_shocks_only'};
br = {'LA', 'BA', 'LA', 'BA'};
hmptitles = {'LA population, CS reponse', 'BA population, CS response', 'LA population, US reponses', 'BA population, US response'};
for hmp = 1:4 
    ax = nexttile(t,mploc(hmp),[1 1]);
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
    psth_spx = zscore(psth_spx,0,2);   
    if any(g.smoothvalue)
        psth_spx = smoothdata(psth_spx,2,'movmean',g.smoothvalue);
    end
    [order, respData] = calc_conbinedSocre(psth_spx,g);
    idx_PN = strcmp(g.cell_metrics.brainRegion(order),br{hmp}) & strcmp(g.cell_metrics.putativeCellType(order),'PN');
    idx_IN = strcmp(g.cell_metrics.brainRegion(order),br{hmp}) & strcmp(g.cell_metrics.putativeCellType(order),'IN');
    order = [order(idx_PN); order(idx_IN)]; % locally selected neurosn
    psth_spx = psth_spx(order,:);
    matrix = psth_spx(:,(g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
    imagesc(g.timeaxis,1:size(matrix,1),matrix);
    title(hmptitles{hmp}, 'FontSize', g.fontSize2)
    clim(g.clim)
    disp([-max(max(psth_spx,[],1))/2.5 max(max(psth_spx,[],1))])
    colormap(g.colors.Heatmap); 
    idx_PN = strcmp(g.cell_metrics.brainRegion(order),br{hmp}) & strcmp(g.cell_metrics.putativeCellType(order),'PN');
    yline(sum(idx_PN)+0.5, 'k', 'LineWidth', 1);
    xlabel('Time (s)', 'FontSize', g.fontSize2)
    ylabel('Cell number', 'FontSize', g.fontSize2)
    % hold on
    % respDataOrd = respData(order,:);
    % for ii = 1:size(matrix,1)
    %     plot([respDataOrd(ii,1)*g.bin_time, respDataOrd(ii,2)*g.bin_time], [ii, ii], 'Color', 'k', 'LineWidth', 1)
    %     plot(respDataOrd(ii,3)*g.bin_time, ii, 'x', 'Color', 'r', 'MarkerSize', 4)
    % end
    % hold off    
    hmpData.(['hmp_' num2str(hmp)])(:,1) = num2cell(respData(order,4));
    hmpData.(['hmp_' num2str(hmp)])(:,2) = g.cell_metrics.putativeCellType(order)';
    if hmp == 1
        cb = colorbar('westoutside', 'FontSize',15);
    end

end

clearvars -except g t hmpData

%% (1, [3 6]) - Barplots
for ii = 1:4
    T = cell2table(hmpData.(['hmp_' num2str(ii)]), 'VariableNames', {'respDir', 'cellType'});
    PN_tot(1,ii) = sum(strcmp(T.cellType,'PN'));
    IN_tot(1,ii) = sum(strcmp(T.cellType,'IN'));
    PN_exc(1,ii) = sum(strcmp(T.cellType,'PN')&T.respDir==1);
    IN_exc(1,ii) = sum(strcmp(T.cellType,'IN')&T.respDir==1);
end

t2 = tiledlayout(t,2,1, 'TileSpacing', 'tigh');
t2.Layout.Tile = 3;
t2.Layout.TileSpan = [1 1];

ax = nexttile(t2, [1 1]);
b1 = bar(1, PN_exc(1)/PN_tot(1));  % LA PN exc
hold on
b2 = bar(2,  (PN_tot(1)-PN_exc(1))/PN_tot(1));  % LA PN inh
b3 = bar(4, PN_exc(2)/PN_tot(2));  % BA PN exc
b4 = bar(5,  (PN_tot(2)-PN_exc(2))/PN_tot(2));  % BA PN inh
b1.FaceColor = g.colors.c3;  
b2.FaceColor = g.colors.c1; 
b3.FaceColor = g.colors.c3;  
b4.FaceColor = g.colors.c1;  
ylim([0 1.5])
h = legend([b1, b2], {'Increased', 'Decreased'}, 'FontSize', g.fontsize3);
h.Location = 'northeast';    % or southwest, northwest, southeast
h.Box = 'off'; 
title(ax, 'Firing rate change, PN', 'FontSize', g.fontSize2)
xticks([1.5 4.5])
set(gca, 'XTickLabel', {'LA', 'BA'}, 'FontSize', g.fontSize2)

ax = nexttile(t2, [1 1]);
b1 = bar(1, IN_exc(1)/IN_tot(1));  % LA PN exc
hold on
b2 = bar(2,  (IN_tot(1)-IN_exc(1))/IN_tot(1));  % LA PN inh
b3 = bar(4, IN_exc(2)/IN_tot(2));  % BA PN exc
b4 = bar(5,  (IN_tot(2)-IN_exc(2))/IN_tot(2));  % BA PN inh
b1.FaceColor = g.colors.c3;  
b2.FaceColor = g.colors.c1; 
b3.FaceColor = g.colors.c3;  
b4.FaceColor = g.colors.c1;  
ylim([0 1.5])
h = legend([b1, b2], {'Increased', 'Decreased'}, 'FontSize', g.fontsize3);
h.Location = 'northeast';    % or southwest, northwest, southeast
h.Box = 'off'; 
title(ax, 'Firing rate change, IN', 'FontSize', g.fontSize2)
xticks([1.5 4.5])
set(gca, 'XTickLabel', {'LA', 'BA'}, 'FontSize', g.fontSize2)

t3 = tiledlayout(t,2,1, 'TileSpacing', 'tigh');
t3.Layout.Tile = 6;
t3.Layout.TileSpan = [1 1];

ax = nexttile(t3, [1 1]);
b1 = bar(1, PN_exc(3)/PN_tot(3));  % LA PN exc
hold on
b2 = bar(2,  (PN_tot(3)-PN_exc(3))/PN_tot(3));  % LA PN inh
b3 = bar(4, PN_exc(4)/PN_tot(4));  % BA PN exc
b4 = bar(5,  (PN_tot(4)-PN_exc(4))/PN_tot(4));  % BA PN inh
b1.FaceColor = g.colors.c3;  
b2.FaceColor = g.colors.c1; 
b3.FaceColor = g.colors.c3;  
b4.FaceColor = g.colors.c1;  
ylim([0 1.5])
h = legend([b1, b2], {'Increased', 'Decreased'}, 'FontSize', g.fontsize3);
h.Location = 'northeast';    % or southwest, northwest, southeast
h.Box = 'off'; 
title(ax, 'Firing rate change, PN', 'FontSize', g.fontSize2)
xticks([1.5 4.5])
set(gca, 'XTickLabel', {'LA', 'BA'}, 'FontSize', g.fontSize2)

ax = nexttile(t3, [1 1]);
b1 = bar(1, IN_exc(3)/IN_tot(3));  % LA PN exc
hold on
b2 = bar(2,  (IN_tot(3)-IN_exc(3))/IN_tot(3));  % LA PN inh
b3 = bar(4, IN_exc(4)/IN_tot(4));  % BA PN exc
b4 = bar(5,  (IN_tot(4)-IN_exc(4))/IN_tot(4));  % BA PN inh
b1.FaceColor = g.colors.c3;  
b2.FaceColor = g.colors.c1; 
b3.FaceColor = g.colors.c3;  
b4.FaceColor = g.colors.c1; 
ylim([0 1.5])
h = legend([b1, b2], {'Increased', 'Decreased'}, 'FontSize', g.fontsize3);
h.Location = 'northeast';    % or southwest, northwest, southeast
h.Box = 'off'; 
title(ax, 'Firing rate change, IN', 'FontSize', g.fontSize2)
xticks([1.5 4.5])
set(gca, 'XTickLabel', {'LA', 'BA'}, 'FontSize', g.fontSize2)

clearvars -except g t hmpData


%% (2,1:3) - Population spikes with LFP
ttl = {'triptest_sound_only', 'triptest_sound_only',  'triptest_shocks_only', 'triptest_shocks_only'};
br = {'LA', 'BA', 'LA', 'BA'};
hmptitles = {'LA population, CS reponse', 'LA population, US response', 'BA population, CS reponses', 'BA population, US response'};
layoutTile = [7 8 10 11];

for ii = 1:4
    ts = tiledlayout(t,2,1);
    ts.Layout.Tile = layoutTile(ii);
    ts.Layout.TileSpan = [1 1];
    
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{ii}, 'pre_time', g.spxwin(1), 'post_time', g.spxwin(2), 'bin_time', g.spxbin);
    time_axis = linspace(-g.spxwin(1),g.spxwin(2),size(psth_spx,2))*1000;
    
    ax1 = nexttile(ts,[1 1]);
    idx =  strcmp(g.cell_metrics.brainRegion,br{ii}) & strcmp(g.cell_metrics.putativeCellType,'PN');
    spx = sum(psth_spx(idx,:),1);
    %spx = smoothdata(sum(psth_spx(idx,:),1), 'movmean', 3);
    bar(time_axis,spx, 'k'); % 'k' for black bars
    ylim([0 50])
    xlabel('Time (ms)', 'FontSize', g.fontSize2)
    ylabel('Number of spikes', 'FontSize', g.fontSize2)
    
    ax2 = nexttile(ts,[1 1]);
    idx =  strcmp(g.cell_metrics.brainRegion,br{ii}) & strcmp(g.cell_metrics.putativeCellType,'IN');
    spx = sum(psth_spx(idx,:),1);
    %spx = smoothdata(sum(psth_spx(idx,:),1), 'movmean', 3);
    bar(time_axis,spx, 'k'); % 'k' for black bars
    ylim([0 150])
    xlabel('Time (ms)', 'FontSize', g.fontSize2)
    ylabel('Number of spikes', 'FontSize', g.fontSize2)
end
    








function [order, respData] = calc_conbinedSocre(psth_spx,g)
    roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
    psth_zscore = zeros(size(psth_spx));
    respData = zeros(size(psth_spx,1),3);
    means = zeros(size(psth_spx,1),1);
    for ii = 1:size(psth_spx,1)
        psth_zscore(ii,:) = (psth_spx(ii,:) - mean(psth_spx(ii,1:g.pre_time/g.bin_time)))/std(psth_spx(ii,1:g.pre_time/g.bin_time));
        % if the mean response is negative during test_time, and neved above g.exc_iflarger zscore, than finding negative peak
        if mean(psth_zscore(ii,roi)) < 0 && ~any(psth_zscore(ii,roi)>=g.exc_iflarger) 
            datasegment = -(psth_zscore(ii,roi));
        else
            datasegment = (psth_zscore(ii,roi));
        end
        if all(datasegment == datasegment(1))
            respData(ii,1) = 1;
            respData(ii,2) = length(datasegment);
            respData(ii,3) = 1;
        else
            resp_peak = max(datasegment);
            respData(ii,3) = find(datasegment==resp_peak,1,'first'); % peak of the response
            respData(ii,1) = find(datasegment>=resp_peak/2,1,'first'); % begining of the response
            datasegment(1:respData(ii,3)) = NaN;
            if isempty(find(datasegment<=resp_peak/2,1,'first'))
                respData(ii,2) = length(datasegment);
            else
                respData(ii,2) = find(datasegment<=resp_peak/2,1,'first'); % end of response
            end
        end
        means(ii) = mean(psth_spx(ii,g.pre_time/g.bin_time+respData(ii,1):g.pre_time/g.bin_time+respData(ii,2)))';
    end
    means = round(means*10)/10;
    respData(means<=0,4) = -1; % direction of response
    respData(means>0,4) = 1;  
    [~, order] = sortrows([respData(:,4) 1-respData(:,3)], [1 2], 'descend');
end
