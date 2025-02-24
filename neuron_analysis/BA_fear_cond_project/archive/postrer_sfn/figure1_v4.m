cell_metrics = BAfc_load_neurons;
[responses_CSplus_habit, FirstSpxEst_CSplus_habit, SmoothedZscore_habit, psth_spx_habit] =  BAfc_find_sound_v03('TTL_tone_habit_first');
[responses_CSplus_recall, FirstSpxEst_CSplus_recall, SmoothedZscore_recall, psth_spx_recall]  =  BAfc_find_sound_v03('TTL_tone_recall_first');
[responses_shock, FirstSpxEst_shock] =  BAfc_find_shocks_v03('TTL_shocks');
[pRanksum_plasticity, resRanksum_plasticity] = BAfc_find_plasticity; 

bR_order= {'Astria','CeA','BA_lat','BA_med','LA'};
bR_order_plotNames = {'Astria','CeA','BA Lat','BA Med','LA'};

for ii = 1:numel(bR_order)
    CellNum_order(ii) =  numel(find(contains(cell_metrics.brainRegion, bR_order{ii})));
end

colormap =  orderedcolors('reef');
colormap = colormap(1:numel(bR_order), :);

plotsubtitles = {'Excited neurons', 'Inhibited neurons', 'Non-responseive neurons'};
% Plot shock responses
shock_resp_types = categories(categorical(responses_shock));
figure(1)
tiledlayout(numel(shock_resp_types),1)
for aa = 1:length(shock_resp_types)
    idx = find(contains(responses_shock,shock_resp_types{aa}));
    nexttile
    shock_regions = cell_metrics.brainRegion(idx);
    shock_regions_uni = unique(shock_regions,'stable');
    shock_regions_num = cellfun(@(x) sum(ismember(shock_regions,x)),shock_regions_uni,'un',0);

    for bb = 1:numel(bR_order)
        if any(contains(shock_regions_uni, bR_order{bb}))
            respShockOrd.(shock_resp_types{aa})(bb) = cell2mat(shock_regions_num(find(contains(shock_regions_uni, bR_order(bb)))));
        else 
            respShockOrd.(shock_resp_types{aa})(bb) = 0;
        end
    end

    bar_shock = barh(respShockOrd.(shock_resp_types{aa}));
    bar_shock.FaceColor = 'flat';
    for ii = 1:numel(bR_order)
        bar_shock.CData(ii,:) = colormap((ii),:);
    end
    ax = gca;
    ax.XLim = [0, 130];
    yticklabels(bR_order_plotNames)
    title(plotsubtitles(aa))
    xlabel('Number of cells')
    ylabel('Brain region')

end
set(findall(gcf,'-property','FontSize'),'FontSize',20)

% Plot sound responses
sound_resp_types = categories(categorical(responses_CSplus_habit));
figure(2)
tiledlayout(numel(sound_resp_types),1)
for aa = 1:length(sound_resp_types)
    idx2 = find(contains(responses_CSplus_habit,sound_resp_types{aa}));
    nexttile
    sound_regions = cell_metrics.brainRegion(idx2);
    sound_regions_uni = unique(sound_regions,'stable');
    sound_regions_num = cellfun(@(x) sum(ismember(sound_regions,x)),sound_regions_uni,'un',0);

    for bb = 1:numel(bR_order)
        if any(contains(sound_regions_uni, bR_order{bb}))
            respSoundOrd.(sound_resp_types{aa})(bb) = cell2mat(sound_regions_num(find(contains(sound_regions_uni, bR_order(bb)))));
        else 
            respSoundOrd.(sound_resp_types{aa})(bb) = 0;
        end
    end    
   

    bar_sound = barh(respSoundOrd.(sound_resp_types{aa}));
    bar_sound.FaceColor = 'flat';
    for ii = 1:numel(bR_order)
        bar_sound.CData(ii,:) = colormap((ii),:);
    end
    ax = gca;
    ax.XLim = [0, 160];
    yticklabels(bR_order_plotNames)
    title((plotsubtitles(aa)))
    xlabel('Number of cells')
    ylabel('Brain region')
end
set(findall(gcf,'-property','FontSize'),'FontSize',20)

% Plot CS and US response percentages per region
% Must run fig(1) and (2) first
% Shock responses
plotsubtitles2 = {'Exc.', 'Inh.', 'Non-resp.'};
figure(3)
tiledlayout(1,2)
nexttile
for ee = 1:numel(shock_resp_types)
    respShockPct(ee,:) = respShockOrd.(shock_resp_types{ee})./CellNum_order*100;
end
barh_sh = barh(respShockPct', 'stacked', 'FaceColor','flat');
barh_sh(1).CData(1:5,:) = repelem([0.8500, 0.3250, 0.0980],5,1);
barh_sh(2).CData(1:5,:) = repelem([0, 0.4470, 0.7410],5,1);
barh_sh(3).CData(1:5,:) = repelem([0.85 0.85 0.85],5,1);
yticklabels(bR_order_plotNames)
title('Shock responses')
xlabel('Percentage of cells')
ylabel('Brain region')
Lgd = legend(plotsubtitles2, 'NumColumns',3, 'Location','south');
fontsize(Lgd,14,'points')
% Sound responses
nexttile
for ff = 1:numel(sound_resp_types)
    respSoundPct(ff,:) = respSoundOrd.(sound_resp_types{ff})./CellNum_order*100;
end
barh_so = barh(respSoundPct', 'stacked', 'FaceColor','flat');
barh_so = barh(respShockPct', 'stacked', 'FaceColor','flat');
barh_so(1).CData(1:5,:) = repelem([0.8500, 0.3250, 0.0980],5,1);
barh_so(2).CData(1:5,:) = repelem([0, 0.4470, 0.7410],5,1);
barh_so(3).CData(1:5,:) = repelem([0.85 0.85 0.85],5,1);
yticklabels(bR_order_plotNames)
title('Sound responses')
xlabel('Percentage of cells')
ylabel('Brain region')
Lgd = legend(plotsubtitles2, 'NumColumns',3, 'Location','south');
fontsize(Lgd,14,'points')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

% Find plasticity showing neurons
for ii = 1:numel(cell_metrics.cellID)
    if anynan(SmoothedZscore_habit(ii,10001:14500)) || anynan(SmoothedZscore_recall(ii,10001:14500))
        pRanksum_plasticity(ii) = 0;
        resRanksum_plasticity(ii) = 0;
    else
        [pRanksum_plasticity(ii), resRanksum_plasticity(ii)] = signrank(SmoothedZscore_habit(ii,10001:14500), SmoothedZscore_recall(ii,10001:14500), 'alpha', 0.05);
    end
end

soundSign = unique([...
    find(contains(responses_CSplus_habit, 'exc')) ...
    find(contains(responses_CSplus_habit, 'inh')) ...
    find(contains(responses_CSplus_recall, 'exc')) ...
    find(contains(responses_CSplus_recall, 'inh'))]);

soundSignPlast = intersect(soundSign, find(resRanksum_plasticity == 1));

cell_metrics.labels(1:end) = {'neutral'};
cell_metrics.labels(soundSignPlast) = {'plast'};
cell_metrics = CellExplorer('metrics',cell_metrics);








clearvars -except cell_metrics responses_CSplus_habit responses_CSplus_recall responses_shock bR_order bR_order_plotNames colormap CellNum_order...
    FirstSpxEst_CSplus_habit FirstSpxEst_CSplus_recall FirstSpxEst_shock 
