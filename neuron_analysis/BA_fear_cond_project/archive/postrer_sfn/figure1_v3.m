cell_metrics = BAfc_load_neurons;
[responses_CSplus_habit, FirstSpxEst_CSplus_habit] =  BAfc_find_sound_v03('TTL_tone_habit_first');
[responses_CSplus_recall, FirstSpxEst_CSplus_recall]  =  BAfc_find_sound_v03('TTL_tone_recall_first');
[responses_shock, FirstSpxEst_shock] =  BAfc_find_shocks_v03('TTL_shocks');

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

    % ax = nexttile;
    % piedata1 = respShockOrd.(shock_resp_types{aa});
    % labels1 = bR_order_plotNames;
    % pie(ax, piedata1, labels1);
    % ax.Colormap = colormap;
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
    % 
    % ax = nexttile;
    % piedata2 = respSoundOrd.(sound_resp_types{aa});
    % labels2 = bR_order_plotNames;
    % pie(ax, piedata2, labels2);
    % ax.Colormap = colormap;
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
barh_so = barh(respSoundPct', 'stacked', 'FaceColor','flat');
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



% % Plot US response of CS excited neurons
% figure(4)
% tiledlayout(1,3)
% cs_resps = {'exc', 'inh', 'neutral'};
% for bb = 1:numel(cs_resps)
%     CS_exc =  find(contains(responses_CSplus_habit, cs_resps(bb)));
%     CS_exc_US_resp = responses_shock(CS_exc);
% 
%     c_CS_exc_US = categorical(CS_exc_US_resp);
%     cattype_c_CS_exc_US  = categories(c_CS_exc_US);
%     catnum_c_CS_exc_US  = countcats(c_CS_exc_US);
% 
%     ax = nexttile;
%     piedata3 = catnum_c_CS_exc_US;
%     for cc = 1:numel(cattype_c_CS_exc_US)
%         labels3{cc} = [cattype_c_CS_exc_US{cc} ' = ' num2str(catnum_c_CS_exc_US(cc))];
%     end
%     pie(ax, piedata3, labels3);
%     title(['US responses of CS ' cs_resps{bb} ' neurons'])
% 
% end
% 
% figure(5)
% tiledlayout(1,9)
% resps = {'exc', 'inh', 'neutral'};
% iter = 1;
% for ff = 1:3
    for gg = 1:3
        labels4 = [];
        idx_recall{iter} = intersect(find(contains(responses_CSplus_habit,resps{ff})), find(contains(responses_shock, resps{gg})));
        cat_recall{iter} = categorical(responses_CSplus_recall(idx_recall{iter}));

        catType_rec = categories(cat_recall{iter});
        catNum_rec = countcats(cat_recall{iter});
        ax = nexttile;
        piedata4 = catNum_rec;
        for hh = 1:numel(catType_rec)
            labels4{hh} = [catType_rec{hh} ' = ' num2str(catNum_rec(hh))];
        end
        pie(ax, piedata4, labels4);
        iter = iter+1;
    end
end




% Find plasticity showing neurons
figure(6)
respDiff = ~strcmp(responses_CSplus_habit,responses_CSplus_recall);
idx_respDiff = find(respDiff == 1);
plastic_regions = cell_metrics.brainRegion(idx_respDiff);
plastic_regions_uni = unique(plastic_regions,'stable');
plastic_regions_num = cellfun(@(x) sum(ismember(plastic_regions,x)),plastic_regions_uni,'un',0);

for ll = 1:numel(bR_order)
    if any(contains(plastic_regions_uni, bR_order{ll}))
        plasticOrd(ll) = cell2mat(plastic_regions_num(find(contains(plastic_regions_uni, bR_order(ll)))));
    else 
        plasticOrd(ll) = 0;
    end
end    
bar_plastic = barh(plasticOrd);
bar_plastic.FaceColor = 'flat';
for ii = 1:numel(bR_order)
    bar_plastic.CData(ii,:) = colormap((ii),:);
end
ax = gca;
ax.XLim = [0, 80];
yticklabels(bR_order_plotNames)
title('Plasticity showing neurons')
xlabel('Number of cells')
ylabel('Brain region')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

% plotting percentages of plasticity per region
figure(7)
plastic_change_types = {'increased', 'decreased', 'unchanged'};
for mm = 1:numel(bR_order)
    plasticity.idx.(bR_order{mm}) = intersect(idx_respDiff, find(contains(cell_metrics.brainRegion, bR_order{mm}))); % indices of plastic neurons per region
    plasticity.respHabit.(bR_order{mm}) = responses_CSplus_habit(plasticity.idx.(bR_order{mm}));
    plasticity.respRecall.(bR_order{mm}) = responses_CSplus_recall(plasticity.idx.(bR_order{mm}));
end
resps = {'exc', 'inh', 'neutral'};
for nn = 1:numel(bR_order)
    for oo = 1:numel(plasticity.respHabit.(bR_order{nn}))
        if strcmp(plasticity.respHabit.(bR_order{nn})(oo), resps{1}) & strcmp(plasticity.respRecall.(bR_order{nn})(oo), resps{1})
            plasticity.respChange.(bR_order{nn})(oo) = {'unchanged'};
        elseif strcmp(plasticity.respHabit.(bR_order{nn})(oo), resps{1}) & strcmp(plasticity.respRecall.(bR_order{nn})(oo), resps{2})
            plasticity.respChange.(bR_order{nn})(oo) = {'decreased'};
        elseif strcmp(plasticity.respHabit.(bR_order{nn})(oo), resps{1}) & strcmp(plasticity.respRecall.(bR_order{nn})(oo), resps{3})
            plasticity.respChange.(bR_order{nn})(oo) = {'decreased'};
        elseif strcmp(plasticity.respHabit.(bR_order{nn})(oo), resps{2}) & strcmp(plasticity.respRecall.(bR_order{nn})(oo), resps{1})
            plasticity.respChange.(bR_order{nn})(oo) = {'increased'};
        elseif strcmp(plasticity.respHabit.(bR_order{nn})(oo), resps{2}) & strcmp(plasticity.respRecall.(bR_order{nn})(oo), resps{2})
            plasticity.respChange.(bR_order{nn})(oo) = {'unchanged'};
        elseif strcmp(plasticity.respHabit.(bR_order{nn})(oo), resps{2}) & strcmp(plasticity.respRecall.(bR_order{nn})(oo), resps{3})
            plasticity.respChange.(bR_order{nn})(oo) = {'increased'};
        elseif strcmp(plasticity.respHabit.(bR_order{nn})(oo), resps{3}) & strcmp(plasticity.respRecall.(bR_order{nn})(oo), resps{1})
            plasticity.respChange.(bR_order{nn})(oo) = {'increased'};
        elseif strcmp(plasticity.respHabit.(bR_order{nn})(oo), resps{3}) & strcmp(plasticity.respRecall.(bR_order{nn})(oo), resps{2})
            plasticity.respChange.(bR_order{nn})(oo) = {'decreased'};
        elseif strcmp(plasticity.respHabit.(bR_order{nn})(oo), resps{3}) & strcmp(plasticity.respRecall.(bR_order{nn})(oo), resps{3})
            plasticity.respChange.(bR_order{nn})(oo) = {'unchanged'};
        end
    end
end
for pp = 1:numel(bR_order)
    for qq = 1:numel(plastic_change_types)
        plasticity.bardataPct(qq,pp) = numel(find(contains(plasticity.respChange.(bR_order{pp}), plastic_change_types{qq})))...
            /numel(plasticity.respChange.(bR_order{pp}))*100;
    end
end
bar_plasticPct = barh(plasticity.bardataPct', 'stacked', 'FaceColor','flat');
bar_plasticPct(1).CData(1:5,:) = repelem([0.8500, 0.3250, 0.0980],5,1);
bar_plasticPct(2).CData(1:5,:) = repelem([0, 0.4470, 0.7410],5,1);
bar_plasticPct(3).CData(1:5,:) = repelem([0.85 0.85 0.85],5,1);


yticklabels(bR_order_plotNames)
title('Plasticity showing neurons')
xlabel('Percentage of cells')
ylabel('Brain region')
Lgd = legend(plastic_change_types(1:2), 'NumColumns',2, 'Location','south');
fontsize(Lgd,14,'points')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

% Plotting dynamic changes of short latency US responses
blockStarts = [1 5 9 12 16];
idx_fastRespBlock1 = find(contains(responses_shock_class, 'excEarly'))';
idx_fastRespBlock5 = find(contains(responses_shock_class_block5, 'excEarly'))';
idx_fastRespBoth = intersect(idx_fastRespBlock1, idx_fastRespBlock5);
idx_fastRespLoss = setdiff(idx_fastRespBlock1, idx_fastRespBlock5);
idx_fastRespGain = setdiff(idx_fastRespBlock5, idx_fastRespBlock1);
[~, response_data_fastResp] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_shocks'},...
   'TTLselect', [1:20], 'psth_bin', 3000, 'int', [-0 0.1]); % each ttl is a column.
response_data_fastResp = cat(1,response_data_fastResp{:});
for ss = 1:numel(blockStarts) % number of blocks   
    response_data_fastResp_blocks(:,ss) = sum(response_data_fastResp(:, blockStarts(ss):blockStarts(ss)+3),2);
end

response_data_fastResp_blocks_norm = response_data_fastResp_blocks./response_data_fastResp_blocks(:,1);
figure(8)
plot(response_data_fastResp_blocks_norm([idx_fastRespBlock1;idx_fastRespBlock5],:)')
ylim([0 10])


% Plotting dynamic changes of long latency US responses
blockStarts = [1 5 9 12 16];
idx_longRespBlock1 = find(contains(responses_shock_class, 'excLate'))';
idx_longRespBlock5 = find(contains(responses_shock_class_block5, 'excLate'))';
idx_lateRespBoth = intersect(idx_longRespBlock1, idx_longRespBlock5);
idx_lateRespLoss = setdiff(idx_longRespBlock1, idx_longRespBlock5);
idx_lateRespGain = setdiff(idx_longRespBlock5, idx_longRespBlock1);
[~, response_data_longResp] = BAfc_calc_Zscore_prenorm('Stims', {'TTL_shocks'},...
   'TTLselect', [1:20], 'psth_bin', 60000, 'int', [-0 2]); % each ttl is a column.
response_data_longResp = cat(1,response_data_longResp{:});
for ss = 1:numel(blockStarts) % number of blocks   
    response_data_longResp_blocks(:,ss) = sum(response_data_longResp(:, blockStarts(ss):blockStarts(ss)+3),2);
end

response_data_longResp_blocks_norm = response_data_longResp_blocks./response_data_longResp_blocks(:,1);
figure(9)
plot(response_data_longResp_blocks_norm([idx_longRespBlock1;idx_longRespBlock5],:)')
ylim([0 10])


plot(response_data_longResp_blocks_norm(1,:)')

clearvars -except cell_metrics responses_CSplus_habit responses_CSplus_recall responses_shock bR_order colormap CellNum_order
PSTHall = (zscore(response_data_longResp_blocks_norm([idx_longRespBlock1;idx_longRespBlock5],:),[],2));