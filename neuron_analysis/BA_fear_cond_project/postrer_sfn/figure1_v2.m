cell_metrics = BAfc_load_neurons;
responses_sound =  BAfc_find_habit_v01';
responses_shock =  BAfc_find_shocks_v01';

bR_order= flip({'LA', 'BA_med', 'BA_lat', 'CeA', 'Astria'});
colormap =  orderedcolors('reef');
colormap = colormap(1:numel(bR_order), :);

% Plot shock responses
shock_resp_types = categories(categorical(responses_shock));
figure(1)
tiledlayout(numel(shock_resp_types),2)
for aa = 1:length(shock_resp_types)
    idx = find(contains(responses_shock,shock_resp_types{aa}));
    nexttile
    shock_regions = cell_metrics.brainRegion(idx);
    shock_regions_uni = unique(shock_regions,'stable');
    shock_regions_num = cellfun(@(x) sum(ismember(shock_regions,x)),shock_regions_uni,'un',0);

    for bb = 1:numel(bR_order)
        if any(contains(shock_regions_uni, bR_order{bb}))
            respShockOrd(bb) = cell2mat(shock_regions_num(find(contains(shock_regions_uni, bR_order(bb)))));
        else 
            respShockOrd(bb) = 0;
        end
    end

    bar_shock = barh(respShockOrd);
    bar_shock.FaceColor = 'flat';
    for ii = 1:numel(bR_order)
        bar_shock.CData(ii,:) = colormap((ii),:);
    end
    ax = gca;
    ax.XLim = [0, 130];
    yticklabels(bR_order)
    title([shock_resp_types{aa} ' shock responses'])
    xlabel('Number of cells')
    ylabel('Brain region')

    ax = nexttile;
    piedata1 = respShockOrd;
    labels1 = bR_order;
    pie(ax, piedata1, labels1);
    ax.Colormap = colormap;
end

% Plot sound responses
sound_resp_types = categories(categorical(responses_sound));
figure(2)
tiledlayout(numel(sound_resp_types),2)
for aa = 1:length(sound_resp_types)
    idx2 = find(contains(responses_sound,sound_resp_types{aa}));
    nexttile
    sound_regions = cell_metrics.brainRegion(idx2);
    sound_regions_uni = unique(sound_regions,'stable');
    sound_regions_num = cellfun(@(x) sum(ismember(sound_regions,x)),sound_regions_uni,'un',0);

    for bb = 1:numel(bR_order)
        if any(contains(sound_regions_uni, bR_order{bb}))
            respSoundOrd(bb) = cell2mat(sound_regions_num(find(contains(sound_regions_uni, bR_order(bb)))));
        else 
            respSoundOrd(bb) = 0;
        end
    end    
   

    bar_sound = barh(respSoundOrd);
    bar_sound.FaceColor = 'flat';
    for ii = 1:numel(bR_order)
        bar_sound.CData(ii,:) = colormap((ii),:);
    end
    ax = gca;
    ax.XLim = [0, 160];
    yticklabels(bR_order)
    title([sound_resp_types{aa} ' sound responses'])
    xlabel('Number of cells')
    ylabel('Brain region')
    
    ax = nexttile;
    piedata2 = respSoundOrd;
    labels2 = bR_order;
    pie(ax, piedata2, labels2);
    ax.Colormap = colormap;
end



% Plot US response of CS excited neurons
figure(4)
tiledlayout(1,3)
cs_resps = {'exc', 'inh', 'neutral'};
for bb = 1:numel(cs_resps)
    CS_exc =  find(contains(responses_sound, cs_resps(bb)));
    CS_exc_US_resp = responses_shock(CS_exc);
    
    c_CS_exc_US = categorical(CS_exc_US_resp);
    cattype_c_CS_exc_US  = categories(c_CS_exc_US);
    catnum_c_CS_exc_US  = countcats(c_CS_exc_US);
       
    ax = nexttile;
    piedata3 = catnum_c_CS_exc_US;
    for cc = 1:numel(cattype_c_CS_exc_US)
        labels3{cc} = [cattype_c_CS_exc_US{cc} ' = ' num2str(catnum_c_CS_exc_US(cc))];
    end
    pie(ax, piedata3, labels3);
    title(['US responses of CS ' cs_resps{bb} ' neurons'])

end



clearvars -except cell_metrics responses_sound responses_shock bR_order colormap