cell_metrics = BAfc_load_neurons;
responses_sound =  BAfc_find_habit_v01';
responses_shock =  BAfc_find_shocks_v01';


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
    
    c_shock = categorical(shock_regions);
    cattype_shock = categories(c_shock);
    catnum_shock = countcats(c_shock);
    skip_shock = find(cattype_shock == "SKIP");
    cattype_shock(skip_shock) = [];
    catnum_shock(skip_shock) = [];
    [~,idx_shock] = sort(catnum_shock);
    
    colormap =  orderedcolors('reef');
    colormap = colormap(1:numel(catnum_shock), :);

    bar_shock = barh(catnum_shock(idx_shock));
    bar_shock.FaceColor = 'flat';
    for ii = 1:numel(catnum_shock)
        bar_shock.CData(ii,:) = colormap(idx_shock(ii),:);
    end
    ax = gca;
    ax.XLim = [0, 130];
    yticklabels(cattype_shock(idx_shock))
    title([shock_resp_types{aa} ' shock responses'])
    xlabel('Number of cells')
    ylabel('Brain region')
    
    ax = nexttile;
    piedata = [catnum_shock(idx_shock)];
    labels1 = cattype_shock(idx_shock);
    pie(ax, piedata, labels1);
    ax.Colormap = colormap(idx_shock,:);
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
    
    c_sound = categorical(sound_regions);
    cattype_sound = categories(c_sound);
    catnum_sound = countcats(c_sound);
    skip_sound = find(cattype_sound == "SKIP");
    cattype_sound(skip_sound) = [];
    catnum_sound(skip_sound) = [];
    [~,idx_sound] = sort(catnum_sound);
    
    colormap2 =  orderedcolors('reef');
    colormap2 = colormap2(1:numel(catnum_sound), :);

    bar_sound = barh(catnum_sound(idx_sound));
    bar_sound.FaceColor = 'flat';
    for ii = 1:numel(catnum_sound)
        bar_sound.CData(ii,:) = colormap2(idx_sound(ii),:);
    end
    ax = gca;
    ax.XLim = [0, 160];
    yticklabels(cattype_sound(idx_sound))
    title([sound_resp_types{aa} ' sound responses'])
    xlabel('Number of cells')
    ylabel('Brain region')
    
    ax = nexttile;
    piedata2 = [catnum_sound(idx_sound)];
    labels2 = cattype_sound(idx_sound);
    pie(ax, piedata2, labels2);
    ax.Colormap = colormap2(idx_sound,:);
end



% Plot US response of CS excited neurons
figure(3)
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


