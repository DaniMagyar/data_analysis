% function BAfc_figure_5

% ha a NEM PN-eket akarom a heatmapre plottolni, akkor csak az strcmp('PN')
% ele kell egy ~ jel.


clear all

recordings = {...
    'MD307_kilosort',...
    'MD309_kilosort'};

g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', {'triptest_sound_only', 'triptest_sound_only_light', 'triptest_shocks_only', 'triptest_shocks_only_light', 'triptest_both', 'triptest_both_light'});

clearvars -except g
g.cell_metrics = BAfc_putative_cellTypes('cell_metrics', g.cell_metrics);
g.mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
g.colors = BAfc_colors;
g.fontSize1 = 13;
g.fontSize2 = 12;
g.pre_time = 5;
g.post_time = 0.5;
g.bin_time = 0.02; % 0.01 volt sokat 5 os smoothal
g.timeaxis = -g.pre_time:g.bin_time:g.post_time;
g.clustnum = 4;
g.smoothvalue = 11;
g.test_time = 0.5; % 1 sec is good, beacuse captures better than 0.5
g.plotwin = [g.pre_time g.post_time];
g.clim = [-2 5];
g.xlinewidth = 2;
g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time;
g.roi_pca = round((g.pre_time-0.5)/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);


%% Initiate figure
fig = figure('Position', [400, 100, 1800, 1200]);
t = tiledlayout(fig,16,12,'TileSpacing', 'tight', 'Padding', 'none');

%% (1,1) - Schematic figure
ax1 = nexttile(t,1,[4 4]);
[img, cmap] = imread([g.mainFolder '\drawed_injection.png']);
if ~isempty(cmap)
    img = ind2rgb(img, cmap);  % Convert to RGB
end
imshow(img, 'Parent', ax1);
title(ax1, 'Experimental setup', 'FontSize', g.fontSize1)
clearvars -except t g % clear variables


%% Example rasters

cellID_IN = intersect(find(strcmp(g.cell_metrics.animal, 'MD307')), find(g.cell_metrics.cluID == 18)); 
cellID_PN = intersect(find(strcmp(g.cell_metrics.animal, 'MD307')), find(g.cell_metrics.cluID == 163)); 

stimulus_times{1} = g.cell_metrics.general.triptest_sound_only{cellID_IN};
stimulus_times{2} = g.cell_metrics.general.triptest_sound_only{cellID_PN};
stimulus_times{3} = g.cell_metrics.general.triptest_sound_only_light{cellID_IN};
stimulus_times{4} = g.cell_metrics.general.triptest_sound_only_light{cellID_PN};
spike_times{1} = g.cell_metrics.spikes.times{cellID_IN}; 
spike_times{2} = g.cell_metrics.spikes.times{cellID_PN}; 
spike_times{3} = g.cell_metrics.spikes.times{cellID_IN}; 
spike_times{4} = g.cell_metrics.spikes.times{cellID_PN}; 

% Plot raster
time_bins = -g.pre_time:g.bin_time:g.post_time;    
spike_counts = zeros(1, length(time_bins) - 1);   
rasterloc = [5 29 7 31];
for ii = 1:4
    ax = nexttile(t,rasterloc(ii),[2 2]);
    hold on;    
    for trial = 1:length(stimulus_times{ii})
        aligned_spikes = spike_times{ii} - stimulus_times{ii}(trial);
        valid_spikes = aligned_spikes(aligned_spikes >= -g.pre_time & aligned_spikes <= g.post_time);
        scatter(valid_spikes, trial * ones(size(valid_spikes)), 10, 'k', 'filled');
        spike_counts = spike_counts + histcounts(valid_spikes, time_bins);
    end
end

%% inhibition effect
% plot baseline firing rate change of photoinhibited neurons
ttl_all_inh = {'triptest_sound_only_light', 'triptest_shocks_only_light', 'triptest_both_light'};
ttl_all = {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'};
pretime_fr = 1;
bintime_fr = 0.005;
smoothvalue_fr = 5;
inhibvalue = 0.5; % zscore. 
excitvalue = 3;
for ii = 1:3
    psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl_all_inh{ii}, 'pre_time', pretime_fr, 'post_time', 0, 'bin_time', bintime_fr);
    if exist('psth_spx_all', 'var')
        psth_spx_all = psth_spx_all + psth_spx;
    else
        psth_spx_all = psth_spx;
    end
    num_ttl(:,ii) = cellfun(@numel,g.cell_metrics.general.(ttl_all_inh{ii}))';
end
fr_baseline = sum(psth_spx_all(:,1:(pretime_fr-0.5)/bintime_fr),2)./(sum(num_ttl,2)*pretime_fr);
fr_inhib = sum(psth_spx_all(:,(pretime_fr-0.5)/bintime_fr+1:end),2)./(sum(num_ttl,2)*0.5);
psth_spx_all_zscore = zscore(psth_spx_all,0,2);
psth_spx_all_zscore = smoothdata(psth_spx_all_zscore,2,'sgolay',smoothvalue_fr);
n_inhibbins = sum(psth_spx_all_zscore(:,(pretime_fr-0.5)/bintime_fr+1:end)<-inhibvalue,2);
idx_inhibited = n_inhibbins > (0.5/bintime_fr)/2; % 0.5 is the inhibition time, it's divided by bintime. If half of the bins are inhibited, cell inhibited

ax = nexttile(t,9,[4 1]);
jitter_amount = 0.1;
jitter_offset = (rand(sum(idx_inhibited), 1) - 0.5) * jitter_amount;
scatter(1 + jitter_offset, fr_baseline(idx_inhibited))
hold on
scatter(2 + jitter_offset, fr_inhib(idx_inhibited))
x_coords = [1 + jitter_offset, 2 + jitter_offset]';
y_coords = [fr_baseline(idx_inhibited), fr_inhib(idx_inhibited)]';
plot(x_coords, y_coords, 'k-')

% plot response magnitude change of photoinhibited neurons
posttime_fr = 0.5;
psth_spx_sound = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'triptest_sound_only', 'pre_time', pretime_fr, 'post_time', posttime_fr, 'bin_time', bintime_fr);
psth_spx_sound_light = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'triptest_sound_only_light', 'pre_time', pretime_fr, 'post_time', posttime_fr, 'bin_time', bintime_fr);
psth_spx_shock = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'triptest_shocks_only', 'pre_time', pretime_fr, 'post_time', posttime_fr, 'bin_time', bintime_fr);
psth_spx_shock_light = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'triptest_shocks_only_light', 'pre_time', pretime_fr, 'post_time', posttime_fr, 'bin_time', bintime_fr);

psth_spx_sound_zscore = (psth_spx_sound-mean(psth_spx_sound(:,1:(pretime_fr-0.5)/bintime_fr),2))./std(psth_spx_sound(:,1:(pretime_fr-0.5)/bintime_fr),0,2);
psth_spx_sound_light_zscore = (psth_spx_sound_light-mean(psth_spx_sound_light(:,1:(pretime_fr-0.5)/bintime_fr),2))./std(psth_spx_sound_light(:,1:(pretime_fr-0.5)/bintime_fr),0,2);
psth_spx_shock_zscore = (psth_spx_shock-mean(psth_spx_shock(:,1:(pretime_fr-0.5)/bintime_fr),2))./std(psth_spx_shock(:,1:(pretime_fr-0.5)/bintime_fr),0,2);
psth_spx_shock_light_zscore = (psth_spx_shock_light-mean(psth_spx_shock_light(:,1:(pretime_fr-0.5)/bintime_fr),2))./std(psth_spx_shock_light(:,1:(pretime_fr-0.5)/bintime_fr),0,2);

psth_spx_sound_zscore = smoothdata(psth_spx_sound_zscore,2,'sgolay', smoothvalue_fr);
psth_spx_sound_light_zscore = smoothdata(psth_spx_sound_light_zscore,2,'sgolay', smoothvalue_fr);
psth_spx_shock_zscore = smoothdata(psth_spx_shock_zscore,2,'sgolay', smoothvalue_fr);
psth_spx_shock_light_zscore = smoothdata(psth_spx_shock_light_zscore,2,'sgolay', smoothvalue_fr);

mean_resp_sound = mean(psth_spx_sound_zscore(:,pretime_fr/bintime_fr+1:end),2);
mean_resp_sound_light = mean(psth_spx_sound_light_zscore(:,pretime_fr/bintime_fr+1:end),2);
mean_resp_shock = mean(psth_spx_shock_zscore(:,pretime_fr/bintime_fr+1:end),2);
mean_resp_shock_light = mean(psth_spx_shock_light_zscore(:,pretime_fr/bintime_fr+1:end),2);

ax = nexttile(t,10,[4 1]);
jitter_amount = 0.1;
jitter_offset = (rand(sum(idx_inhibited), 1) - 0.5) * jitter_amount;
scatter(1 + jitter_offset, mean_resp_sound(idx_inhibited))
hold on
scatter(2 + jitter_offset, mean_resp_sound_light(idx_inhibited))
x_coords = [1 + jitter_offset, 2 + jitter_offset]';
y_coords = [mean_resp_sound(idx_inhibited), mean_resp_sound_light(idx_inhibited)]';
plot(x_coords, y_coords, 'k-')
%ylim([-1 4])

ax = nexttile(t,11,[4 1]);
jitter_amount = 0.1;
jitter_offset = (rand(sum(idx_inhibited), 1) - 0.5) * jitter_amount;
scatter(1 + jitter_offset, mean_resp_shock(idx_inhibited))
hold on
scatter(2 + jitter_offset, mean_resp_shock_light(idx_inhibited))
x_coords = [1 + jitter_offset, 2 + jitter_offset]';
y_coords = [mean_resp_shock(idx_inhibited), mean_resp_shock_light(idx_inhibited)]';
plot(x_coords, y_coords, 'k-')
%ylim([-1 4])

n_excitbins_sound = sum(psth_spx_sound_zscore(:,pretime_fr/bintime_fr+1:end)>excitvalue,2);
n_excitbins_shock = sum(psth_spx_shock_zscore(:,pretime_fr/bintime_fr+1:end)>excitvalue,2);

% idx_excited_sound = n_excitbins_sound > (posttime_fr/bintime_fr)/10; % at least half of the bins (/2) are above excitvalue
% idx_excited_shock = n_excitbins_shock > (posttime_fr/bintime_fr)/10; 
idx_excited_sound = n_excitbins_sound > 0; % at least half of the bins (/2) are above excitvalue
idx_excited_shock = n_excitbins_shock > 0; 



% trying to compare cell responses individually
% sound responses 

[~, preAP_sound, postAP_sound] = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'triptest_sound_only', 'pre_time', 0.5, 'post_time', 0.5, 'bin_time', bintime_fr);
[~, preAP_sound_light, postAP_sound_light] = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'triptest_sound_only_light', 'pre_time', 0.5, 'post_time', 0.5, 'bin_time', bintime_fr);
[~, preAP_sound_light_shifted, postAP_sound_light_shifted] = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'triptest_sound_only_light', 'pre_time', 0.5, 'post_time', 0.5, 'bin_time', bintime_fr, 'TTLshift', -0.5);

for ii = 1:size(preAP_sound,2)
    [~, h_preAP_sound(ii,1)] = signrank(preAP_sound{ii}, postAP_sound{ii});
    [~, h_preAP_sound_light(ii,1)] = signrank(preAP_sound_light{ii}, postAP_sound_light{ii});
    [~, h_inhibited_sound(ii,1)] = signrank(preAP_sound_light_shifted{ii}, postAP_sound_light_shifted{ii});
    [~, h_postAP_diff_sound(ii,1)] = signrank(postAP_sound{ii}, postAP_sound_light{ii});
end

sound_exc_exc_increased = h_preAP_sound & h_preAP_sound_light & h_postAP_diff_sound & ~idx_inhibited;





% shock responses
[~, preAP_shock, postAP_shock] = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'triptest_shocks_only', 'pre_time', 0.5, 'post_time', 0.5, 'bin_time', bintime_fr);
[~, preAP_shock_light, postAP_shock_light] = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'triptest_shocks_only_light', 'pre_time', 0.5, 'post_time', 0.5, 'bin_time', bintime_fr);
[~, preAP_shock_light_shifted, postAP_shock_light_shifted] = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'triptest_shocks_only_light', 'pre_time', 0.5, 'post_time', 0.5, 'bin_time', bintime_fr, 'TTLshift', -0.5);

for ii = 1:size(preAP_shock,2)
    [~, h_preAP_shock(ii,1)] = signrank(preAP_shock{ii}, postAP_shock{ii});
    [~, h_preAP_shock_light(ii,1)] = signrank(preAP_shock_light{ii}, postAP_shock_light{ii});
    [~, h_inhibited_shock(ii,1)] = signrank(preAP_shock_light_shifted{ii}, postAP_shock_light_shifted{ii});
    [~, h_postAP_diff_shock(ii,1)] = signrank(postAP_shock{ii}, postAP_shock_light{ii});
end

shock_exc_exc_increased = h_preAP_shock & h_preAP_shock_light & h_postAP_diff_shock & ~idx_inhibited;




% % idx_excited_sound = sound_exc_exc_increased;
% idx_excited_shock = shock_exc_exc_increased;

% n_spx_resp_sound = sum(psth_spx_sound(:,pretime_fr/bintime_fr+1:end),2);
% n_spx_resp_sound_light = sum(psth_spx_sound_light(:,pretime_fr/bintime_fr+1:end),2);
% n_spx_resp_shock = sum(psth_spx_shock(:,pretime_fr/bintime_fr+1:end),2);
% n_spx_resp_shock_light = sum(psth_spx_shock_light(:,pretime_fr/bintime_fr+1:end),2);
% 
% 

% baseline spikes
for ii = 1:3
    [~, preAP_light, postAP_light] =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl_all_inh{ii}, 'pre_time', 0.5, 'post_time', 0.5, 'bin_time', bintime_fr, 'TTLshift', -0.5);
    [~, preAP, postAP] =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl_all{ii}, 'pre_time', 0.5, 'post_time', 0.5, 'bin_time', bintime_fr, 'TTLshift', -0.5);
    if exist('preAP_light_all', 'var')
        preAP_light_all = cellfun(@plus, preAP_light_all, preAP_light, 'UniformOutput', false);
        postAP_light_all = cellfun(@plus, postAP_light_all, postAP_light, 'UniformOutput', false);
        preAP_all = cellfun(@plus, preAP_all, preAP, 'UniformOutput', false);
        postAP_all = cellfun(@plus, postAP_all, postAP, 'UniformOutput', false);
    else
        preAP_light_all = preAP_light;
        postAP_light_all = postAP_light;
        preAP_all = preAP;
        postAP_all = postAP;
    end
end





idx_manual_sound = [1, 10, 17, 19, 23, 25, 35, 36, 56];
idx_excited_sound = logical(zeros(size(preAP_sound)))';
idx_excited_sound(idx_manual_sound) = 1;

idx_manual_shock = [13, 35, 36, 44, 51, 56, 57];
idx_excited_shock = logical(zeros(size(preAP_shock)))';
idx_excited_shock(idx_manual_shock) = 1;



for ii = 1:size(preAP_sound,2)
    diff_preAP_sound(ii) = mean(postAP_sound{ii}-preAP_sound{ii});
    diff_preAP_sound_light(ii) = mean(postAP_sound_light{ii}-preAP_sound_light{ii});
    diff_preAP_shock(ii) = mean(postAP_shock{ii}-preAP_shock{ii});
    diff_preAP_shock_light(ii) = mean(postAP_shock_light{ii}-preAP_shock_light{ii});
    diff_preAP_light_all(ii,1) = mean(postAP_light_all{ii}-preAP_light_all{ii});
    diff_preAP_all(ii,1) = mean(postAP_all{ii}-preAP_all{ii});
end

n_spx_resp_sound = diff_preAP_sound';
n_spx_resp_sound_light = diff_preAP_sound_light';
n_spx_resp_shock = diff_preAP_shock';
n_spx_resp_shock_light = diff_preAP_shock_light';








% 
% ax = nexttile(t,49,[4 1]);
% jitter_amount = 0.1;
% jitter_offset = (rand(sum(idx_excited_sound), 1) - 0.5) * jitter_amount;
% scatter(1 + jitter_offset, mean_resp_sound(idx_excited_sound))
% hold on
% scatter(2 + jitter_offset, mean_resp_sound_light(idx_excited_sound))
% x_coords = [1 + jitter_offset, 2 + jitter_offset]';
% y_coords = [mean_resp_sound(idx_excited_sound), mean_resp_sound_light(idx_excited_sound)]';
% plot(x_coords, y_coords, 'k-')
% 
% 
% ax = nexttile(t,50,[4 1]);
% jitter_amount = 0.1;
% jitter_offset = (rand(sum(idx_excited_shock), 1) - 0.5) * jitter_amount;
% scatter(1 + jitter_offset, mean_resp_shock(idx_excited_shock))
% hold on
% scatter(2 + jitter_offset, mean_resp_shock_light(idx_excited_shock))
% x_coords = [1 + jitter_offset, 2 + jitter_offset]';
% y_coords = [mean_resp_shock(idx_excited_shock), mean_resp_shock_light(idx_excited_shock)]';
% plot(x_coords, y_coords, 'k-')


% 
idx_allexc = idx_excited_sound | idx_excited_shock;


ax = nexttile(t,49,[4 1]);
jitter_amount = 0.1;
jitter_offset = (rand(sum(idx_allexc), 1) - 0.5) * jitter_amount;
scatter(1 + jitter_offset, diff_preAP_all(idx_allexc))
hold on
scatter(2 + jitter_offset, diff_preAP_light_all(idx_allexc))
x_coords = [1 + jitter_offset, 2 + jitter_offset]';
y_coords = [diff_preAP_all(idx_allexc), diff_preAP_light_all(idx_allexc)]';
plot(x_coords, y_coords, 'k-')
hold off

ax = nexttile(t,50,[4 1]);
jitter_amount = 0.1;
jitter_offset = (rand(sum(idx_excited_sound), 1) - 0.5) * jitter_amount;
scatter(1 + jitter_offset,n_spx_resp_sound(idx_excited_sound))
hold on
scatter(2 + jitter_offset, n_spx_resp_sound_light(idx_excited_sound))
x_coords = [1 + jitter_offset, 2 + jitter_offset]';
y_coords = [n_spx_resp_sound(idx_excited_sound), n_spx_resp_sound_light(idx_excited_sound)]';
plot(x_coords, y_coords, 'k-')
hold off

ax = nexttile(t,51,[4 1]);
jitter_amount = 0.1;
jitter_offset = (rand(sum(idx_excited_shock), 1) - 0.5) * jitter_amount;
scatter(1 + jitter_offset, n_spx_resp_shock(idx_excited_shock))
hold on
scatter(2 + jitter_offset, n_spx_resp_shock_light(idx_excited_shock))
x_coords = [1 + jitter_offset, 2 + jitter_offset]';
y_coords = [n_spx_resp_shock(idx_excited_shock), n_spx_resp_shock_light(idx_excited_shock)]';
plot(x_coords, y_coords, 'k-')
hold off


signrank(n_spx_resp_sound(idx_excited_sound), n_spx_resp_sound_light(idx_excited_sound))
signrank(n_spx_resp_shock(idx_excited_shock), n_spx_resp_shock_light(idx_excited_shock))





%% (1,1:4) - Heatmaps
clearvars -except g t
ttl_all = {'triptest_sound_only', 'triptest_sound_only_light', 'triptest_shocks_only', 'triptest_shocks_only_light', 'triptest_both', 'triptest_both_light'};
hmptitles_all = {'CS', 'CS_light', 'US', 'US_light', 'CS + US', 'CS + US_light'};
mploc_all = [97 99 101 103 105 107];
mploc_all_2 = [145 145 149 149 153 153];

for ii = 1:3
    idx = (ii-1)*2 + (1:2);   % gives [1 2], [3 4], [5 6]
    ttl = ttl_all(idx);
    hmptitles = hmptitles_all(idx);
    mploc = mploc_all(idx);
    
    PSTHall = [];
    for hmp = 1:size(ttl,2)
        psth_spx_og =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
        psth_spx = zscore(psth_spx_og,0,2);   
        psth_spx = smoothdata(psth_spx,2,'sgolay', g.smoothvalue);
        idx_PN = strcmp(g.cell_metrics.brainRegion,'LA') & strcmp(g.cell_metrics.putativeCellType,'PN');
        psth_spx_og_PN = psth_spx_og(idx_PN,:);
        psth_spx_PN = psth_spx(idx_PN,:);
        hmpdata.(['psth_' num2str(hmp)]) = psth_spx_PN(:,g.roi_pca);
        PSTHall = [PSTHall hmpdata.(['psth_' num2str(hmp)])]; 
    end
    
    [~,PCA1]  = pca(PSTHall); % running pricipal component analysis
    PCA2      = PCA1(:,1:2); % extracting the first PCA_num pricinpal components 
    inNAN = nanmedian(PCA2,1); % to avoid NaN values, set them to the median
    for jj = 1:size(PSTHall,1)
        if isnan(PCA2(jj,1)) == 1
            PCA2(jj,1:3) = inNAN;
        end
    end
    Dend      = linkage(PCA2,'complete','mahalanobis'); % calculating the dendrogram
    Clusters  = cluster(Dend,'maxclust',g.clustnum); % clustering based on the tree
    D         = pdist(PCA2); %euclidean distrance between point in the artifical space
    leafOrder = optimalleaforder(Dend,D)'; % optimal order for plotting    
    PCA        = PCA2;
    Dendrogram = Dend;
    leafOrderfl = leafOrder;
    %leafOrderfl([54 57]) = leafOrderfl([57 54]);
    
    for hmp = 1:size(ttl,2)
        ax = nexttile(t,mploc(hmp),[4 2]);
        psth_spx =  BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', ttl{hmp}, 'pre_time', g.pre_time, 'post_time', g.post_time, 'bin_time', g.bin_time);
        psth_spx = zscore(psth_spx,0,2);   
        psth_spx = smoothdata(psth_spx,2,'sgolay', g.smoothvalue);
        idx_PN = strcmp(g.cell_metrics.brainRegion,'LA') & strcmp(g.cell_metrics.putativeCellType,'PN');
        psth_spx_PN = psth_spx(idx_PN,:);
        psth_spx_sorted = [psth_spx_PN(leafOrderfl,:)];
        matrix = psth_spx_sorted(:,(g.pre_time-g.plotwin(1))/g.bin_time+1:(g.pre_time+g.plotwin(2))/g.bin_time);
        imagesc(g.timeaxis,1:size(matrix,1),matrix);
        clim(g.clim)
        disp([-max(max(psth_spx_sorted,[],1))/2.5 max(max(psth_spx_sorted,[],1))])   
        colormap(g.colors.Heatmap); 
        yline(sum(idx_PN)+0.5, 'k', 'LineWidth', 1);
        xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha',1);
        ylabel('Cell number')
        set(gca, 'FontSize', g.fontSize2);
        title(hmptitles{hmp}, 'FontSize', g.fontSize1)
        if hmp == 1
            %cb = colorbar('westoutside', 'FontSize', g.fontSize2);
        end
        hold on
        n_clu = find(diff(Clusters(leafOrderfl))~=0);
        yline(n_clu+0.5, 'Color', 'k', 'LineWidth', 1);
        hold off
        mtxdata.(['psth_' num2str(hmp)]) = matrix;
    end



% Assume psth_spx_sorted_sound and psth_spx_sorted_sound_light have been created from the original code.

% Example: Calculate the difference and plot a new heatmap
ax = nexttile(t, mploc_all_2(idx(1)), [4 2]);
difference_matrix = (mtxdata.psth_2 - mtxdata.psth_1);
imagesc(g.timeaxis, 1:size(difference_matrix, 1), difference_matrix);
clim(g.clim)
colormap(g.colors.Heatmap);
%colorbar;
title('Difference: Sound vs. Sound + Light', 'FontSize', g.fontSize1);
xlabel('Time (s)');
ylabel('Cell number');
xline(0, '--k', 'LineWidth', g.xlinewidth, 'Alpha', 1);







end




