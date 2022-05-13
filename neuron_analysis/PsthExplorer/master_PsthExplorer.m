function master_PsthExplorer(Recordings, Stim, preferences)

% Input
% -Recordings: vertical list of included recordings in cell format, first 5 char
% -Stim: variable to use from TTLsKS.mat, string
% GUI related
% Time axis setup
preferences.time = (preferences.int(1)*preferences.fs:preferences.psth_bin:preferences.int(2)*preferences.fs)';
preferences.time = preferences.time(1:end-1);

% Colormap blue:white:red -------------------------------------------------
c1 = 1/255*[0,115,185];
c2 = 1/255*[239,239,239];
c3 = 1/255*[217,84,26];
preferences.mycolormap(1:20,1)  = linspace(c1(1),c2(1),20);
preferences.mycolormap(21:64,1) = linspace(c2(1),c3(1),44);
preferences.mycolormap(1:20,2)  = linspace(c1(2),c2(2),20);
preferences.mycolormap(21:64,2) = linspace(c2(2),c3(2),44);
preferences.mycolormap(1:20,3)  = linspace(c1(3),c2(3),20);
preferences.mycolormap(21:64,3) = linspace(c2(3),c3(3),44);

% preferences = preferences_PsthExplorer;
PSTHall=[];
Wilcoxon_results = [];

for ii = 1:length(Recordings)
    cd([preferences.mainFolder '\' Recordings{ii} '_kilosort\kilosort3preprocess'])

    cluster_group = tdfread('cluster_group.tsv');
    cluster_info = tdfread('cluster_info.tsv');
%     if length(cluster_group.group) ~= length(cluster_info.group)
%         error('Not all clusters have group label in Phy')
%     else
%         disp('All groups labeled')
%     end
    
    spike_times = readNPY('spike_times.npy');
    spike_clusters = readNPY('spike_clusters.npy');
    
    groupCell = cellstr(cluster_info.group);
    goodIdx = find(strcmp(groupCell, 'good')==1); % to include MUA: find(strcmp(groupCell, 'noise')==0)
    
    for kk = 1:length(goodIdx)
        ttl = load_TTL(Stim); % Load TTLs & Choose variable to use as stimulus
        spk_idx = find(spike_clusters==cluster_info.cluster_id(goodIdx(kk)));
        spk_t = double(spike_times(spk_idx));
        spikes = spk_t/preferences.fs;
        psth_spikes =  calculate_psth(spikes, ttl, preferences);
        Wilcoxon_resCurr{kk,2} = calculate_Wilcoxon(spikes, ttl, preferences);
        if preferences.norm == 1
            newpsth = zscore(psth_spikes);
        else
            newpsth = psth_spikes;
        end
        PSTHcurr(kk,:) = newpsth';
    end
    PSTHall = vertcat(PSTHall, PSTHcurr);
    Wilcoxon_results =  vertcat(Wilcoxon_results, Wilcoxon_resCurr);
    clear Wilcoxon_resCurr PSTHcurr;
end

switch preferences.plottingMethod
    case 'Zscore'
        [SortIDX, testMean] = calculate_ZscoreOrder(PSTHall, Stim, preferences);
        plot_PsthZscore(Stim, PSTHall, Wilcoxon_results, SortIDX, testMean, preferences)
    case 'PCA'
        [leafOrder, Dend, Clusters] = calculate_PCAOrder(PSTHall, preferences);
        plot_PsthPCA(PSTHall, leafOrder, Dend, Clusters, preferences)
end
