function PostProcessAfterKilosort(varargin)

% This function removes artifact contamination from clusters

prs =  inputParser;
addOptional(prs,'pLength',2,@isnumeric) % where to look for contamination
addOptional(prs,'extracut',1,@isnumeric)
addOptional(prs, 'experiment', '', @ischar)
addOptional(prs, 'within_unit_overlap_window',0,@isnumeric)
addOptional(prs, 'between_unit_overlap_window',0,@isnumeric)
addOptional(prs, 'fs', 30000,@isnumeric)
parse(prs,varargin{:})
g = prs.Results;


%% load shock timestamps 
info = load_open_ephys_binary([cd '\structure.oebin'],'continuous',1, 'mmap');
zeroTime = info.Timestamps(1);

channel_states = readNPY([cd '\events\Rhythm_FPGA-100.0\TTL_1\channel_states.npy']);
channels = readNPY([cd '\events\Rhythm_FPGA-100.0\TTL_1\channels.npy']);
timestamps = double(readNPY([cd '\events\Rhythm_FPGA-100.0\TTL_1\timestamps.npy']));

for ii = 1:max(channels)
    TTL_ON = find(channel_states==ii);
    TTL_channels{1,ii} = ['CH' num2str(ii)];
    TTL_channels{2,ii} = timestamps(TTL_ON); % Here only ON-s are collected
end

switch g.experiment
    case 'M2_shock_only'
        manipulation = (TTL_channels{2,1}-double(zeroTime)); % from 0 time start point
    case 'PFC_ChETA'
        manipulation = (TTL_channels{2,3}-double(zeroTime)); % from 0 time start point. 
    case 'PFC_shock_laser'
        manipulation = (TTL_channels{2,2}-double(zeroTime));

end

cd ([cd '\kilosort3preprocess']) % move into folder 
amplitudes = readNPY('amplitudes.npy');
spike_times = readNPY('spike_times.npy'); % from 0 time start point
spike_clusters = readNPY('spike_clusters.npy');
spike_templates = readNPY('spike_templates.npy');

%% remove artifacts from cluster
for jj = 1:numel(manipulation)
    idx_1{:,jj} = find(ismember(spike_times, spike_times(spike_times>(manipulation(jj)-g.extracut*30)...
        & spike_times<(manipulation(jj)+(g.pLength+g.extracut)*30)))); % find artifact timestamps
end

artifact_idx = vertcat(idx_1{1:numel(idx_1)}); % artifact indices in NPY files

spike_times(artifact_idx) = [];
amplitudes(artifact_idx) = [];
spike_clusters(artifact_idx) = [];
spike_templates(artifact_idx) = [];

%% remove remaining high amplitude events from clusters
[~, ClusterIDs] = groupcounts(spike_clusters);
for kk = 1:numel(ClusterIDs)
    idx_2 = find(spike_clusters == ClusterIDs(kk)); % spike indices of current cluster in .npy files
    ClusterAmps = amplitudes(idx_2);
    ClusterMean = mean(ClusterAmps);
    ClusterStd =  std(ClusterAmps);
    idx_3 = find(ismember(ClusterAmps, ClusterAmps(ClusterAmps<(ClusterMean-3*ClusterStd) | ...
        ClusterAmps>(ClusterMean+3*ClusterStd)))); % events outside 3std removed
    EventIdx{:,kk} = idx_2(idx_3); % high amplitude event indices in .npy files
end

event_idx =  vertcat(EventIdx{1:numel(EventIdx)});

spike_times(event_idx) = [];
amplitudes(event_idx) = [];
spike_clusters(event_idx) = [];
spike_templates(event_idx) = [];

%% remove double detected spikes

templates = readNPY('templates.npy');
for mxch = 1:length(templates(:,1,1))
    [~, maxChannels(mxch,1)] = max(max(templates(mxch,:,:),[], 2));
end
maxChannels = maxChannels-1; % because channel_map starts from 0. Ez a peak_channels a pyton scriptbe. 
                             % Ott nincs preordered chanel map, ezert mas kicsit. 


% A Phy-ben elmeletben az elso cluster 0 ID-t kap, gyakorlatban ez
% valtozhat. Ami a Phy-ben a id, ahhoz +1-et hozza kell adni, mivel a
% matlabban az 1. index 1. Az igy kapott szam adja meg hogy a maxChannels
% hanyadik soraban talaljuk, hogy melik csatornan latszik a sejt a
% legnagyobb amplitudoval. 




%% remove within unit overlapping spikes
if g.within_unit_overlap_window > 0
    within_unit_overlap_samples = g.within_unit_overlap_window*(g.fs/1000); % /1000 because in miliseconds
    for WUS = 1:length(ClusterIDs)
        idx_4 = find(spike_clusters==ClusterIDs(WUS));
        cluSpkTimes = spike_times(idx_4);
        cluSpkDiff = cluSpkTimes(2:end) - cluSpkTimes(1:end-1); 
        idx_5 = find(cluSpkDiff<within_unit_overlap_samples); 
        WUremove_idx{:,WUS} = idx_4(idx_5);
    end
    
    WU_double_idx = vertcat(WUremove_idx{1:numel(WUremove_idx)});
    
    spike_times(WU_double_idx) = [];
    amplitudes(WU_double_idx) = [];
    spike_clusters(WU_double_idx) = [];
    spike_templates(WU_double_idx) = [];
end

%% remove between unit overlapping spikes
%between_unit_overlap_samples = g.between_unit_overlap_window*(g.fs/1000);

writeNPY(spike_times, 'spike_times.npy');
writeNPY(amplitudes, 'amplitudes.npy');
writeNPY(spike_clusters, 'spike_clusters.npy');
writeNPY(spike_templates, 'spike_templates.npy');

disp('done')



