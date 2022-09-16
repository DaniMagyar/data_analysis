function PostProcessAfterKilosort(varargin)

% This function removes artifact contamination from clusters

prs =  inputParser;
addOptional(prs,'pLength',6,@isnumeric) % where to look for contamination
addOptional(prs,'extracut',1,@isnumeric)
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
    TTL_channels{2,ii} = timestamps(TTL_ON);
end
shocks = (TTL_channels{2,1}-double(zeroTime)); % from 0 time start point


cd ([cd '\kilosort3preprocess']) % move into folder 
amplitudes = readNPY('amplitudes.npy');
spike_times = readNPY('spike_times.npy'); % from 0 time start point
spike_clusters = readNPY('spike_clusters.npy');
spike_templates = readNPY('spike_templates.npy');

%% remove artifacts from cluster
for jj = 1:numel(shocks)
    idx_1{:,jj} = find(ismember(spike_times, spike_times(spike_times>(shocks(jj)-g.extracut*30)...
        & spike_times<(shocks(jj)+(g.pLength+g.extracut)*30)))); % find artifact timestamps
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
        ClusterAmps>(ClusterMean+3*ClusterStd))));
    EventIdx{:,kk} = idx_2(idx_3); % high amplitude event indices in .npy files
end

event_idx =  vertcat(EventIdx{1:numel(EventIdx)});

spike_times(event_idx) = [];
amplitudes(event_idx) = [];
spike_clusters(event_idx) = [];
spike_templates(event_idx) = [];


writeNPY(spike_times, 'spike_times.npy');
writeNPY(amplitudes, 'amplitudes.npy');
writeNPY(spike_clusters, 'spike_clusters.npy');
writeNPY(spike_templates, 'spike_templates.npy');

disp('done')



