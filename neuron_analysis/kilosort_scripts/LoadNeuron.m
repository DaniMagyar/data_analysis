function [TS] = LoadNeuron(cellID)

% input
% -cellID: ID in CellExplorer

fs = 30000; %sampling rate
spike_times = readNPY('spike_times.npy');
spike_clusters = readNPY('spike_clusters.npy');
load('temp_wh.cell_metrics.cellinfo.mat');

cluID = cell_metrics.cluID(find(cell_metrics.cellID==cellID));
spk_idx = find(spike_clusters==cluID);
spk_t = double(spike_times(spk_idx));
TS = spk_t/fs;
disp([num2str(numel(TS)) ' spikes loaded'])
disp(['ClusterID in Phy: ' num2str(cluID)])