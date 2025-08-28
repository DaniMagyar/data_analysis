% calculating
clear all

cd('C:\Users\dmagyar\My Drive\Papers\Zsofi')

for ii = 1:size(cell_metrics.cellID,2)
    cell_metrics.nSpikes(ii) = numel(cell_metrics.spikes.times{ii});
    cell_metrics.isi{ii} = diff(cell_metrics.spikes.times{ii});
    cell_metrics.t_active(ii) = round(cell_metrics.spikes.times{ii}(end) - cell_metrics.spikes.times{ii}(1),2);
    cell_metrics.mean_fr(ii) = round(cell_metrics.nSpikes(ii)/cell_metrics.t_active(ii),2);
    cell_metrics.mean_isi(ii) = round(mean(cell_metrics.isi{ii}),2);
    cell_metrics.std_isi(ii) = round(std(cell_metrics.isi{ii}),2);
    cell_metrics.cv_isi(ii) = round(cell_metrics.std_isi(ii)/cell_metrics.mean_isi(ii),2);
end
