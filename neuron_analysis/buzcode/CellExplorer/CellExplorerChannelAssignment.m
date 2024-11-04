% modify 2 times

for ii = 1:length(cell_metrics.brainRegion)
    idx = find(channel_assignments.MD278.channels == cell_metrics.maxWaveformCh1(ii));
    cell_metrics.brainRegion(ii) = channel_assignments.MD278.regions(idx);
end
save('temp_wh.cell_metrics.cellinfo.mat', 'cell_metrics')
disp('done')
clearvars -except channel_assignments