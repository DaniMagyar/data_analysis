LFPdata = loadLFPdata('resdir', [cd '\myresdir'], 'CHspec', 64, 'TTspec', 1:16,...
    'CHmap', 'Cambridge64_H2', 'reference', 'common_avg', 'lastTS', 250);

% downsampling nem kell, de CAR igen. Erdemes inkabb 200-250 pulzust plottolni


BA_500=BA_500-LFPdata.timestamps(1); % ttl timestamps normalized to data timestamps