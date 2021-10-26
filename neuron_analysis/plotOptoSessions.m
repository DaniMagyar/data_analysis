savedFreqs = who('-file','optoTTL.mat');

for ii = 1:numel(savedFreqs)
    load('optoTTL.mat', savedFreqs{ii})
    groups=(1:16);
    allneurons_psth_dani_dat(groups, eval(savedFreqs{ii}), 'pre', 0.05, 'post', 0.05, 'bin', 60, 's', 'figname', savedFreqs{ii})
end

