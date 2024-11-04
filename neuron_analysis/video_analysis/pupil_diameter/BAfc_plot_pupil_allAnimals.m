animalID = {'MD254', 'MD268', 'MD277', 'MD278'};

for ii = 1:numel(animalID)
    group_averages.(animalID{ii}) = BAfc_plot_pupil_changes_group(animalID{ii});
end


% at kell irni a BAfc_plot_pupil_changes_group scriptet, hogya rez file
% helyett a .mat filebol vegyek ki az est. fps-t. Ezutan mehet mar ez a
% loop gyorsabban