function BAfc_import_brainRegion

% importes brainRegion data from channel_assignments.mat to cell_metrics.
% The channel_assignments.mat file should be in the mainFolder, that 
% contains subfolders of different experiments 
% (e.g. mainFolder/MD243_kilosort/kilosort25preprocess/)

mainFolder = 'C:\Users\dmagyar\Desktop\BA_fear_cond';
cd(mainFolder)
load('channel_assignments.mat')
folders = dir;
folders = folders([folders.isdir]); % Keep only directories
folders = folders(~ismember({folders.name}, {'.', '..'})); % Remove '.' and '..'
folderNames = {folders.name}; % Extract folder names

for jj = 1:size(folderNames,2)
    cd([mainFolder '/' folderNames{jj} '/kilosort25preprocess'])
    load('temp_wh.cell_metrics.cellinfo.mat')
    for ii = 1:length(cell_metrics.brainRegion)
        idx = find(channel_assignments.(folderNames{jj}(1:5)).channels == cell_metrics.maxWaveformCh1(ii));
        cell_metrics.brainRegion(ii) = channel_assignments.(folderNames{jj}(1:5)).regions(idx);
    end
    save('temp_wh.cell_metrics.cellinfo.mat', 'cell_metrics')
    disp([folderNames{jj}(1:5) ': done'])
    clearvars -except channel_assignments mainFolder folderNames
end