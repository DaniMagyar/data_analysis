function [cell_metrics] = DM_load_LFP(varargin)

% Default params
prs = inputParser;
addRequired(prs,'g',@isstruct)
addRequired(prs,'cell_metrics',@isstruct)
addParameter(prs,'notch',false,@islogical)
parse(prs, varargin{:})

g = prs.Results.g;
g.notch = prs.Results.notch;
cell_metrics = prs.Results.cell_metrics;

channel_assignments = load([g.mainFolder '\channel_assignments.mat'], 'channel_assignments');
channel_assignments = channel_assignments.channel_assignments;
loadBar = waitbar(0,'Loading theta...');

for ii = 1:size(g.recordings,2)
    waitbar(ii/size(g.recordings,2), loadBar);
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.(g.ttl);
    mDatafilename = dir([cell_metrics.general.basepaths{ii} '\LFP\LFP*.mat']);
    mData = load([cell_metrics.general.basepaths{ii} '\LFP\' mDatafilename.name]);
    lfpFilename = dir([cell_metrics.general.basepaths{ii} '\LFP\LFP*.dat']);
    lfpFile = fopen([cell_metrics.general.basepaths{ii} '\LFP\' lfpFilename.name], 'r');
    lfpData = fread(lfpFile, [mData.numChannels, Inf], 'int16');
    fclose(lfpFile);
    if g.notch
        f0 = 60;        % Frequency to notch (60 Hz)
        Q = 35;         % Quality factor     
        % Design a bandstop IIR filter
        notchFilter = designfilt('bandstopiir', ...
            'FilterOrder', 2, ...
            'HalfPowerFrequency1', f0*(1 - 1/Q), ...
            'HalfPowerFrequency2', f0*(1 + 1/Q), ...
            'SampleRate', mData.fs);
    
        % Apply the notch filter to each channel both ways (filtfilt)
        filteredData = filtfilt(notchFilter, lfpData')';
    
        % Store filtered LFP data
        lfp.data = filteredData(1:64, :)';
    else
        lfp.data = lfpData(1:64, :)';
    end
    lfp.samplingRate = mData.fs;
    lfp.timestamps = (1:size(lfp.data,1))';
    lfp.timestamps = lfp.timestamps/mData.fs;
    [~, lfpAvg] = bz_eventCSD(lfp, TTL, 'spat_sm', 0, 'temp_sm',0, 'twin', g.twin, 'plotCSD', false, 'plotLFP', false);
    cell_metrics.LFP.(g.recordings{ii}).(g.ttl).data = lfpAvg.data';
    cell_metrics.LFP.(g.recordings{ii}).(g.ttl).timestamps = {lfpAvg.timestamps/lfpAvg.samplingRate};
    cell_metrics.LFP.(g.recordings{ii}).(g.ttl).channelID = (1:64)';
    cell_metrics.LFP.(g.recordings{ii}).(g.ttl).channelBR = channel_assignments.(g.recordings{ii}(1:5)).regions';
end

close(loadBar)