function [outgoing_data] = Preprocess_OEPv063_v7(varargin)

% This function removes artefacts from recording and applies channel map.
% v3: CAR and filtering added
% v4: filtfilt added
% v6: channel mapping loop is removed, much faster
% v7: lowpass, highpass and bandpass filter available, 'none' channel map
% avaialable, DM_load_channel_map, DM_load_artefacts added

%  INPUT: 
% - pLength: length of the artefact in ms. Default is 5ms.
% - extracut: extra time cutted before and after the artefact in ms. Default 1ms
% - firstSec: first extracted second from the recording. Counting starts at 0. 
% - lastSec: last extracted second from the recording. Counting starts at 0.
% - removePeriod: to remove long, noisy periods. These can confuse Kilosort. 
% - experiment: 'PFC_BA_doubleLaser' / 'M2_shock_laser'
% - source: source file, usually 'oebin'
% - CHmap: channel map

prs =  inputParser;
addParameter(prs,'pLength',0,@isnumeric) 
addParameter(prs,'extracut',0,@isnumeric)
addParameter(prs,'firstSec',0) % start = 0sec
addParameter(prs,'lastSec',0,@isnumeric) % start = 0 sec
addParameter(prs,'experiment','',@ischar)
addParameter(prs,'source','oebin',@ischar)
addParameter(prs,'CAR', [],@isnumeric) % common average referencing. 
addParameter(prs,'filter',[],@ischar) % filter type: 'low', 'high', 'bandpass'
addParameter(prs,'cutoff',[],@isnumeric) % single value for highpass and lowpass, and two value vector for bandpass (e.g.[300 6000])
addParameter(prs,'fs', 30000,@isnumeric) %sampling rate
addParameter(prs,'CHmap','none',@ischar) % channel map
addParameter(prs,'downsamplefactor',[],@isnumeric) % downsampling
addParameter(prs,'savename',[],@ischar) % savename e.g. 'mapped_CAR'
addParameter(prs,'savepath',[],@ischar) %savepath
addParameter(prs,'detrend',[],@isnumeric) % two time points (in ms) between which the detrending will be performed. (e.g. [0.015 0.03]) 
parse(prs,varargin{:})
g = prs.Results;
disp(['OEP v0.6.3'])
disp(['Pulse length = ' num2str(g.pLength)])
disp(['Extracut = ' num2str(g.extracut)])
%% Load file and channel map
switch g.source
    case 'oebin'
        info = load_open_ephys_binary([cd '\structure.oebin'],'continuous',1, 'mmap');
        disp('Selected source: OEBIN')
    case 'dat' % in this case, timestamps files must be overwrited and 'events' folder copied into cd
        folder = cd;
        info.Header.num_channels = 64;
        contFile=fullfile(folder,'continuous.dat');
        file=dir(contFile);
        samples=file.bytes/2/info.Header.num_channels;
        info.Data=memmapfile(contFile,'Format',{'int16' [info.Header.num_channels samples] 'mapped'});
        info.Timestamps = readNPY(fullfile(folder,'timestamps.npy'));
        info.Timestamps = info.Timestamps(1:samples);
        disp('Selected source: DAT')
    case 'none'
        error('read_openephys: Unknown source.')
end
channel_map = DM_load_channel_map(g.CHmap);
g.numChannels = numel(channel_map);
g.OriginalHeader = info.Header;
g.filename = info.Data.Filename;
%% Defining time window for preprocessing
if g.firstSec > 0
    firstEL = g.firstSec*30000;
else
    firstEL = 1;
end
if g.lastSec > 0
    lastEL = g.lastSec*30000;
else
    lastEL = length(info.Timestamps);
end
tic
disp('Loading data...')
incoming_data = info.Data.Data(1).mapped(channel_map,firstEL:lastEL);
%% Detrending before anything else
if g.detrend
    disp('Detrending data...')
    incoming_data = DM_detrend_data(g, incoming_data);
end
%% Apply filters
if g.filter
    disp('Filtering data...')
    [b_high, a_high] = butter(2,g.cutoff/g.fs*2,g.filter);
    loadBar = waitbar(0,'Filtering data...');
    for fn = 1:info.Header.num_channels
        incoming_data(fn,:) = int16(filtfilt(b_high,a_high,double(incoming_data(fn,:)))); % filter for one direction, filtfilt for zero-phase
        waitbar(fn/info.Header.num_channels,loadBar);
    end  
    close(loadBar)
end
%% Common average referencing
if g.CAR
    disp('Applying Common Average Referencing...')
    incoming_data = incoming_data - int16(mean(incoming_data, 1));
end
%% Remove stimulation artefacts
if any(g.pLength)
    disp('Removing artefacts...')
    artefacts_sn = DM_load_artefacts(g);
    parfor ll = 1:info.Header.num_channels
        tmp = incoming_data(ll,:);
        for kk = 1:numel(artefacts_sn)
            pStamp = artefacts_sn(kk); % sample_num of current pulse            
            tmp(pStamp-g.extracut*30:pStamp+(g.pLength+g.extracut)*30) = ...
                int16(linspace(double(tmp(pStamp-g.extracut*30)), ...
                double(tmp(pStamp+(g.pLength+g.extracut)*30)), ...
                numel(tmp(pStamp-g.extracut*30:pStamp+(g.pLength+g.extracut)*30))));
        end
        outgoing_data(ll,:) = tmp;
    end
else 
    outgoing_data = incoming_data;
end
%% Downsampling
if any(g.downsamplefactor)
    disp('Downsampling data...')
    outgoing_data = arrayfun(@(row) downsample(outgoing_data(row, :), g.downsamplefactor), 1:size(outgoing_data, 1), 'UniformOutput', false);
    g.fs = g.fs/g.downsamplefactor;
    outgoing_data = cell2mat(outgoing_data');
end

%% initialize a new .dat file, write data and close. Save 
if g.savepath
    cd (g.savepath)
end
if any(g.savename)
    savename = g.savename;
    fid=fopen([savename '.dat'],'w+');
    fwrite(fid,outgoing_data(:,:),'int16');
    fclose(fid);
    save([savename '.mat'], '-struct', 'g')
    disp(['Saved as: ' g.savename])
end
toc
disp('Data loaded.')