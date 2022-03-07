function ClearAndRemapForKilosort(varargin)

% This function removes artefacts from recording and applies channel map.

% Optional INPUT: 
% - timestamps: timestamp where artefact begins
% - pLength: length of the artefact in ms. Default is 5ms.
% - extracut: extra time cutted before and after the artefact in ms. Default 1ms

prs =  inputParser;
addOptional(prs,'pLength',5,@isnumeric)
addOptional(prs,'extracut',1,@isnumeric)
addParameter(prs,'CHmap','none',@ischar)
parse(prs,varargin{:})
g = prs.Results;

switch g.CHmap
    case 'none'
        channel_map = [1:64];
        disp('No channel map selected')
    case 'A1x32'
        channel_map = [1 17 16 32 3 19 14 30 9 25 10 20 8 24 2 29 7 26 15 21 11 23 12 28 6 18 13 22 5 27 4 31];
        pairsNeed = 1; % we need all pairs on a linear probe
        disp('Selected channel map: A1x32')
    case 'Buzs32'
        channel_map = [31 26 27 21 22 23 18 28 29 17 24 32 20 19 25 30 1 2 16 8 3 10 14 9 7 4 15 5 11 13 12 6];
        pairsNeed = 2; % we need only every second pair, because others would be from different shanks
        disp('Selected channel map: Buzs32')
    case 'Cambridge64_H2'
        channel_map = [3 1 31 28 30 23 7 5 25 6 8 29 27 32 2 4 9 11 14 16 18 20 22 24 26 21 19 17 15 13 12 10 56 54 51 49 47 45 43 41 39 44 46 48 50 52 53 55 62 64 34 37 35 42 58 60 40 59 57 36 38 33 63 61];
        pairsNeed = 1;
        disp('Selected channel map: Cambridge64_H2')
    case 'Cambridge64_H7'
        channel_map = [3 10 1 12 31 13 28 15 30 17 23 19 7 21 5 26 25 24 6 22 8 20 29 18 27 16 32 14 2 11 4 9 56 61 54 63 51 33 49 38 47 36 45 57 43 59 41 40 39 60 44 58 46 42 48 35 50 37 52 34 53 64 55 62];
        pairsNeed = 1;
        disp('Selected channel map: Cambridge64_H7')        
    otherwise
        error('read_openephys: Unknown channel map.')
end  

info = load_open_ephys_binary([cd '\structure.oebin'],'continuous',1, 'mmap');

% initialize a big matrix for the data
incoming_data = zeros(info.Header.num_channels,length(info.Timestamps),'int16');
% load data
for ii = 1:info.Header.num_channels
    [next_channel] = info.Data.Data(1).mapped(channel_map(ii),:);
    incoming_data(ii,:)=int16(next_channel);
    disp(['CH ' num2str(channel_map(ii)) ' loaded'])
end

% find artefacts in folder
channel_states = readNPY([cd '\events\Rhythm_FPGA-100.0\TTL_1\channel_states.npy']);
channels = readNPY([cd '\events\Rhythm_FPGA-100.0\TTL_1\channels.npy']);
timestamps = double(readNPY([cd '\events\Rhythm_FPGA-100.0\TTL_1\timestamps.npy']));
for jj = 1:max(channels)
    TTL_ON = find(channel_states==jj);
    TTL_channels{1,jj} = ['CH' num2str(jj)];
    TTL_channels{2,jj} = timestamps(TTL_ON);
end
pulses = (TTL_channels{2,3});

for kk = 1:numel(pulses)
    pStamp = find(info.Timestamps==pulses(kk)); % timestamp of current pulse
    incoming_data(:,pStamp-g.extracut*30:pStamp+(g.pLength+g.extracut)*30) = 0;
    disp(['artefact ' num2str(kk) '/' num2str(numel(pulses)) ' removed'])
end

%initialize a new .dat file
fid=fopen(['continuous_Clear_Map' g.CHmap '.dat'],'w+');
% write data to raw .dat file
fwrite(fid,incoming_data(:,:),'int16');
% close that file
fclose(fid);
disp('done')