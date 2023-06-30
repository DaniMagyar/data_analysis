function Preprocess_OEPv063(varargin)

% This function removes artefacts from recording and applies channel map.

%  INPUT: 
% - pLength: length of the artefact in ms. Default is 5ms.
% - extracut: extra time cutted before and after the artefact in ms. Default 1ms
% - lastTS: Last TimeStamp included (default: end of the recording)
% - experiment: 'PFC_BA_doubleLaser' / 'M2_shock_laser'

prs =  inputParser;
addOptional(prs,'pLength',0,@isnumeric) 
addOptional(prs,'extracut',0,@isnumeric)
%addOptional(prs,'lastTS',0,@isnumeric)
addOptional(prs,'lastSec',0,@isnumeric) % start = 0 sec
addOptional(prs, 'removePeriod', [0 0],@isnumeric) % in seconds from 0 sec
addOptional(prs, 'experiment', '', @ischar)
addOptional(prs, 'source', 'none', @ischar)
addParameter(prs,'CHmap','none',@ischar)
parse(prs,varargin{:})
g = prs.Results;
fs = 30000;
disp(['OEP v0.6.3'])
disp(['Pulse length = ' num2str(g.pLength)])
disp(['Extracut = ' num2str(g.extracut)])

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

switch g.CHmap
    case 'none'
        error('No channel map selected')
    case 'preordered'
        channel_map = [1:64];
        disp('Selected channel map: PREORDERED')
    case 'A1x32'
        channel_map = [1 17 16 32 3 19 14 30 9 25 10 20 8 24 2 29 7 26 15 21 11 23 12 28 6 18 13 22 5 27 4 31];
        disp('Selected channel map: A1x32')
    case 'Buzs32'
        channel_map = [31 26 27 21 22 23 18 28 29 17 24 32 20 19 25 30 1 2 16 8 3 10 14 9 7 4 15 5 11 13 12 6];
        disp('Selected channel map: Buzs32')
    case 'Cambridge64_H2_Adpt_NNx'
        channel_map = [3 1 31 28 30 23 7 5 25 6 8 29 27 32 2 4 9 11 14 16 18 20 22 24 26 21 19 17 15 13 12 10 56 54 51 49 47 45 43 41 39 44 46 48 50 52 53 55 62 64 34 37 35 42 58 60 40 59 57 36 38 33 63 61];
        disp('Selected channel map: Cambridge64_H2_Adpt_NNx')
    case 'Cambridge64_H2_Adpt_Cambridge'
        channel_map = [14 13 12 7 6 57 1 2 9 15 16 11 10 5 4 3 49 50 51 52 53 54 55 56 8 58 59 60 61 62 63 64 48 47 46 45 44 43 42 41 25 39 38 37 36 35 34 33 19 20 21 26 27 40 32 31 24 18 17 22 23 28 29 30];
        disp('Selected channel map: Cambridge64_H2_Adpt_Cambridge')
    case 'Cambridge64_H7_Adpt_NNx'
        channel_map = [3 10 1 12 31 13 28 15 30 17 23 19 7 21 5 26 25 24 6 22 8 20 29 18 27 16 32 14 2 11 4 9 56 61 54 63 51 33 49 38 47 36 45 57 43 59 41 40 39 60 44 58 46 42 48 35 50 37 52 34 53 64 55 62];
        disp('Selected channel map: Cambridge64_H7_Adpt_NNx')
    case 'Cambridge64_H7_Adpt_Cambridge'
        channel_map = [14 64 13 63 12 62 7 61 6 60 57 59 1 58 2 8 9 56 15 55 16 54 11 53 10 52 5 51 4 50 3 49 48 30 47 29 46 28 45 23 44 22 43 17 42 18 41 24 25 31 39 32 38 40 37 27 36 26 35 21 34 20 33 19];
        disp('Selected channel map: Cambridge64_H7_Adpt_Cambridge')
    case 'Cambridge64_P1'
        channel_map = [23 3 7 1 5 31 25 28 6 30 8 29 27 32 2 4 20 9 22 11 24 14 26 16 21 18 19 17 15 13 12 10 45 56 43 54 41 51 39 49 44 47 46 48 50 52 53 55 42 62 58 64 60 34 40 37 59 35 57 36 38 33 63 61];
        disp('Selected channel map: Cambridge64_P1')
    case 'Cambridge64_H5_Adpt_Cambridge'
        channel_map = [50 14 51 13 52 12 53 7 54 6 55 57 56 1 8 2 58 9 59 15 60 16 61 11 62 10 63 5 64 4 48 3 47 49 46 45 44 43 42 41 25 39 38 37 36 35 34 33 19 20 21 26 27 40 32 31 24 18 17 22 23 28 29 30];
        disp('Selected channel map: Cambridge64_H5_Adpt_Cambridge')
    case 'Cambridge64_H10_Adpt_Cambridge'
        channel_map = [49 12 62 3 13 63 50 14 64 51 7 61 52 6 60 53 57 59 54 1 58 4 2 8 5 9 56 10 15 55 11 16 19 46 28 33 47 29 20 48 30 21 45 23 26 44 22 27 43 17 40 42 18 34 41 24 35 25 31 36 39 32 37 38];
        disp('Selected channel map: Cambridge64_H10_Adpt_Cambridge')
    otherwise
        error('read_openephys: Unknown channel map.')
        % [~,idx]=ismember(mapH2nnx,mapH7nnx) h2nnx  to h7nnx konverzio idx
        % mapH2cmb = mapH7cmb(idx) h2cmb to h7cmb with same idx

end  

if g.lastSec > 0
%     g.lastTS = (info.Timestamps(1)+g.lastSec*30000)/30000;
%     lastEL = find(info.Timestamps==g.lastTS*30000);
    lastEL = g.lastSec*30000;
else
    lastEL = length(info.Timestamps);
end

% initialize a big matrix for the data
incoming_data = zeros(info.Header.num_channels,length(info.Timestamps(1:lastEL)),'int16');
% load data
loadBar = waitbar(0,'Loading channels...');
for ii = 1:info.Header.num_channels
    [next_channel] = info.Data.Data(1).mapped(channel_map(ii),1:lastEL);
    incoming_data(ii,:)=int16(next_channel);
    waitbar(ii/64,loadBar);
end
close(loadBar)
% % to improve: try to load with parfor, and use tic toc
% pwb_load = parwaitbar(info.Header.num_channels);
% parfor ii = 1:info.Header.num_channels
%     incoming_data(ii,:) = int16(info.Data.Data(1).mapped(channel_map(ii),1:lastEL));
%     pwb_load.progress();
% end



%% Remove stimulation artefacts
tic
if g.pLength ~= 0
    channel_states = readNPY([cd '\events\Acquisition_Board-100.Rhythm Data\TTL\states.npy']);
    event_sample_numbers = readNPY([cd '\events\Acquisition_Board-100.Rhythm Data\TTL\sample_numbers.npy']);
    data_sample_numbers = readNPY([cd '\continuous\Acquisition_Board-100.Rhythm Data\sample_numbers.npy']);
    
    for jj = 1:max(channel_states)
        TTL_ON = find(channel_states==int16(jj));
        TTL_channels{1,jj*2-1} = ['CH' num2str(jj) 'ON_sample_num'];
        TTL_channels{2,jj*2-1} = event_sample_numbers(TTL_ON);
        
        TTL_OFF = find(channel_states==-int16(jj));
        TTL_channels{1,jj*2} = ['CH' num2str(jj) 'OFF_sample_num'];
        TTL_channels{2,jj*2} = event_sample_numbers(TTL_OFF);
    end
    
    switch g.experiment
        case 'PFC_shock_run_rest'
            pulses_sn = sort([TTL_channels{2,1}]); % number of samples where stimulation starts
        case 'PFC_BAopto_run_rest'
             pulses_sn = sort([TTL_channels{2,1}]);
        otherwise 
            error('unknown experiment')
    end
    pulses_sn = (pulses_sn - data_sample_numbers(1)) +1; % normalize sample number; start from 1 instead of 0
    if g.lastSec > 0
        pulses_sn = pulses_sn(pulses_sn<g.lastSec*30000);
    end
    pwb = parwaitbar(info.Header.num_channels);
    parfor ll = 1:info.Header.num_channels
        tmp = incoming_data(ll,:);
        for kk = 1:numel(pulses_sn)
            pStamp = pulses_sn(kk); % sample_num of current pulse
            
            tmp(pStamp-g.extracut*30:pStamp+(g.pLength+g.extracut)*30) = ...
                int16(linspace(double(tmp(pStamp-g.extracut*30)), ...
                double(tmp(pStamp+(g.pLength+g.extracut)*30)), ...
                numel(tmp(pStamp-g.extracut*30:pStamp+(g.pLength+g.extracut)*30))));
        end
        incoming_data2(ll,:) = tmp;
        pwb.progress();
    end
else 
    incoming_data2 = incoming_data;
end
toc

%% Remove defined time periods (long artefacts)
if g.removePeriod ~= [0 0]

    rPLength = g.removePeriod(2) - g.removePeriod(1);
    incoming_data2(:,(g.removePeriod(1)*fs):(g.removePeriod(2)*fs)) = incoming_data2(:,((g.removePeriod(1)-rPLength)*fs):((g.removePeriod(2)-rPLength)*fs));
    disp('Artefact period removed')
end


%% initialize a new .dat file
if g.removePeriod ~= [0 0]
    fid=fopen(['continuous_Preprocessed_' g.CHmap '_pLength' num2str(g.pLength) '_extracut' num2str(g.extracut) 'OverwritedFrom' num2str(g.removePeriod(1)) 'To' num2str(g.removePeriod(2)) '.dat'],'w+');
else 
    fid=fopen(['continuous_Preprocessed_' g.CHmap '_pLength' num2str(g.pLength) '_extracut' num2str(g.extracut) '.dat'],'w+');
end
% write data to raw .dat file
fwrite(fid,incoming_data2(:,:),'int16');
% close that file
fclose(fid);
disp('done')
