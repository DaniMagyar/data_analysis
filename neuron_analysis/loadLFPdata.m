function LFPdata = loadLFPdata(varargin)
 
prs = inputParser;

addOptional(prs,'datadir',cd,@(s)isempty(s)|isdir(s))   % data directory
addOptional(prs,'resdir','',@(s)isempty(s)|isdir(s))   % results directory
addOptional(prs,'CHspec',32,@isnumeric)   % Number of channels (default: 32 channels)
addOptional(prs,'TTspec',1:8,@isnumeric)   % Number of tetrodes (default, 8 tetrodes)
addOptional(prs,'rawdatafiletag','',@ischar)   % switch for filtering
addOptional(prs,'processor',101,@isnumeric) % Processor number, default 101
addOptional(prs,'lastTS',0,@isnumeric) % Last TimeStamp included (default: end of the recording)
addOptional(prs,'sampling',1,@isnumeric) % downsample the original 30kHz data 
addParameter(prs,'CHmap','',@ischar)
addParameter(prs,'reference','common_avg',@(s)ischar(s)|isempty(s))   % switch for referencing
addOptional(prs,'exclude_channels',[],@isnumeric)   % exclude these channels from averaging
parse(prs,varargin{:})
g = prs.Results;

Th = 25;  % Threshold for the spike detection
disp(['Selected threshold: ' num2str(Th)])

switch g.CHmap 
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

if g.lastTS > 0
    lastEL = find(info.Timestamps==g.lastTS*30000);
else
    lastEL = length(info.Timestamps);
end
 

 % Tetrode organization of the channels
if nargin < 3 || isempty(g.CHspec)
    g.CHspec = 32;
    NumTetrodes = g.CHspec / 4;
    g.TTspec = 1:NumTetrodes;
end

% Directories
if isequal(g.datadir(end),'\')
    g.datadir = g.datadir(1:end-1);
end

SaveFeatures = true;   % save MClust feature files

% Common average reference
switch g.reference
    case 'common_avg'
        common_avg = common_avg_ref_probe_dat(g.datadir,lastEL,info,'channel_number', g.CHspec); 
    case {'','none'}
        common_avg = 0;
    otherwise
        error('read_openephys: Unknown reference option.')
end

% Load LFPdata

%LFPdata.data = zeros(g.CHspec,find(info.Timestamps==g.lastTS*30000));
LFPdata.timestamps = double(info.Timestamps(1:find(info.Timestamps==g.lastTS*30000)))/30000;
LFPdata.timestamps = LFPdata.timestamps(1:g.sampling:end);
LFPdata.samplingRate = 30000/g.sampling;
for iX = 1:g.CHspec
    if ~ismember(iX,g.exclude_channels)   % exclude channels
        data = info.Data.Data(1).mapped(channel_map(iX),1:lastEL);
        data = double(data)'*0.195; %intan
        data = data - common_avg(:,1);
        data = data(1:g.sampling:end); % downsampling
        LFPdata.data(:,iX) = data;
        disp(['LFPdata - ' num2str(iX)])
    end
end

