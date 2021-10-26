function TTL = load_ttl_dani_dat(varargin)

% This function reads a file (.continuous from OpenEphys) or data+time
% variables (.mat generated from an .rhd from Intan) or data+time in the
% Workspace (loaded audio file) and detects values above a threshold as
% timestamps for the TTL.  
% 
% INPUTS: 
%   Varargin: 
%   - file: signal to analyse, string. e. g. '101_ADC1.continuous'
%   - 'origin': where the  signal comes from. CASE 1: origin = 'OE',
%   recorded with OpenEphys; CASE 2: origin = 'IN', recorded with Intan
%   Board, channels saved from the .rhd file; CASE 3: origin = 'WS', for
%   data in the Workspace. DEFAULT 'OE'.
%   - 'n': number of TTLs. DEFAULT 1000 
%   - 'mPP': MinPeakProminence, theshold for the detection. DEFAULT 0.5.
%   The data is normalised in a 0 to 1 range. If you are not sure about the
%   threshold: 1) read the data to Workspace; 2) normalize it: data2 =
%   normalize(data,'range',[0 1]); 3) plot(data2); 4) Choose value.
%   - 'mPD': MinPeakDistance, min time between peaks, in seconds. DEFAULT
%   0.1 
%   - 'mPW': MinPeakWidth, min width of the peak, in seconds. DEFAULT 0.01 
%   - 's': The TTL is automatically saved in the current folder. DEFAULT
%   not save.
%   - 'ttl_dur': Duration of the ttl in sec. Calculates also the last point
%   of the peak, to double-check that the length is correct. DEFAULT 11s.
%
% OUTPUTS: 
%   - TTL: a vector list of the timestamps of the TTLs.
% 
% Examples: 
% TTL = load_ttl('101_ADC1.continuous');
% TTL = load_ttl('TTL.mat','origin','IN','n',60,'mPP',0.5,'mPD',0.1,'s');       
% TTL_light_video = load_ttl('file.mat','origin','WS','n',60,'mPP',0.5,...
%       'mPD',0.1,'mPW',10);
% TTL_light_video = load_ttl('file.mat','origin','WS','n',60,'mPP',0.5,...
%       'mPD',2,'mPW',2,'ttl_dur',11);
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% MATLAB toolbox: Signal Processing Toolbox.
% -------------------------------------------------------------------------

% Params, defalut values:
n = 1000;
mPP = 0.5;     
mPD = 0.05;
mPW = 0.001;
s = 1;
origin = 'OE';
ttl_dur = 11; 

info = load_open_ephys_binary([cd '\structure.oebin'],'continuous',1, 'mmap'); % for .dat files

% Params, read from varargin
if nargin
  for i=1:size(varargin,2)
    switch varargin{i}
      case 'origin'
        origin = varargin{i+1};
      case 'n'
        n = varargin{i+1};
      case 'mPP'
        mPP = varargin{i+1};
      case 'mPD'
        mPD = varargin{i+1};
      case 'mPW'
        mPW = varargin{i+1};
      case 's'
        s = 1;
      otherwise
        file = varargin{1};
    end
  end
else
end

switch origin
    case 'OE'
        %[data, timestamps, ~] = load_open_ephys_data(file); % Load .continuous files
        switch file
            case 'ADC1'
                data = info.Data.Data(1).mapped(33,1:end); % Ch33 is ADC1
                data = double(data)'*0.195; %intan
                timestamps = double(info.Timestamps)/30000; %in seconds
            case 'ADC2'
                data = info.Data.Data(1).mapped(35,1:end); % Ch34 is ADC2
                data = double(data)'*0.195; %intan
                timestamps = double(info.Timestamps)/30000; %in seconds
        end
        data2 = normalize(data,'range',[0 1]);
        data2 = data2>mPP;

        [~,locs] = findpeaks(double(data2),timestamps,'MinPeakProminence'...
            ,mPP,'MinPeakDistance',mPD,'MinPeakWidth',mPW,'NPeaks',n);
        TTL = locs;

        if s == 1
        save('TTL.mat','TTL');
        elseif s ~= 1
        end
        
    case 'IN'
        TTL_id = input('TTL? ENTER NAME (of the variable containing the TTL):','s');
        data = load(file,TTL_id);
        feval(@()assignin('caller','data',  data.(TTL_id)));
        t_dig = 't_dig';
        timestamps = load(file,'t_dig');
        feval(@()assignin('caller','timestamps',  timestamps.(t_dig)));
        data2 = normalize(data,'range',[0 1]);
        data2 = data2>mPP;
        
        [~,locs] = findpeaks(double(data2),timestamps,'MinPeakProminence',mPP,...
            'MinPeakDistance',mPD,'MinPeakWidth',mPW,'NPeaks',n);
        TTL = locs;
        
        if s == 1
        save('TTL.mat','TTL');
        elseif s ~= 1
        end
        
    case 'WS'
        load(file);
        data = input('VARIABLE TO USE AS DATA:'); % e. g. If you tipe audio, 
        % it looks for a variable in the WS with this name to use as your data.
        timestamps = input('VARIABLE TO USE AS TIME:'); %e. g. If you tipe time, 
        % it looks for a variable in the WS with this name to use as your time.
        data2 = normalize(data,'range',[0 1]);
        % This part turns NaN into 0 (min value)
        for ii = 1:length(data2)
            if isnan(data2(ii)) == 1
                data2(ii) = 0;
            elseif  isnan(data2(ii)) == 0
                data2(ii) = data(ii);
            end
        end
        data2 = data2>mPP;
        
        [~,locs] = findpeaks(double(data2),timestamps,'MinPeakProminence'...
            ,mPP,'MinPeakDistance',mPD,'MinPeakWidth',mPW,'NPeaks',n);
        TTL = locs;
        
        if logical(max(strcmp(varargin, 'ttl_dur')) == 1)
            negadata = data .* -1;
            data3 = normalize(negadata,'range',[0 1]);
            for ii = 1:length(data3)
                if isnan(data3(ii)) == 1
                data3(ii) = 1;
                elseif  isnan(data3(ii)) == 0
                data3(ii) = negadata(ii);
                end
            end
            data3 = data3>-mPP;
            data3(end) = 0;

            [~,ttl] = findpeaks(double(data3),timestamps,'MinPeakProminence'...
                    ,mPP,'MinPeakDistance',mPD,'MinPeakWidth',mPW,'NPeaks',n);
             diff = ttl - TTL;
             diff = round(diff);
             
             TTL2 = nan(1,length(diff));
             for dd = 1:length(diff)-1
                if diff(dd) == ttl_dur
                    TTL2(dd) = TTL(dd);
                elseif  diff(dd) < ttl_dur
                    TTL2(dd) = ttl(dd) - ttl_dur;
                else
                    disp('Plot and check the TTL, there is something strange.');
                end
             end
             TTL2(end) = TTL(end);
             TTL = TTL2;
        end

        if s == 1
        save('TTL.mat','TTL');
        elseif s ~= 1
        end

end
end
