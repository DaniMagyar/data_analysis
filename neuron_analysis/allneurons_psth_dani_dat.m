function psth = allneurons_psth_dani_dat(ttl, varargin)

% allneurons_psth_dani_dat(BA_50, 'pre', 0.1, 'post', 0.1, 'bin', 150, 's', 'figname', 'BA_50') 



% Get data from psth and save to mat for each neuron
% IMPORTANT: For the analysis of pairs use ttl_psth.
% 
% INPUTS: 
%   - ngroup: groups of neurons to perform the analysis. It can be a single
%   group (ngroup = 2, for GR2) or several (ngroup = [1:8], for GR1 to
%   GR8). 
%   - ttl: a vector containing the ttl timestamps.
% Varargin
%   - 'fs': sampling rate. Dafault 30000.
%   - 'Dt': Delay time, between the recording and the ttl. Ex. in case of
%   getting the ttl from a video. Default NOT used. In SECONDS!
%   - 'pre': time before trigger to include in psth. Default 1s. 
%   - 'post': time after trigger to include in psth. Default 1s.
%   - 'bin': width of the histogram column. Adjust the 'bin', 1000 bin for
%   0.5s; 3000 bin for 2s. 
%   - 'excel': to save an excel of the psth. Default NOT save excel 
%   - 's': to save the images. Default NOT save figs
%   - 'name': Name of Valiable to be saved. Default ['psth', num2str(floor(post/fs)),
%       's'].
% 
% OUTPUTS: 
%   - psth: a structure with the spikes around the event (all and up to
%   1st), and the mean values for the firing times of the neuron.
%   
% Examples: 
% psth = allneurons_psth([1:8], ttl); 
% psth = allneurons_psth(2, ttl, 'pre',0.5,'post',0.5,'fr',30000); 
% psth = allneurons_psth([1:8], sound(:,1), 'name', 'psth1s_sound', 'pre', 0.5, 'post', 1, 'Dt', -5,'s');
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2020
% Laboratory of Network Neurophysiology
% Instinute of Experimantal Medicine, Hungary.
%
% MATLAB toolboxes: - 
% -------------------------------------------------------------------------
  
% Default params ----------------------------------------------------------
fs = 30000;
pre = 0.5;
post = 0.5;
bin = 1000;
excel = 0; % Default NOT save excel
s = 1; % Default NOT save figs
name = ['psth_',num2str(floor(post/fs))]; % Name of Valiable to be saved
figname = '';
ngroup = (1:10);
% Params introduced in the varargin ---------------------------------------
if nargin
    for ii = 1:length(varargin)
        switch varargin{ii}
            case 'fs'
                fs = varargin{ii+1};
            case 'Dt'
                Dt = varargin{ii+1};
            case 'pre'
                pre = varargin{ii+1};
                if pre >= 5000
                elseif pre < 5000
                    pre = pre * fs;
                end
            case 'post'
                post = varargin{ii+1};
                if post >= 5000
                elseif post < 5000
                    post = post * fs;
                end
                if max(strcmp(varargin, 'name')) == 1
                elseif max(strcmp(varargin, 'name')) == 0
                    name = ['psth',num2str(floor(post/fs)),'s'];
                end
            case 'bin'
                bin = varargin{ii+1};
            case 'excel'
                excel = 1;                
            case 's'
                s = 1;
            case 'name'
                name = varargin{ii+1};
                B = strfind(name,'_');
                figname = name(B:end);
            case 'figname'
                figname = varargin{ii+1};
            otherwise
        end
    end
end

for nn = ngroup
% tmp = dir(['*',num2str(nn),'_','*.mat']); % all .mat that belong to the group
tmp = dir([ cd '\GR_35\' '*',num2str(nn),'_','*.mat']); % find the same .mat files, but in a subfolder, in order to access structure.oebin simultaneously
files = {tmp.name}'; 

% This reorders the files so after 1 comes 2 and not 10. Thank you MATLAB.
if length(files)>9
    A = length(files)-9;
    for kk = 1:A
        B = files{1+kk};
        [files{1+kk}] = [];
        files{end+1} = B;
    end
    
    for kk = 1:length(files)
        file(kk) = ~isempty(files{kk});
    end
    files = files(file==1);
end

% Preallocate vars --------------------------------------------------------
psth = struct;
xls = zeros(length(files),4); % To save an excel with relevant values
Gr_N = cell(length(files),1);
Time = zeros(length(files),1);

% Let's calculate things --------------------------------------------------
for ii = 1:length(files)
    
    FstSpk_means = zeros(1,4);
    AllSpk_means = [];
    AllSpk = [];
    FstSpk = [];

    neuron = files{ii};
    load([cd '\GR_35\' neuron],'TS');
    
    if exist('Dt', 'var') == 0 
        TT = ((TS(:,1)/10000)); 
    elseif exist('Dt', 'var') == 1
        if exist('ts','var') == 0
            A = dir('*.continuous');
            B = char({A.name});
            [~, ts, ~] = load_open_ephys_data(B(1,:));
%             TT = ((TS(:,1)/10000)-(min(ts)+Dt)/60); % Normal way
            TT = ((TS(:,1)/10000)-(min(ts))+Dt); % Unsusual way
        elseif exist('ts','var')== 1
%             TT = ((TS(:,1)/10000)-(min(ts)+Dt)/60); % Normal way
            TT = ((TS(:,1)/10000)-(min(ts))+Dt); % Unsusual way
        end
    end
    
    % Calculate spike time around ttl 
    [psth1,ts1,~,ts2] = ttl_psth_dani(TT*fs, ttl*fs, bin, 'fs', fs, 'pre',...
        pre/fs, 'post', post/fs, 'chart', 2);
   
    % All spikes before & after ttl            
    for jj = 1:length(ts1)
        A = cell2mat(ts1(jj,1));
        if isempty(A) == 1
            A = NaN;
        end
        n = max(size(AllSpk,1),size(A,1));
        AllSpk(end+1:n,:) = nan;
        A(end+1:n,1) = nan;
        AllSpk = [AllSpk, A];
    end 
    AllSpk(:,1) = [];
    AllSpk = AllSpk/fs;
    
    % Spikes before & 1st after ttl 
    for kk = 1:length(ts2)
        A = cell2mat(ts2(kk,1));
        if isempty(A) == 1
            A = NaN;
        end
        n = max(size(FstSpk,1),size(A,1));
        FstSpk(end+1:n,:) = nan;
        A(end+1:n,1) = nan;
        FstSpk = [FstSpk, A]; 
    end 
    FstSpk(:,1) = [];
    FstSpk = FstSpk/fs;
    
    % Find means of 1st spike after ttl
    FstSpk_means(1,1) = nanmean(FstSpk(FstSpk>0));
    FstSpk_means(1,2) = nanmedian(FstSpk(FstSpk>0));
    FstSpk_means(1,3) = nanstd(FstSpk(FstSpk>0));
    FstSpk_means(1,4) = length(FstSpk(FstSpk>0));
        
    % Calculate the firing peaks (and times) for all-spikes
    [psth_spx, psth_t] = psth_hist(psth1, bin);
    [AllSpk_means(:,1), AllSpk_means(:,2)] = findpeaks(psth_spx,psth_t/fs,...
            'MinPeakHeight',0.8);
    if isempty(AllSpk_means) == 1
        AllSpk_means(1,1) = NaN;
        AllSpk_means(1,2) = NaN;
    end
    
    psth = struct('FstSpk', FstSpk, 'AllSpk', AllSpk,...
        'FstSpk_means', FstSpk_means, 'AllSpk_means', AllSpk_means);
    eval([name '= psth']);
    save([cd '\GR_35\' neuron],name,'-append');
    
    if s == 1 % Save figure -----------------------------------------------
        saveas(gcf, matlab.lang.makeValidName([neuron(1:end-4),...
            '_', num2str(floor(post/fs)),'s' figname]), 'jpg')
%         saveas(gcf, matlab.lang.makeValidName([group(1:end-4),'_',...
%             num2str(ii),'_',num2str(floor(post/fs)),'s']), 'svg')
        close gcf
    else 
        close gcf
    end
    
    % Generate variable with the FstSpk_means to later save --------------- 
    Gr_N{ii,:}  = neuron(1:end-4);
    xls(ii,1:4) = FstSpk_means;
    Time(ii,:)  = floor(post/fs);
    clearvars -except GR Dt ttl files fs pre post bin ts TTL...
        name sound excel sec Gr_N s xls group psth Time figname
%     close gcf;
end

if excel == 1 % Save excel ------------------------------------------------
    
    % Save neurons FstSpk_means in the GR .mat file -----------------------
    PSTH_psth = table('Size',[length(files) 6],'VariableTypes',...
            {'string','double','double','double','double','double'},...
            'VariableNames',{'Neuron','T_psth','Mean','Median','Std','numSpk'});
    PSTH_psth.T_psth = Time;
    PSTH_psth.Neuron = Gr_N;
    PSTH_psth.Mean = xls(:,1);
    PSTH_psth.Median = xls(:,2);
    PSTH_psth.Std = xls(:,3);
    PSTH_psth.numSpk = xls(:,4);

    varInfo = who('-file',[Gr_N{1}(1,1:end-2), '.mat']);
    if ismember('PSTH_tab', varInfo) == 0
        save([Gr_N{1}(1,1:end-2), '.mat'], 'PSTH_tab', '-append');
    elseif ismember('PSTH_tab', varInfo) == 1
        PSTH_tab2 = PSTH_psth;
        load([Gr_N{1}(1,1:end-2), '.mat'], 'PSTH_tab');
        PSTH_psth = [PSTH_psth;PSTH_tab2];
        save([Gr_N{1}(1,1:end-2), '.mat'], 'PSTH_tab', '-append');
    end
    
    writetable(PSTH_psth, [Gr_N{1}(1,1:end-2),'_', num2str(floor(post/fs)),'.xls']);
else
end

end

end