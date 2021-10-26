
function [PCA, Dendrogram] = PCA_analysis_psth_inter(File, nucleus, Stim, varargin)

% Function PCA_analysis_psth:
% Get data from psth of a neuron and perform a PCA analysis.
% 
% INPUTS: 
%   - File: TABLE. List of the neurons.
%   - nucleus: STRING. Nucleus from which you want to analyse the neurons.
%   More than one nucleus can be selected. Ex. {'LA','BA'}.
%   - TTL: STRING. Stimulus used as TTL.  IMPORTANT, the name of the ttl
%   should be saved in a variable called TTL.mat. Review section "Load
%   TTLs" to check the name of these variables.
% Varargin
%   - 'fs': sampling rate. Dafault 30000.
%   - 'interval': specify the time before and after the ttl to analyse.
%   Default [-1 1] in s. Format [-time time] in s.
%   - 'numPCA': number of principal components to use. Default 3. 
%   - 'numClusers': number of clusters. Default 4, assuming 1) Activated,
%   2) Inactivated, 3) Not affected, 4) Others.
%   - 'normality': normalise waveform data. Default 1 = Yes. No = 0.
%   Default method 'zscore', to change this go to section "Generate the
%   PSTHall matrix".
%   - 'method': creates the hierarchical tree using the specified method.
%   Default 'ward'. For other OPTIONS see mathworks documentation: LINKAGE
%   - 'metric': performs clustering using metric and computes the distance
%   between the PCA elements. Default 'euclidean'. For other OPTIONS see
%   mathworks documentation: LINKAGE
%   - 'bin': bining of the psth. Default 600, corresponding to 20ms.
%
% OUTPUTS: 
%   - PCA: principal component scores of the numPCA calculated.
%   - Dendrogram: tree containing hierarchical clusters of the PCA.
%   - a figure showing the PCA analysis.
%   
% Examples: 
% [PCA, Dendrogram] = PCA_analysis_psth(Neurons_List, 'LA', 'SK',); 
% [PCA, Dendrogram] = PCA_analysis_psth(Neurons_List, 'LA', 'SK','numPCA',3,'numClusers',4); 
% [PCA, Dendrogram] = PCA_analysis_psth(Neurons_List, 'LA', 'SK','method','complete','metric','mahalanobis'); 
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Instinute of Experimental Medicine, Hungary.
%
% MATLAB toolboxes: - Statistics and Machine Learning Toolbox 
% -------------------------------------------------------------------------

% Default params ----------------------------------------------------------
fs       = 30000; % Firing rate
int      = [-0.025 0.025]; % psth interval
PCA_num  = 2; % Number of principal components to use
clustnum = 4; % Number of clusters
norm     = 1; % Normalise data
LKmethod = 'complete';
LKmetric = 'mahalanobis';
psth_bin = 60; % 600 = 20ms
mainFolder = 'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR';

% User defined params -----------------------------------------------------
if nargin
    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'fs'
                fs       = varargin{ii+1};
            case 'interval'
                int      = varargin{ii+1};
            case 'numPCA'
                PCA_num  = varargin{ii+1};
            case 'numClusers'
                clustnum = varargin{ii+1};
            case 'normality'
                norm     = varargin{ii+1};
            case 'method'
                LKmethod = varargin{ii+1};
            case 'metric'
                LKmetric = varargin{ii+1};
            case 'bin'
                psth_bin = varargin{ii+1};
        end
    end
end

% Select neurons assigned to the selected NUCLEUS and Set the folders -----
Tab = File;
for ii = size(File,1):-1:1
    if startsWith(File.Location(ii), nucleus) == 0
        Tab(ii,:) = [];
    end
end


folder     = Tab.Recording;
nNeurons   = size(folder,1);

% Generate the PSTHall matrix ---------------------------------------------
if isempty(int) == 1
    int = [-1 1];
end

% PSTHall =
% zeros(size(Tab,1),ceil((fs*(abs(int(1))+abs(int(2)))+1)/psth_bin));
% %flexible
PSTHall = zeros(size(Tab,1),ceil((fs*(abs(int(1))+abs(int(2))))/psth_bin));
for kk = 1:nNeurons
    % Go to folder
    cd([mainFolder filesep folder{kk,:}]);
    
    % Load TTLs & Choose variable to use as stimulus
    if kk == 1
        load('TTLs.mat'); %#ok<LOAD>
        switch Stim
            case 'BA_25'
                ttl = BA_25;
            case 'TO_25'
                ttl = TO_25;
            case 'BA_250'
                ttl = BA_250;
            case 'TO_250'
                ttl = TO_250;
            case 'SK'
                ttl = TTL; % The SK TTL is already called TTL.
            case 'WN'
                if exist('NewWN','var')
                    ttl = NewWN;
                elseif exist('NewWS','var')
                    ttl = NewWS;
                elseif exist('TTL_WN','var')
                    ttl = TTL_WN;
                elseif exist('TTL_WN_video','var')
                    ttl = TTL_WN_video;
                elseif exist('sound','var')
                    ttl = sound;
                end
            case 'PT'
                if exist('NewPT','var')
                    ttl = NewPT;
                elseif exist('TTL_PT','var')
                    ttl = TTL_PT;
                elseif exist('TTL_PT_video','var')
                    ttl = TTL_PT_video;
                elseif exist('sound','var')
                    ttl = sound;
                end
        end
    elseif min(folder{kk,:} == folder{kk-1,:}) == 1
        % Use the same ttl
    elseif min(folder{kk,:} == folder{kk-1,:}) == 0
        clearvars ttl TTL NewWN NewWS TTL_WN TTL_WN_video NewPT TTL_PT...
            TTL_PT_video sound
        load('TTLs.mat'); %#ok<LOAD>
        switch Stim
            case 'BA_25'
                ttl = BA_25;
            case 'TO_25'
                ttl = TO_25;
            case 'BA_250'
                ttl = BA_250;
            case 'TO_250'
                ttl = TO_250;
            case 'SK'
                ttl = TTL; % The SK TTL is already called TTL.
            case 'WN'
                if exist('NewWN','var')
                    ttl = NewWN;
                elseif exist('NewWS','var')
                    ttl = NewWS;
                elseif exist('TTL_WN','var')
                    ttl = TTL_WN;
                elseif exist('TTL_WN_video','var')
                    ttl = TTL_WN_video;
                elseif exist('sound','var')
                    ttl = sound;
                end
            case 'PT'
                if exist('NewPT','var')
                    ttl = NewPT;
                elseif exist('TTL_PT','var')
                    ttl = TTL_PT;
                elseif exist('TTL_PT_video','var')
                    ttl = TTL_PT_video;
                elseif exist('sound','var')
                    ttl = sound;
                end
        end
    end
    
    % Select the neuron
    NeuronID = dir(['*',Tab.Group{kk},'_',num2str(Tab.Neuron(kk)),...
        '.mat']);
    NeuronID = NeuronID.name;
    
    load(NeuronID,'TS')
    TT = TS(:,1)/10000;
    
    psth1 = ttl_psth (TT*fs, ttl*fs, psth_bin, 'pre', abs(int(1)),...
        'post', abs(int(2)));
    [psth_spx, psth_t] = psth_hist(psth1, psth_bin);
    
    if norm == 1
%         BL = median(psth_spx(1:sum(double(psth_t < 0))));
%         newpsth = psth_spx(:) - BL;
%         newpsth = zscore(psth_spx); flexible
         newpsth = zscore(psth_spx(1:(end-1)));
    else
        newpsth = psth_spx;
    end
    PSTHall(kk,:) = newpsth';
end

time = psth_t;

% Colormap blue:white:red -------------------------------------------------
c1 = 1/255*[0,115,185];
c2 = 1/255*[239,239,239];
c3 = 1/255*[217,84,26];
mycolormap(1:20,1)  = linspace(c1(1),c2(1),20);
mycolormap(21:64,1) = linspace(c2(1),c3(1),44);
mycolormap(1:20,2)  = linspace(c1(2),c2(2),20);
mycolormap(21:64,2) = linspace(c2(2),c3(2),44);
mycolormap(1:20,3)  = linspace(c1(3),c2(3),20);
mycolormap(21:64,3) = linspace(c2(3),c3(3),44);

% PCA analysis ------------------------------------------------------------
% In "PSTHall" the rows are the neurons
[~,PCA1]  = pca(PSTHall); % running pricipal component analysis
PCA2      = PCA1(:,1:PCA_num); % extracting the first PCA_num pricinpal components 
inNAN = nanmedian(PCA2,1); % to avoid NaN values, set them to the median
for jj = 1:nNeurons
    if isnan(PCA2(jj,1)) == 1
        PCA2(jj,1:3) = inNAN;
    end
end
Dend      = linkage(PCA2,LKmethod,LKmetric); % calculating the dendrogram
Clusters  = cluster(Dend,'maxclust',clustnum); % clustering based on the tree
D         = pdist(PCA2); %euclidean distrance between point in the artifical space
leafOrder = optimalleaforder(Dend,D); % optimal order for plotting

PCA        = PCA2;
Dendrogram = Dend;

% PCA Figure --------------------------------------------------------------
figure; 
subplot(1,4,1:2)
imagesc(time/fs, 1:size(PSTHall,1), PSTHall(leafOrder,:)); %plotting sorted psth matrix 
clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
colormap(mycolormap); 
hold on;
plot([-0.00125 -0.00125],[1 size(PSTHall,1)],'r');
hold off;
ylabel('# Cell');
xlabel('Time (s)')

subplot(1,4,3:4)
cutoff = Dend(end-clustnum+2,3); % cutting the tree to get 'clustnum' clusters
h = dendrogram(Dend,0,'Orientation','right','reorder',leafOrder,...
    'ColorThreshold',cutoff); %plotting dendrogram
set(h,'LineWidth',1)
set(gca,'Ydir','reverse');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'visible','off')

end