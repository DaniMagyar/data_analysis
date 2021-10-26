function [PCA, Dendrogram] = PCA_analysis_reworked(File, nucleus, Stim, varargin)
% Default params ----------------------------------------------------------
fs       = 30000; % Firing rate
int      = [-0.03 0.03]; % psth interval, even numbers are better for plot
PCA_num  = 2; % Number of principal components to use
clustnum = 3; % Number of clusters
norm     = 1; % Normalise data
average  = 0; % plotting averages
LKmethod = 'complete';
LKmetric = 'mahalanobis';
psth_bin = 90; % 600 = 20ms
mainFolder = 'C:\Users\dmagyar\Desktop\Chrimson_stGtACR';

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

PSTHall = zeros(size(Tab,1),ceil((fs*(abs(int(1))+abs(int(2))))/psth_bin));

for kk = 1:nNeurons
    % Go to folder
    cd([mainFolder filesep folder{kk,:}]);
    
    % Load TTLs & Choose variable to use as stimulus
    load('TTLs.mat'); %#ok<LOAD>
    switch Stim
        case 'TTL_500'
            ttl = TTL;
        case 'TTL_50'
            ttl = TTL(1:10:end);
        case 'BA_25'
            ttl = BA_25;
        case 'TO_25'
            ttl = TO_25;
        case 'BA_250'
            ttl = BA_250;
%             ttl = BA_250(1:10:end);
%             ttl = [ttl; BA_250(2:10:end)];
%             ttl = [ttl; BA_250(3:10:end)];
%             ttl = sort(ttl);
        case 'TO_250'
            ttl = TO_250;
        case 'SK'
            ttl = TTL; % The SK TTL is already called TTL.
    end
    % Select the neuron
    NeuronID = dir(['*',Tab.Group{kk},'_',num2str(Tab.Neuron(kk)),...
        '.mat']);
    NeuronID = NeuronID.name;
    
    load(NeuronID,'TS')
    TT = TS(:,1)/10000;
    
    spktimes = TT*fs;
    trigtimes = ttl*fs;
    bins = psth_bin;
    pre = abs(int(1))*fs;
    post = abs(int(2))*fs;
    

    % Pre-allocate for speed --------------------------------------------------
    psth = zeros(pre+post,2);
    psth (:,1) = (-1*pre:1:post-1); % time base
    % Construct psth & trialspx -----------------------------------------------
    ts = cell(numel(trigtimes),1);
    for ii = 1:numel(trigtimes)
      clear spikes
      spikes = spktimes - trigtimes(ii); % all spikes relative to current trigtime
      ts{ii} = round(spikes(spikes>=-pre & spikes<(post-1))); % spikes close to current trigtime
      psth(ts{ii}+pre+1,2) = psth(ts{ii}+pre+1,2)+1; 
    end

        psth1 = psth;
        [psth_spx, psth_t] = psth_hist(psth1, psth_bin);
        
        %remove artefact
        if size(PSTHall,2)== fs*(abs(int(1))+abs(int(2)))/psth_bin
            PSTHall(:,fs*(abs(int(1))+abs(int(2)))/psth_bin) = [];
        end
        
        psth_spx(abs(int(1)*fs/psth_bin)+1)=[]; % delete  to remove artefact
        
       
        
        if norm == 1
            newpsth = zscore(psth_spx);
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
    MyNewOrder = Tab.Recording(leafOrder);
    MyNewOrder(:,2) = Tab.Group(leafOrder);
    MyNewOrder(:,3) = num2cell(Tab.Neuron(leafOrder));
    MyNewOrder(:,4) = Tab.Type(leafOrder);
    
    save ('BAparams.mat', 'leafOrder', 'Dend', 'Clusters', 'MyNewOrder')
    

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
    
    
    % The first cluster is the one who is the 'deepest'
    if average == 1
        for ll = 1:clustnum
            Clustmeans(ll,:) = mean(PSTHall((find(Clusters==ll)),:),1);
        end
    % 
    %     averageFig = tiledlayout (clustnum,1);
    %     for mm = 1:clustnum
    %         nexttile
    %         plot(time/fs, Clustmeans(manualclusterorder(mm),:))
    %     end
        for mm = 1:clustnum
            figure(mm+1)
            hold on
            plot(time/fs, Clustmeans(mm,:), 'LineWidth', 1.5)
            title(['n neurons = ' num2str(numel(find(Clusters==mm)))])
            hold off
        end
    end

end



