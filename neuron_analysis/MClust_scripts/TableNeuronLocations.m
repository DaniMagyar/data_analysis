% xlsx file must be in current folder
clear all;
name = 'PFC_Chrimson_stGtACR_2021_dec'; 
filename = ([name '.xlsx']);
mainFolder = cd;
opts = detectImportOptions(name);
opts = setvartype(opts, 'Group', 'char');
Tab = readtable(name, opts);
%Settings
findContactSites = 0;
globalLayerBorders = 0;
%% Finding contact sites and fast layer separation
if findContactSites == 1
    for ii = 1:length(Tab.Recording) % finding contact site of each neuron
        %shank#1
        if strcmp(Tab.Group{ii},'1')
            Tab.ContactSite{ii} = 0 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter1')
            Tab.ContactSite{ii} = 2 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'2')
            Tab.ContactSite{ii} = 4 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter2')
            Tab.ContactSite{ii} = 6 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'3')
            Tab.ContactSite{ii} = 8 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter3')
            Tab.ContactSite{ii} = 10 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'4')
            Tab.ContactSite{ii} = 12 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter4')
            Tab.ContactSite{ii} = 14 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'5')
            Tab.ContactSite{ii} = 16 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter5')
            Tab.ContactSite{ii} = 18 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'6')
            Tab.ContactSite{ii} = 20 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter6')
            Tab.ContactSite{ii} = 22 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'7')
            Tab.ContactSite{ii} = 24 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter7')
            Tab.ContactSite{ii} = 26 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'8')
            Tab.ContactSite{ii} = 28 + Tab.Ch_inGr_(ii);
        %shank#2
        elseif strcmp(Tab.Group{ii},'9')
            Tab.ContactSite{ii} = 0 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter9')
            Tab.ContactSite{ii} = 2 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'10')
            Tab.ContactSite{ii} = 4 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter10')
            Tab.ContactSite{ii} = 6 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'11')
            Tab.ContactSite{ii} = 8 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter11')
            Tab.ContactSite{ii} = 10 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'12')
            Tab.ContactSite{ii} = 12 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter12')
            Tab.ContactSite{ii} = 14 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'13')
            Tab.ContactSite{ii} = 16 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter13')
            Tab.ContactSite{ii} = 18 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'14')
            Tab.ContactSite{ii} = 20 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter14')
            Tab.ContactSite{ii} = 22 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'15')
            Tab.ContactSite{ii} = 24 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'inter15')
            Tab.ContactSite{ii} = 26 + Tab.Ch_inGr_(ii);
        elseif strcmp(Tab.Group{ii},'16')
            Tab.ContactSite{ii} = 28 + Tab.Ch_inGr_(ii);
        end
    end
end

if globalLayerBorders == 1 % finding estimated layer of each  neuron
    for jj = 1:length(Tab.Recording) 
        if Tab.ContactSite{jj} >= 5 && Tab.ContactSite{jj} <= 16
            Tab.Layer{jj} = 'L5b';
        elseif Tab.ContactSite{jj} >= 17 && Tab.ContactSite{jj} <= 22
            Tab.Layer{jj} = 'L5a';
        elseif Tab.ContactSite{jj} >= 23 && Tab.ContactSite{jj} <= 32
            Tab.Layer{jj} = 'L2/3';
        end
    end
else 
    Tab.Layer = Tab.ValidLayer;
end

%% Selection of included neurons

for kk = 1:length(Tab.Recording) % selecting included recordings
    if any(strcmp(Tab.Recording{kk}(1:5), {'MD098', 'MD099', 'MD101', 'MD108'}))
        Tab.Location{kk} = Tab.Layer{kk};
    else
        Tab.Location{kk} = 'excluded_recording';
    end
end
folder     = Tab.Recording;      % calculating AvgFR
nNeurons   = size(folder,1);
for ll = 1:nNeurons
    % Go to folder
    cd([mainFolder filesep folder{ll,:}]);
    NeuronID = dir(['GR',Tab.Group{ll},'_',num2str(Tab.Neuron(ll)),'.mat']);
    NeuronID = NeuronID.name;
    load(NeuronID,'TS')
    TT = TS(:,1)/10000;
    % Calculate avg firing rate
    firstSpk = TT(1);
    lastSpk = TT(end);
    firingTime = lastSpk - firstSpk;
    firingRate = numel(TT)/firingTime;
    Tab.AvgFR{ll} = firingRate;
end
for mm = 1:length(Tab.Recording) % excluding low firing-rate neurons
    if Tab.AvgFR{mm} < 0.1
        Tab.Location{mm} = 'excluded_FR';
    end
end
for nn = 1:length(Tab.Recording) % excluding low quality neurons
    if any(strcmp(Tab.Type{nn}, {'MUA', 'good'}))
        Tab.Location{nn} = 'excluded_quality';
    end
end
cd (mainFolder)
save ([name '.mat'], 'Tab')
clear all;