function [iFR,tFR] = instFiringRate_allneurons(File, nucleus)

Tab = File;
for ii = size(File,1):-1:1
    if startsWith(File.Location(ii), nucleus) == 0
        Tab(ii,:) = [];
    end
end

mainFolder = 'L:\Magyar_Daniel\experiments\M2_shock_response';
folder     = Tab.Recording;
nNeurons   = size(folder,1);

for kk = [38,7,117,52,169,49,111,160,28]% 1:nNeurons
    % Go to folder
    cd([mainFolder filesep folder{kk,:}]);
     % Select the neuron
    NeuronID = dir(['*',Tab.Group{kk},'_',num2str(Tab.Neuron(kk)),...
        '.mat']);
    NeuronID = NeuronID.name;
    
    load(NeuronID,'TS')
    timestamps = TS(:,1)/10000;
    time = [timestamps(1) timestamps(end)];
    bins = 1;
    [iFR,tFR] = instFiringRate(timestamps,time,bins);
    plot(tFR, iFR)
end