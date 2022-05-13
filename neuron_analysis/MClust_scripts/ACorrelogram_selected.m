function ACorrelogram_selected(File)

Tab = File;
mainFolder = 'L:\Magyar_Daniel\experiments\M2_shock_response';
saveFolder = 'L:\Magyar_Daniel\experiments\M2_shock_response\Autocorr\';
folder     = Tab.Recording;
nNeurons   = size(folder,1);

for kk = 1:nNeurons
    % Go to folder
    cd([mainFolder filesep folder{kk,:}]);
    
    NeuronID = dir(['*',Tab.Group{kk},'_',num2str(Tab.Neuron(kk)),...
        '.mat']);
    NeuronID = NeuronID.name;
    
    ACorrelogram(NeuronID)
    saveas(gcf,[saveFolder num2str(kk) '_' num2str(folder{kk})...
        '_' NeuronID(1:end-4)], 'jpg')
end