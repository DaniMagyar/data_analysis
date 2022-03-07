function rasterPlotAllneurons(File, Stim)
% Generates a raster plot for the specified unit for a specific time
% window.
% File : MyNewOrder table from Z_score_analysis_Dani

% Default params ----------------------------------------------------------
mainFolder = 'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december';

Tab = File;
folder     = Tab.Recording;
nNeurons   = size(folder,1);

for kk = 1:nNeurons
    % Go to folder
    cd([mainFolder filesep folder{kk,:}]);
    
    % Load TTLs & Choose variable to use as stimulus
    load('TTLs.mat'); %#ok<LOAD>
        switch Stim
            case 'BA_25_5Hz'
                ttl = BA_25_5Hz;
            case 'TO_25_5Hz'
                ttl = TO_25_5Hz;
        end

    % Select the neuron
    NeuronID = (['GR',Tab.Group{kk},'_',num2str(Tab.Neuron(kk))]);
    % Plot raster
    raster_plot(NeuronID, ttl, 'window', [-5 5])
    if Stim(1:2) == 'BA'
        saveas(gcf, [mainFolder '\rasters\MyNewOrder_' num2str(kk) 'BA.png'])
    elseif Stim(1:2) == 'TO'
        saveas(gcf, [mainFolder '\rasters\MyNewOrder_' num2str(kk) 'TO.png'])
    end
    close
end