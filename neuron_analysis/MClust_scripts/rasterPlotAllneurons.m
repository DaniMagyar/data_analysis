function rasterPlotAllneurons(File, Stim)
% Generates a raster plot for the specified unit for a specific time
% window.
% File : MyNewOrder table from Z_score_analysis_Dani

% Default params ----------------------------------------------------------
mainFolder = 'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december\MClust';

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
            case 'shock_only'
                ttl = shock_only;
            case 'shock_inh'
                ttl = shock_inh;
        end

    % Select the neuron
    NeuronID = (['GR',Tab.Group{kk},'_',num2str(Tab.Neuron(kk))]);
    % Plot raster
    raster_plot(NeuronID, ttl, 'window', [-0.5 0.5])
%     switch Stim
%         case {'BA_25_5Hz', 'shock_only'}
%             saveas(gcf, [mainFolder '\rasters\MyNewOrder_' num2str(kk) '_' Stim '.png'])
%         case {'TO_25_5Hz', 'shock_inh'}
%             saveas(gcf, [mainFolder '\rasters\MyNewOrder_' num2str(kk) '_' Stim '.png'])
%     end
     close
end