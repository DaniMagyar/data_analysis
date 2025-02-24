% Function to plot the selected animals' data
    function plotSelectedAnimals(~, ~)
        selectedAnimals = animalID([checkboxes.Value] == 1);  % Get selected animals
        if isempty(selectedAnimals)
            disp('No animals selected.');
            return;
        end
        
        % Load the data for the selected animals
        epoch_averages = struct();
        fpsArray = struct();
        stimStarts = struct();
        for aa = 1:numel(stim)
            for ii = 1:numel(selectedAnimals)
                epoch_averages.(selectedAnimals{ii}) = BAfc_plot_pupil_changes_group(selectedAnimals{ii}, secPre, secPost, stim{aa});
                fpsArray.(selectedAnimals{ii}) = load([selectedAnimals{ii} '_data.mat'], 'sampling_rate');
                stimStarts.(selectedAnimals{ii}) = round(secPre*fpsArray.(selectedAnimals{ii}).sampling_rate)+1;
            end

            for ii = 1:size(epoch_averages.(selectedAnimals{1}), 1)
                epochData = [];
                for jj = 1:numel(selectedAnimals)
                    currEpoch = epoch_averages.(selectedAnimals{jj})(ii, stimStarts.(selectedAnimals{jj})-2000:stimStarts.(selectedAnimals{jj})+3999);
                    if anynan(currEpoch)
                        nanIndices = find(isnan(currEpoch));
                    end
                    if zscore == 1
                        epochData = [epochData; currEpoch - mean(currEpoch(stimStarts.(selectedAnimals{jj})-2000:stimStarts.(selectedAnimals{jj})))];
                    elseif zscore == 0
                        epochData = [epochData; currEpoch];
                    end
                end
                epochMean.(stim{aa})(ii, :) = nanmean(epochData, 1);
            end
        end

        % Plot the data for the selected animals
        dataMatrix1 = epochMean.tone;  % 8x3000 matrix for 'tone'
        dataMatrix2 = epochMean.noise;  % 8x3000 matrix for 'noise'
        timeVector = linspace(-20, 40, 6000);  % Adjust time range as needed

        figure;
        numPlots = numel(selectedAnimals);  % Number of plots based on selected animals
        for i = 1:numPlots
            subplot(2, numPlots / 2, i);
            % Plot tone and noise data for the selected animals
            idx = find(strcmp(animalID, selectedAnimals{i}));
            plot(timeVector, dataMatrix1(idx, :));  % Plot tone data for the animal
            hold on;
            plot(timeVector, dataMatrix2(idx, :));  % Plot noise data for the animal
            hold off;
            title(selectedAnimals{i});
            xlabel('Time (s)');
            ylabel('Response');
        end
    end