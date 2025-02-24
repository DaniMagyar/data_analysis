clear all;
animalID = { 'MD243', 'MD250', 'MD251', 'MD252', 'MD253', 'MD254', 'MD268', 'MD277', 'MD278'};

secPre = 16;
secPost = 31;
zscore = 0;
stim = {'tone', 'noise'};


for aa = 1:numel(stim)



    for ii = 1:numel(animalID)
        epoch_averages.(animalID{ii}) = BAfc_plot_pupil_changes_group(animalID{ii}, secPre, secPost, stim{aa});
        fpsArray.(animalID{ii}) = load([animalID{ii} '_data.mat'], 'sampling_rate');
        stimStarts.(animalID{ii}) = round(secPre*fpsArray.(animalID{ii}).sampling_rate)+1;
    end
    
    for ii = 1:size(epoch_averages.(animalID{1}),1)
        epochData = [];
        for jj = 1:numel(animalID)
            currEpoch = epoch_averages.(animalID{jj})(ii,stimStarts.(animalID{jj})-1500:stimStarts.(animalID{jj})+2999);
            if anynan(currEpoch)
                nanIndices = find(isnan(currEpoch));
                %error(num2str(nanIndices))
            end
            if zscore == 1
                 epochData = [epochData; (currEpoch-mean(currEpoch(stimStarts.(animalID{jj})-500:stimStarts.(animalID{jj}))))/...
                     std(currEpoch(stimStarts.(animalID{jj})-500:stimStarts.(animalID{jj})))];
                %epochData = [epochData; currEpoch-mean(currEpoch(stimStarts.(animalID{jj})-500:stimStarts.(animalID{jj})))];
            elseif zscore == 0
                epochData = [epochData; currEpoch];
            end
        end
        epochMean.(stim{aa})(ii,:) = nanmean(epochData, 1);

    end
end




dataMatrix1 = epochMean.tone;  % Replace this with your 8x3000 matrix
dataMatrix2 = epochMean.noise;  % Replace this with your 8x3000 matrix
timeVector = linspace(-15, 30, 4500);  % Adjust the time range as needed

figure;
numPlots = size(dataMatrix1, 1);  % Number of plots (rows in the matrix)
for i = 1:numPlots     
    subplot(2, numPlots / 2, i);
    plot(timeVector, dataMatrix1(i, :));  % Plot the i-th row from dataMatrix1
    hold on
    plot(timeVector, dataMatrix2(i, :));  % Plot the i-th row from dataMatrix2
    hold off
    ylim([40 70])
end
