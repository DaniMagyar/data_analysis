function [SortIDX, testMean] = calculate_ZscoreOrder(PSTHall, Stim, preferences)

%% Calculating plotting order based on Z-score change
testWindow_firstBin = abs(preferences.int(1))*preferences.fs/preferences.psth_bin+1;
testWindow_lastBin = testWindow_firstBin + (preferences.testBins-1);    
testWindow = PSTHall(:, testWindow_firstBin:testWindow_lastBin);
testMean = mean(testWindow,2);
switch Stim
    case {'ChETA_50_20Hz', 'BA_25_5Hz','BA_25_10Hz', 'BA_250_5Hz', 'shock_only'}
        [~,SortIDX] = sort(testMean, 'descend');

    case {'TO_25_5Hz', 'TO_25_10Hz', 'TO_250_5Hz', 'shock_inh'}
        load ([preferences.mainFolder '\BAparams.mat'], 'SortIDX')
end

