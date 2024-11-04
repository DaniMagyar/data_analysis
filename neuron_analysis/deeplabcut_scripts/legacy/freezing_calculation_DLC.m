function freezing_calculation_DLC

threshold = 5; % lower value requires more immobility
binsize = 10; % frames binned together (~60FPS, 16.6ms/frame)


files = dir([cd '\*.csv']);

for ii = 1:length(files)
    data = readcell(files(ii).name);
    bodyparts = unique(data(2,2:end));
    xvalues = cell2mat(data(4:end,2:3:end));
    yvalues = cell2mat(data(4:end,3:3:end));
    movement = sum(horzcat(sum(xvalues,2),sum(yvalues,2)),2);
    framediff = abs(diff(movement)); % length decreasing by 1
    freezing(framediff<threshold) = 1;
    freezing(framediff>=threshold) = 0;
    freezing = [0, freezing]; % increase length by 1
    
    RoundFreezing = freezing;
    xRound = ceil(numel(RoundFreezing)/binsize)*binsize; 
    RoundFreezing(numel(RoundFreezing)+1:xRound) = 0; % length rounded up to next number divisible by 'binsize'

    
    freezingBinned =  mean(reshape(RoundFreezing,[binsize,numel(RoundFreezing)/binsize]));
    timestampsBinned = linspace(1,numel(RoundFreezing),numel(freezingBinned));
    for jj =  1:binsize
        resized(jj,:) = freezingBinned;
    end
    resizedVec =reshape(resized, [], 1)';

    plot(movement)
    hold on
    plotlim = ylim;
    freezingBar(resizedVec==1) = plotlim(2);
    bar(freezingBar)
    
    save(['freezing_' files(ii).name(1:end-4) '.mat'], 'freezingBar')
end
    