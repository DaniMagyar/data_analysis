function freezing_calculation_DLC_v2

g.thFreezing = 2.5; % detection threshold between frames
g.thResting = 4; % detection threshold between frames
g.binsize = 30; % frames binned together (~60FPS, 16.6ms/frame)
g.acWeight = 5; % multiply activity value when large changes detected to separate resting from moving
g.acTh = 1; % activity threshold between resting and moving
g.frTh = 0.3;  % activity threshold between freezing and resting


g.files = dir([cd '\*.csv']);

for ii = 1:length(g.files)
    disp(g.files(ii).name)
    data = readcell(g.files(ii).name);
    bodyparts = unique(data(2,2:end));
    xvalues = cell2mat(data(4:end,2:3:end));
    yvalues = cell2mat(data(4:end,3:3:end));
    xDiff = [0; sum(abs(diff(xvalues)),2)]; % absolute movement of all bodyparts on on X axis, zero corrects length diff
    yDiff = [0; sum(abs(diff(yvalues)),2)]; % absolute movement of all bodyparts on on Y axis
    movement = sum(horzcat(xDiff,yDiff),2); % absolute movement of all bodyparts on both axes
    framediff = [0; abs(diff(movement))]; 

    activity(1:numel(framediff)) = 1; % create a vector, where everyting is resting
    activity(framediff<g.thFreezing) = 0; % values under thFreezing coded as 0 (freezing)
    activity(framediff>=g.thResting) = 2*g.acWeight; % values above thResting coded as 2*acWeight (activity)
    activity = [0, activity]; % increase length by 1
    
    RoundFreezing = activity;
    xRound = ceil(numel(RoundFreezing)/g.binsize)*g.binsize; 
    RoundFreezing(numel(RoundFreezing)+1:xRound) = 0; % length rounded up to next number divisible by 'binsize'

    
    freezingBinned =  mean(reshape(RoundFreezing,[g.binsize,numel(RoundFreezing)/g.binsize]));
    timestampsBinned = linspace(1,numel(RoundFreezing),numel(freezingBinned));
    for jj =  1:g.binsize
        resized(jj,:) = freezingBinned;
    end
    resizedVec =reshape(resized, [], 1)';

    plot(movement)
    hold on
    plotlim = ylim;
    b = bar(ones(numel(resizedVec),1)*plotlim(2), 'y');
    b.FaceAlpha = 0.2;
    b = bar((resizedVec==0)*plotlim(2), 'r');
    b.FaceAlpha = 0.2;
    b = bar((resizedVec>g.acTh)*plotlim(2), 'g');
    b.FaceAlpha = 0.2;
    hold off
    save(['Analysis_' g.files(ii).name(1:end-4) '.mat'], 'resizedVec', 'g')
    clearvars -except g ii
end