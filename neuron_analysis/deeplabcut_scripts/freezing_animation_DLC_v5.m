function freezing_animation_DLC_v5

% Use input data generated with BehaviorDEPOT
% Video name should be e.g. MD219_recall   (.avi or .mp4)

% work in CData

files = dir([cd '\*.mp4']);
borderWidth = 20;
resizeRatio = 0.5;
TwinAll = load('TimeWindows.mat').TimeWindows;
TtrainingBase = 180; % baseline period before training

for ii = 1:length(files)
    tic
    inputVideo = VideoReader(files(ii).name);
    fps = inputVideo.FrameRate;
    outputVideo =  VideoWriter(['Analysis_' files(ii).name], 'MPEG-4');
    outputVideo.FrameRate = fps;
    open(outputVideo)
    loadBar = waitbar(0,['Processing video ' num2str(ii) ' ...']);
    load([cd '/' files(ii).name(1:end-4) '_analyzed/Behavior.mat'])
    vidID = files(ii).name(1:5);
    vidType = files(ii).name(7:9);
    if vidType == 'rec'
        vidType = 'recall';
    elseif vidType == 'tra'
        vidType = 'training';

    end
    TwinCurr = TwinAll.(sprintf(vidID)).(sprintf(vidType));
    resizedVec = Behavior.Freezing.Vector;
    resizedVec(1:TwinCurr.startFrameID) = NaN;
    resizedVec(TwinCurr.endFrameID:end) = NaN;   

    for jj = 1:inputVideo.NumFrames       
        vidFrame = readFrame(inputVideo, 'native');  %1 ms  
        vidFrame = imresize(vidFrame, resizeRatio);
        if ~rem(jj,1000)*jj/1000 >= 1
            waitbar((jj)/inputVideo.NumFrames,loadBar);
        end
        if resizedVec(jj) == 1   % Red image border upon freezing
            vidFrame(1:borderWidth,:,1) = 255;
            vidFrame(1:borderWidth,:,2) = 0;
            vidFrame(1:borderWidth,:,3) = 0;
            vidFrame(end-borderWidth:end,:,1) = 255;
            vidFrame(end-borderWidth:end,:,2) = 0;
            vidFrame(end-borderWidth:end,:,3) = 0;
            vidFrame(:,1:borderWidth,1) = 255;
            vidFrame(:,1:borderWidth,2) = 0;
            vidFrame(:,1:borderWidth,3) = 0;
            vidFrame(:,end-borderWidth:end,1) = 255;
            vidFrame(:,end-borderWidth:end,2) = 0;
            vidFrame(:,end-borderWidth:end,3) = 0;
        end    
        writeVideo(outputVideo,vidFrame) %2 ms
    end
    close(outputVideo)
    close(loadBar)
    toc

%     if vidType == 'training'
%             resizedVec(1:TwinCurr.startFrameID) = NaN;
%             resizedVec((TwinCurr.startFrameID+(TtrainingBase*fps)):end) = NaN;
%     end

    pFr = numel(resizedVec(resizedVec==1));
    pNFr = numel(resizedVec(~isnan(resizedVec))) - pFr; % non NaN frames - pFr
    pie([pFr pNFr])
    xlabel(['Time = ' num2str((pFr + pNFr)/fps) ' s']);
    labels = {'pFr','pTot'};
    legend(labels, 'Location', 'eastoutside')
    saveas(gcf,[vidID vidType '.jpg'])
    close
end
