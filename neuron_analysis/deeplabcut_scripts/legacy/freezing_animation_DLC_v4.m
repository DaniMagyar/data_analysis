function freezing_animation_DLC_v4

% Use input data generated with 'freezing_calculation_DLC_v2'

% work in CData

files = dir([cd '\*.avi']);
borderWidth = 20;
resizeRatio = 0.5;
TwinAll = load('TimeWindows.mat').TimeWindows;

for ii = 1:length(files)
    tic
    inputVideo = VideoReader(files(ii).name);
    outputVideo =  VideoWriter(['Analysis_' files(ii).name]);
    open(outputVideo)
    filename = dir(['*' files(ii).name(1:end-4) '*.mat']);
    loadBar = waitbar(0,['Processing video ' num2str(ii) ' ...']);
    load(filename.name);
    vidID = files(ii).name(1:5);
    vidType = files(ii).name(7:9);
    if vidType == 'rec'
        vidType = 'recall';
    end
    TwinCurr = TwinAll.(sprintf(vidID)).(sprintf(vidType));
    resizedVec(1:TwinCurr.startFrameID) = NaN;
    resizedVec(TwinCurr.endFrameID:end) = NaN;   

    for jj = 1:inputVideo.NumFrames
        
        vidFrame = readFrame(inputVideo, 'native');  %1 ms  
        vidFrame = imresize(vidFrame, resizeRatio);
        if ~rem(jj,1000)*jj/1000 >= 1
            waitbar((jj)/inputVideo.NumFrames,loadBar);
        end
        if resizedVec(jj) == 0   % Red image border upon freezing
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
        elseif resizedVec(jj) > g.acTh % Green image border upon activity
            vidFrame(1:borderWidth,:,1) = 0;
            vidFrame(1:borderWidth,:,2) = 128;
            vidFrame(1:borderWidth,:,3) = 0;
            vidFrame(end-borderWidth:end,:,1) = 0;
            vidFrame(end-borderWidth:end,:,2) = 128;
            vidFrame(end-borderWidth:end,:,3) = 0;
            vidFrame(:,1:borderWidth,1) = 0;
            vidFrame(:,1:borderWidth,2) = 128;
            vidFrame(:,1:borderWidth,3) = 0;
            vidFrame(:,end-borderWidth:end,1) = 0;
            vidFrame(:,end-borderWidth:end,2) = 128;
            vidFrame(:,end-borderWidth:end,3) = 0;     
        elseif isnan(resizedVec(jj))  % Black image border while no detection
            vidFrame(1:borderWidth,:,1) = 0;
            vidFrame(1:borderWidth,:,2) = 0;
            vidFrame(1:borderWidth,:,3) = 0;
            vidFrame(end-borderWidth:end,:,1) = 0;
            vidFrame(end-borderWidth:end,:,2) = 0;
            vidFrame(end-borderWidth:end,:,3) = 0;
            vidFrame(:,1:borderWidth,1) = 0;
            vidFrame(:,1:borderWidth,2) = 0;
            vidFrame(:,1:borderWidth,3) = 0;
            vidFrame(:,end-borderWidth:end,1) = 0;
            vidFrame(:,end-borderWidth:end,2) = 0;
            vidFrame(:,end-borderWidth:end,3) = 0;  
        else                                        % Yellow image border upon resting
            vidFrame(1:borderWidth,:,1) = 255;
            vidFrame(1:borderWidth,:,2) = 255;
            vidFrame(1:borderWidth,:,3) = 0;
            vidFrame(end-borderWidth:end,:,1) = 255;
            vidFrame(end-borderWidth:end,:,2) = 255;
            vidFrame(end-borderWidth:end,:,3) = 0;
            vidFrame(:,1:borderWidth,1) = 255;
            vidFrame(:,1:borderWidth,2) = 255;
            vidFrame(:,1:borderWidth,3) = 0;
            vidFrame(:,end-borderWidth:end,1) = 255;
            vidFrame(:,end-borderWidth:end,2) = 255;
            vidFrame(:,end-borderWidth:end,3) = 0;
        end    
        writeVideo(outputVideo,vidFrame) %2 ms
    end
    close(outputVideo)
    close(loadBar)
    toc

    pFr = numel(resizedVec(resizedVec==0));
    pAc = numel(resizedVec(resizedVec>g.acTh));
    pRes = numel(resizedVec)-(pFr+pAc);
    pie([pFr pAc pRes])
    labels = {'pFr','pAc','pRes'};
    legend(labels, 'Location', 'eastoutside')
    saveas(gcf,[vidID vidType '.jpg'])
    close
end
