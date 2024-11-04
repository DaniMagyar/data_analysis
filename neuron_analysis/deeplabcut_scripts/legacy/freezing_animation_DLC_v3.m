function freezing_animation_DLC_v3

% work in CData

files = dir([cd '\*.avi']);
borderWidth = 20;

for ii = 1:length(files)
    inputVideo = VideoReader(files(ii).name);
    outputVideo =  VideoWriter(['Freezing_' files(ii).name]);
    open(outputVideo)
    filename = dir(['*' files(ii).name(1:end-4) '*.mat']);
    resizedVec = load(filename.name);
    for jj = 1:inputVideo.NumFrames
        vidFrame = readFrame(inputVideo, 'native');  %1 ms  
        vidFrame = imresize(vidFrame, 0.5);

        if resizedVec.resizedVec(jj) == 0   % Red image border upon freezing
            vidFrame(1:borderWidth,:,1) = 255;
            vidFrame(1:borderWidth,:,2) = 0;
            vidFrame(1:borderWidth,:,3) = 0;
            vidFrame(581:600,:,1) = 255;
            vidFrame(581:600,:,2) = 0;
            vidFrame(5811:600,:,3) = 0;
            vidFrame(:,1:20,1) = 255;
            vidFrame(:,1:20,2) = 0;
            vidFrame(:,1:20,3) = 0;
            vidFrame(:,781:800,1) = 255;
            vidFrame(:,781:800,2) = 0;
            vidFrame(:,781:800,3) = 0;
        elseif resizedVec.resizedVec(jj) > resizedVec.g.acTh % Green image border upon activity
            vidFrame(1:20,:,1) = 0;
            vidFrame(1:20,:,2) = 128;
            vidFrame(1:20,:,3) = 0;
            vidFrame(581:600,:,1) = 0;
            vidFrame(581:600,:,2) = 128;
            vidFrame(5811:600,:,3) = 0;
            vidFrame(:,1:20,1) = 0;
            vidFrame(:,1:20,2) = 128;
            vidFrame(:,1:20,3) = 0;
            vidFrame(:,781:800,1) = 0;
            vidFrame(:,781:800,2) = 128;
            vidFrame(:,781:800,3) = 0;            
        else                                        % Yellow image border upon resting
            vidFrame(1:20,:,1) = 255;
            vidFrame(1:20,:,2) = 255;
            vidFrame(1:20,:,3) = 0;
            vidFrame(581:600,:,1) = 255;
            vidFrame(581:600,:,2) = 255;
            vidFrame(5811:600,:,3) = 0;
            vidFrame(:,1:20,1) = 255;
            vidFrame(:,1:20,2) = 255;
            vidFrame(:,1:20,3) = 0;
            vidFrame(:,781:800,1) = 255;
            vidFrame(:,781:800,2) = 255;
            vidFrame(:,781:800,3) = 0;
        end    
        tic
        writeVideo(outputVideo,vidFrame) %2 ms
        toc
    end
    close(outputVideo)
    close hFig;
end
