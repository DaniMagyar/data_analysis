function freezing_animation_DLC

files = dir([cd '\*.avi']);

for ii = 1:length(files)
    inputVideo = VideoReader(files(ii).name);
    outputVideo =  VideoWriter(['Freezing_' files(ii).name]);
    open(outputVideo)
    filename = dir(['*' files(ii).name(1:end-4) '*.mat']);
    resizedVec = load(filename.name);
    for jj = 1:inputVideo.NumFrames
        vidFrame = readFrame(inputVideo, 'native');  %1 ms
        imshow(vidFrame) %  imshow ~ 10ms, image ~3 ms
        hold on
        xlimits = xlim;
        ylimits = ylim;
        if resizedVec.resizedVec(jj) == 0       
            plot(xlimits(1)*ones(inputVideo.Width,1),linspace(ylimits(1), ylimits(2), inputVideo.Width), 'LineWidth', 25, 'Color', 'r') % left vertical
            plot(xlimits(2)*ones(inputVideo.Width,1),linspace(ylimits(1), ylimits(2), inputVideo.Width), 'LineWidth', 25, 'Color', 'r') % right vertical
            plot(linspace(xlimits(1), xlimits(2), inputVideo.Height), ylimits(1)*ones(inputVideo.Height,1), 'LineWidth', 25, 'Color', 'r')
            plot(linspace(xlimits(1), xlimits(2), inputVideo.Height), ylimits(2)*ones(inputVideo.Height,1), 'LineWidth', 25, 'Color', 'r')
        elseif resizedVec.resizedVec(jj) > resizedVec.g.acTh
            plot(xlimits(1)*ones(inputVideo.Width,1),linspace(ylimits(1), ylimits(2), inputVideo.Width), 'LineWidth', 25, 'Color', 'g') % left vertical
            plot(xlimits(2)*ones(inputVideo.Width,1),linspace(ylimits(1), ylimits(2), inputVideo.Width), 'LineWidth', 25, 'Color', 'g') % right vertical
            plot(linspace(xlimits(1), xlimits(2), inputVideo.Height), ylimits(1)*ones(inputVideo.Height,1), 'LineWidth', 25, 'Color', 'g')
            plot(linspace(xlimits(1), xlimits(2), inputVideo.Height), ylimits(2)*ones(inputVideo.Height,1), 'LineWidth', 25, 'Color', 'g')
        else 
            plot(xlimits(1)*ones(inputVideo.Width,1),linspace(ylimits(1), ylimits(2), inputVideo.Width), 'LineWidth', 25, 'Color', 'y') % left vertical
            plot(xlimits(2)*ones(inputVideo.Width,1),linspace(ylimits(1), ylimits(2), inputVideo.Width), 'LineWidth', 25, 'Color', 'y') % right vertical
            plot(linspace(xlimits(1), xlimits(2), inputVideo.Height), ylimits(1)*ones(inputVideo.Height,1), 'LineWidth', 25, 'Color', 'y')
            plot(linspace(xlimits(1), xlimits(2), inputVideo.Height), ylimits(2)*ones(inputVideo.Height,1), 'LineWidth', 25, 'Color', 'y')
        end
        hold off
        frame = getframe; %25ms
        writeVideo(outputVideo,frame) %2 ms
        disp(jj) 
    end
    close(outputVideo)
end
