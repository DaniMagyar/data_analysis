function fps_changer

% Use input data generated with BehaviorDEPOT

% work in CData

files = dir([cd '\*.avi']);

for ii = 1:length(files)
    inputVideo = VideoReader(files(ii).name);
    Table = readtable([files(ii).name(1:end-4) '.csv']);
    date = Table.Var17;   
    time = cellfun(@(fun) fun(12:end-6), date, 'UniformOutput', false); % @(fun) fun...: way to write a function  
    hours = str2num(cell2mat(cellfun(@(time) time(1:2), time, 'UniformOutput', false)));
    minutes = str2num(cell2mat(cellfun(@(time) time(4:5), time, 'UniformOutput', false)));
    seconds = str2num(cell2mat(cellfun(@(time) time(7:end), time, 'UniformOutput', false)));
    conv_time = hours*60*60 + minutes*60 + seconds;
    conv_diff = diff(conv_time);
    fps = numel(conv_time) / (conv_time(end) - conv_time(1));
    outputVideo =  VideoWriter(['FPScorr_' files(ii).name(1:end-4)], 'MPEG-4');
    outputVideo.FrameRate = fps;
    
    open(outputVideo)
    loadBar = waitbar(0,['Processing video ' num2str(ii) ' ...']);
    for jj = 1:inputVideo.NumFrames    
        tic
        vidFrame = readFrame(inputVideo);  %1 ms 
        if ~rem(jj,1000)*jj/1000 >= 1
            waitbar((jj)/inputVideo.NumFrames,loadBar);
        end 
        writeVideo(outputVideo,vidFrame) %2 ms
        toc
    end
    close(outputVideo)
    close(loadBar)
    clearvars -except files ii
end
