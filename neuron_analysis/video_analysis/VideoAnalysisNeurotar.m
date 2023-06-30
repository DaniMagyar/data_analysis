function VideoAnalysisNeurotar(varargin)

prs =  inputParser;
addOptional(prs, 'FPS', '',@isnumeric) % camera framerate
addOptional(prs, 'recording', '',@ischar) % ID of recording, e.g. MD130
addOptional(prs, 'TTL', '',@ischar) % variable name from TTLsKS.mat
addOptional(prs, 'TTLselected', 0,@isnumeric) % e.g. if first 5 TTL were not recorded on video
addOptional(prs, 'ROIdetect_refractory', 5,@isnumeric) % after ROI activation, stop detection for x seconds. (trainstart)
addOptional(prs, 'pre_time', 2,@isnumeric)
addOptional(prs, 'post_time', 2,@isnumeric)
parse(prs,varargin{:})
g = prs.Results;

Table = readtable(append(g.recording, '_video_analysis'));
TsMatchROI_Framestate = double(ismember(Table.LaserON, 'True'));
TsMatchROI_ON_idx = find(TsMatchROI_Framestate == 1);
TsMatchROI_ON_idx_2 = TsMatchROI_ON_idx;

TTL = load('TTLsKS.mat', g.TTL);
TTL = TTL.(g.TTL);
if g.TTLselected ~= 0
    TTL = TTL(g.TTLselected);
end

disp(['TsMatchROI detected ' num2str(numel(TsMatchROI_ON_idx)) ' signals'])
disp(['Real TTL Timestamps used for matching: ' num2str(numel(TTL))])

%% Calculate timestamps of each frame

for ii = 2:numel(TsMatchROI_ON_idx)
    if TsMatchROI_ON_idx(ii) < (TsMatchROI_ON_idx(ii-1) + (g.FPS * g.ROIdetect_refractory))
        TsMatchROI_ON_idx_2(ii) = NaN;
    end
end
TsMatchROI_ON_idx_2(find(isnan(TsMatchROI_ON_idx_2))) = [];

Table.Timestamp(TsMatchROI_ON_idx_2(1)) = TTL(1);
Table.Timestamp(1) = TTL(1) - (TsMatchROI_ON_idx_2(1)-1)*(1/g.FPS);
Table.Timestamp(end) = TTL(1) + (numel(TsMatchROI_Framestate) - TsMatchROI_ON_idx_2(1)) * (1/g.FPS);
Table.Timestamp(1:end) = linspace(Table.Timestamp(1), Table.Timestamp(end), numel(Table.Timestamp));
TTL_video = Table.Timestamp(TsMatchROI_ON_idx_2);
MismatchAllTTL = TTL - TTL_video;
disp(MismatchAllTTL)
if max(MismatchAllTTL) > 0.1 || min(MismatchAllTTL) < -0.1
    error('Incorrect FPS')
end


%% Detect movement
MovementROI = double(ismember(Table.Movement, 'True'));
for jj = 1:numel(TTL_video)
    preStim{:,jj} = MovementROI(Table.Timestamp>=(TTL_video(jj)-g.pre_time) & Table.Timestamp < TTL_video(jj));
    postStim{:,jj} = MovementROI(Table.Timestamp>=TTL_video(jj) & Table.Timestamp < (TTL_video(jj)+g.post_time));
end

preAvg = cellfun(@mean,preStim)';
postAvg = cellfun(@mean, postStim)';

labels = cellfun(@num2str,num2cell(1:numel(preAvg)),'UniformOutput',false);
for kk = 1: length(preAvg)
    x = [0.005*kk,4+0.005*kk];
    y = [preAvg(kk), postAvg(kk)];
    plot(x, y, 'rs-',  'LineWidth', 2, 'MarkerSize', 1);
    text(x, y, labels(kk))
    hold on
end
grid on

for ll = 1:length(preAvg)
    if preAvg(ll) >= 0.8
        if postAvg(ll) >= 0.8
            MovementChange{ll} = 'Running';
        elseif postAvg(ll) <= 0.2
            MovementChange{ll} = 'SlowDown';
        else MovementChange{ll} = 'unknown';
        end
    elseif preAvg(ll) <= 0.2
        if postAvg(ll) <= 0.2
            MovementChange{ll} = 'Resting';
        elseif postAvg(ll) >= 0.8
            MovementChange{ll} = 'SpeedUp';
        else MovementChange{ll} = 'unknown';
        end
    else MovementChange{ll} = 'unknown';
    end
end
MovementChange = MovementChange';

% switch g.TTL
%     case 'BA_50_5Hz'
%         BA_50_5Hz_Running = TTL(find(strcmp(MovementChange, 'Running')));
%         BA_50_5Hz_Resting = TTL(find(strcmp(MovementChange, 'Resting')));
%         save('TTLsKS.mat', 'BA_50_5Hz_Running', 'BA_50_5Hz_Resting', '-append')
%     case 'ChETA_50_20Hz'
%         ChETA_50_20Hz_Running = TTL(find(strcmp(MovementChange, 'Running')));
%         ChETA_50_20Hz_Resting = TTL(find(strcmp(MovementChange, 'Resting')));
%         save('TTLsKS.mat', 'ChETA_50_20Hz_Running', 'ChETA_50_20Hz_Resting', '-append')
% end
























