function VideoAnalysisNeurotar_v2(varargin)


% Version 2: when manipulation and frame-TS match TTLs are different.

% Examples:
% -VideoAnalysisNeurotar_v2( 'recording', 'MD098', 'TTLmatch', 'ttl_match', 'TTLtest', 'BA_25_5Hz', 'FPS', 12.5)
% -VideoAnalysisNeurotar_v2( 'recording', 'MD099', 'TTLmatch', 'ttl_match', 'TTLtest', 'BA_25_5Hz', 'TTLtest_recorded', (1:25), 'FPS', 7.143)
% -VideoAnalysisNeurotar_v2( 'recording', 'MD101', 'TTLmatch', 'ttl_match', 'TTLtest', 'BA_25_5Hz','TTLtest_recorded', (1:25), 'ROIdetect_excluded', (1:200), 'FPS', 20)

prs =  inputParser;
addOptional(prs, 'FPS', '',@isnumeric) % camera framerate
addOptional(prs, 'recording', '',@ischar) % ID of recording, e.g. MD130
addOptional(prs, 'TTLmatch', '',@ischar) % variable name from TTLsKS.mat, timepoints for frame-TS matching
addOptional(prs, 'TTLtest', '',@ischar) % variable name from TTLsKS.mat, investigated timepoints
addOptional(prs, 'TTLtest_recorded', 0,@isnumeric) % e.g. TTLs recorded on video (e.g. 5:20) (This is required when first ttls were not recorded)
addOptional(prs, 'ROIdetect_refractory', 3,@isnumeric) % after ROI activation, stop detection for x seconds. (trainstart)
addOptional(prs, 'ROIdetect_excluded', 0, @isnumeric) % e.g. exclude contaminated frames from laser ROI. e.g (1:200)
addOptional(prs, 'pre_time', 2,@isnumeric)
addOptional(prs, 'post_time', 2,@isnumeric)
parse(prs,varargin{:})
g = prs.Results;

Table = readtable(append(g.recording, '_video_analysis'));
if g.ROIdetect_excluded ~= 0 % remove artefact
    Table.LaserON(g.ROIdetect_excluded) = {'False'};
end
TsMatchROI_Framestate = double(ismember(Table.LaserON, 'True'));
TsMatchROI_ON_idx = find(TsMatchROI_Framestate == 1);
TsMatchROI_ON_idx_2 = TsMatchROI_ON_idx;
TTL = load('TTLsKS.mat', g.TTLmatch, g.TTLtest);

TTLmatch = TTL.(g.TTLmatch);
TTLtest = TTL.(g.TTLtest);
if g.TTLtest_recorded ~= 0
    TTLtest = TTLtest(g.TTLtest_recorded);
    TTLmatch = TTLmatch(g.TTLtest_recorded);
end

disp(['TsMatchROI detected ' num2str(numel(TsMatchROI_ON_idx)) ' signals'])
disp(['Real TTL Timestamps used for matching: ' num2str(numel(TTLmatch))])

%% Calculate timestamps of each frame

for ii = 2:numel(TsMatchROI_ON_idx)
    if TsMatchROI_ON_idx(ii) < (TsMatchROI_ON_idx(ii-1) + (g.FPS * g.ROIdetect_refractory)) % all TRUE frames
        TsMatchROI_ON_idx_2(ii) = NaN;
    end
end
TsMatchROI_ON_idx_2(find(isnan(TsMatchROI_ON_idx_2))) = []; % first TRUE frames of each trial
disp(['Video TTL Timestamps used for matching: ' num2str(numel(TsMatchROI_ON_idx_2))])

if numel(TTLmatch) ~= numel(TsMatchROI_ON_idx_2)
    if g.TTLtest_recorded ~= 0
        TsMatchROI_ON_idx_2 = TsMatchROI_ON_idx_2(g.TTLtest_recorded);
    else
        error('Different stimulus number or far too high FPS')
    end
end

Table.Timestamp(TsMatchROI_ON_idx_2(1)) = TTLmatch(1);
Table.Timestamp(1) =  TTLmatch(1) - (TsMatchROI_ON_idx_2(1)-1)*(1/g.FPS);
Table.Timestamp(end) =  TTLmatch(1) + (numel(TsMatchROI_Framestate) - TsMatchROI_ON_idx_2(1)) * (1/g.FPS);
Table.Timestamp(1:end) = linspace(Table.Timestamp(1), Table.Timestamp(end), numel(Table.Timestamp));
MismatchAllTTL =  TTLmatch - Table.Timestamp(TsMatchROI_ON_idx_2);
disp(MismatchAllTTL)
if max(MismatchAllTTL) > 0.1 || min(MismatchAllTTL) < -0.1
    error('Incorrect FPS')
end


%% Detect movement
MovementROI = double(ismember(Table.Movement, 'True'));
for jj = 1:numel(TTLtest)
    preStim{:,jj} = MovementROI(Table.Timestamp>=(TTLtest(jj)-g.pre_time) & Table.Timestamp < TTLtest(jj));
    postStim{:,jj} = MovementROI(Table.Timestamp>=TTLtest(jj) & Table.Timestamp < (TTLtest(jj)+g.post_time));
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

switch g.TTLtest
    case 'BA_25_all'
        BA_25_all_Running = TTLtest(find(strcmp(MovementChange, 'Running')));
        BA_25_all_Resting = TTLtest(find(strcmp(MovementChange, 'Resting')));
        save('TTLsKS.mat', 'BA_25_all_Running', 'BA_25_all_Resting', '-append')
    case 'BA_25_5Hz'
        BA_25_5Hz_Running = TTLtest(find(strcmp(MovementChange, 'Running')));
        BA_25_5Hz_Resting = TTLtest(find(strcmp(MovementChange, 'Resting')));
        save('TTLsKS.mat', 'BA_25_5Hz_Running', 'BA_25_5Hz_Resting', '-append')
    case 'TO_25_5Hz'
        TO_25_5Hz_Running = TTLtest(find(strcmp(MovementChange, 'Running')));
        TO_25_5Hz_Resting = TTLtest(find(strcmp(MovementChange, 'Resting')));
        save('TTLsKS.mat', 'TO_25_5Hz_Running', 'TO_25_5Hz_Resting', '-append')
    case 'BA_50_5Hz'
        BA_50_5Hz_Running = TTLtest(find(strcmp(MovementChange, 'Running')));
        BA_50_5Hz_Resting = TTLtest(find(strcmp(MovementChange, 'Resting')));
        save('TTLsKS.mat', 'BA_50_5Hz_Running', 'BA_50_5Hz_Resting', '-append')
    case 'BA_50_10Hz'
        BA_50_10Hz_Running = TTLtest(find(strcmp(MovementChange, 'Running')));
        BA_50_10Hz_Resting = TTLtest(find(strcmp(MovementChange, 'Resting')));
        save('TTLsKS.mat', 'BA_50_10Hz_Running', 'BA_50_10Hz_Resting', '-append')
    case 'ChETA_50_20Hz'
        ChETA_50_20Hz_Running = TTLtest(find(strcmp(MovementChange, 'Running')));
        ChETA_50_20Hz_Resting = TTLtest(find(strcmp(MovementChange, 'Resting')));
        save('TTLsKS.mat', 'ChETA_50_20Hz_Running', 'ChETA_50_20Hz_Resting', '-append')
    case 'shock_only'
        shock_only_Running = TTLtest(find(strcmp(MovementChange, 'Running')));
        shock_only_Resting = TTLtest(find(strcmp(MovementChange, 'Resting')));
        save('TTLsKS.mat', 'shock_only_Running', 'shock_only_Resting', '-append')
    case 'shock_inh'
        shock_inh_Running = TTLtest(find(strcmp(MovementChange, 'Running')));
        shock_inh_Resting = TTLtest(find(strcmp(MovementChange, 'Resting')));
        save('TTLsKS.mat', 'shock_inh_Running', 'shock_inh_Resting', '-append')
end






















