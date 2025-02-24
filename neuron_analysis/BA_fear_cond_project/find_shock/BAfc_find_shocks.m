function [psth_spx, idx] =  BAfc_find_shocks(varargin)
% v06:   -implementing 'BAfc_psth_spx'
%        -not using version in the function name 
%        -parser input added
% v05:   -removing gaussian filter
%        -using mean zscore in 'time_window' for detection
% v04: using parfor loop, faster


%% Default params
prs =  inputParser;
addRequired(prs,'experiment',@ischar)
addRequired(prs,'ttl',@ischar)
addRequired(prs,'pre_time', @isnumeric)
addRequired(prs,'baseline_time',@isnumeric)
addRequired(prs,'post_time', @isnumeric)
addRequired(prs,'bin_time', @isnumeric)
addRequired(prs,'test_win', @isnumeric)
addParameter(prs,'TTLinclude',0,@isnumeric) 
addParameter(prs,'TTLshift',0,@isnumeric) % shifting the PSTH compared to TTL, negative or positive number
parse(prs,varargin{:})
g = prs.Results;

%% Response calculation
% z-score, Smoothdata, Offset
[psth_spx] =  BAfc_psth_spx(g.experiment, g.ttl, g.pre_time, g.post_time, g.bin_time, 'TTLinclude', g.TTLinclude);

baseline_time = g.baseline_time/g.bin_time;
smoothed_zscore = zeros(size(psth_spx));
parfor ii = 1:size(psth_spx,1)   
    tempdata = smoothdata(psth_spx(ii,:),'movmean', 10);
    smoothed_zscore(ii,:) = (tempdata - mean(tempdata(1,1:baseline_time)))/std(tempdata(1,1:baseline_time));
end


timewin = g.pre_time/g.bin_time+1:(g.pre_time+g.test_win)/g.bin_time;
% Find lowFR neurons in 'test_win'
window_spx = psth_spx(:,timewin);
[idx_lowFR,~] = find(sum(window_spx,2)<10);
idx_lowFR = unique(idx_lowFR);

% Find responses in 'test_win'
window = smoothed_zscore(:,timewin);
[idx_exc, ~] = find(mean(window,2) >= 0.5);
idx_exc = unique(idx_exc);
idx_exc = setdiff(idx_exc, idx_lowFR); % remove lowFr only from excitatiory responses.

% Find inhibitory responses in 'test_win'
[idx_inh,~] = find(mean(window,2)<=-0.5);
idx_inh = unique(idx_inh);


idx = sort([idx_exc; idx_inh]);
