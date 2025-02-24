function [psth_spx, idx] =  BAfc_find_response(varargin)

% v01:  -shock, sound and allsound moved into the same function
% OUTPUT: not smoothed

%% Default params
prs =  inputParser;
addRequired(prs,'experiment',@ischar) % e.g.: 'NP_BAfc' or 'BAfc'
addRequired(prs,'stim',@ischar) % e.g.: 'shock', 'sound' or 'allsound'
addRequired(prs,'ttl',@ischar) % e.g.: 'TTL_shock', 'TTL_tone_habit_first'
addRequired(prs,'result',@ischar) % exc, inn or both
addRequired(prs,'pre_time',@isnumeric) % in sec
addRequired(prs,'base_time',@isnumeric) % usually == pre_time, but not when the baseline is not immediately before the TTL
addRequired(prs,'post_time',@isnumeric)
addRequired(prs,'test_win', @isnumeric) % e.g. look only the first 0.2 sec after stimulus
addRequired(prs,'test_value',@isnumeric)% zscore value for significance
addRequired(prs,'bin_time', @isnumeric)% 0.001 sec usually
addRequired(prs,'smoothvalue',@isnumeric)% 20-50 when bin_time is 0.001 sec
addParameter(prs,'brainRegion','all', @ischar)
addParameter(prs,'TTLinclude',0, @isnumeric) % e.g.: [1:5], [1:10]
addParameter(prs,'TTLshift',0, @isnumeric) % shifting the PSTH compared to TTL, negative or positive number
parse(prs,varargin{:})
g = prs.Results;

%% Z-score, Smoothdata, Offset
[psth_spx, num_ttl, bR] =  BAfc_psth_spx(g.experiment, g.ttl, g.pre_time, g.post_time, g.bin_time, 'TTLinclude', g.TTLinclude);
baseline_time = g.base_time/g.bin_time;
smoothed_zscore = zeros(size(psth_spx));
parfor ii = 1:size(psth_spx,1)   
    tempdata = smoothdata(psth_spx(ii,:),'movmean', g.smoothvalue);
    smoothed_zscore(ii,:) = (tempdata - mean(tempdata(1,1:baseline_time)))/std(tempdata(1,1:baseline_time));
end

%% Find responses
switch g.stim
    case {'shock', 'allsound'}       
        timewin1 = g.pre_time/g.bin_time+1:(g.pre_time+g.test_win)/g.bin_time; % to find lowFR exc responses
        %timewin2 = timewin1 - g.test_win/g.bin_time; % to find lowFR inh responses
        timewin2 = 1:g.pre_time/g.bin_time; % all baseline spikes
        window = smoothed_zscore(:,timewin1);
        % Find exc lowFR neurons in 'timewin1'
        window_spx1 = psth_spx(:,timewin1);
        [idx_lowFR1,~] = find(sum(window_spx1,2)/num_ttl<0.5); % using 0.5 response spike/ttl instead of a discrete number. 
        idx_lowFR1 = unique(idx_lowFR1); 
        % % Find inh lowFR neurons in 'timewin2'    
        % window_spx2 = psth_spx(:,timewin2);
        % [idx_lowFR2,~] = find(sum(window_spx2,2)<20); 
        % idx_lowFR2 = unique(idx_lowFR2);  
        %Find inh loWFR neurons in timewin2 
        window_spx2 = psth_spx(:, timewin2);
        [idx_lowFR2,~] = find(sum(window_spx2,2)/(g.pre_time*num_ttl)<1); % using 1 Hz baseline firing rate instead of a discrete number. 
        idx_lowFR2 = unique(idx_lowFR2); 

        % Find exc responses, remove lowFR from 'timewin1'   
        [idx_exc, ~] = find(mean(window,2)>=g.test_value);
        idx_exc = unique(idx_exc);
        idx_exc = setdiff(idx_exc, idx_lowFR1); 
        % Find inh responses, remove lowFR from 'timewin2'  
        [idx_inh,~] = find(mean(window,2)<=-g.test_value);
        idx_inh = unique(idx_inh);
        idx_inh = setdiff(idx_inh, idx_lowFR2);
    case 'sound'
end
% Pool responses
switch g.result
    case 'exc'
        idx = idx_exc;
    case 'inh'
        idx = idx_inh;
    case 'both'
         idx = sort([idx_exc; idx_inh]);
end
% Selects brainRegion
if ~strcmp(g.brainRegion, 'all')
    idx_region =  find(contains(bR, g.brainRegion));
    idx = intersect(idx, idx_region);
end
