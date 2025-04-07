function [psth_spx, idx, zResp_med] =  BAfc_find_response(varargin)
% OUTPUT:   -psth_spx: binned spikes of all neurons
%           -idx: index of responsive neurons
%           -zResp_med: median zscore responses in 'test_win'
%
% Author:   Daniel Magyar, Indiana University Bloomington

%% Default params
prs =  inputParser;
addRequired(prs,'experiment',@ischar) % e.g.: 'NP_BAfc' or 'BAfc'
addRequired(prs,'stim',@ischar) % e.g.: 'shock', 'sound' or 'allsound'
addRequired(prs,'ttl',@ischar) % e.g.: 'TTL_shock', 'TTL_tone_habit_first'
addRequired(prs,'result',@ischar) % exc, inn or both
addRequired(prs,'pre_time',@isnumeric) % in sec
addRequired(prs,'post_time',@isnumeric)
addRequired(prs,'test_win', @isnumeric) % e.g. [0.012 0.05] calculate mean from 0.012 to 0.05
addRequired(prs,'bin_time', @isnumeric)% 0.001 sec usually
addParameter(prs,'base_time',[], @isnumeric) % usually == pre_time, but not when the baseline is not immediately before the TTL
addParameter(prs,'psth_out',[], @isnumeric) % part of psth_spx to be outputted [post1: post2], e.g. [0.012 0.05] for 12-15ms, or [-0.012 0.05] for 0.012 pre and 0.05 post
addParameter(prs,'smoothvalue',20, @isnumeric)% 20-50 when bin_time is 0.001 sec
addParameter(prs,'test_value',1.96, @isnumeric)% zscore value for significance
addParameter(prs,'artefactLength',0, @isnumeric) % length of the artefacet after ttl in ms
addParameter(prs,'brainRegion','all', @ischar)
addParameter(prs,'TTLinclude',0, @isnumeric) % e.g.: [1:5], [1:10]
addParameter(prs,'TTLshift',0, @isnumeric) % shifting the PSTH compared to TTL, negative or positive number
parse(prs,varargin{:})
g = prs.Results;

%% Z-score, Smoothdata, Offset
[psth_spx, num_ttl, bR] =  BAfc_psth_spx(g.experiment, g.ttl, g.pre_time, g.post_time, g.bin_time, 'TTLinclude', g.TTLinclude);
if ~any(g.base_time)
    g.base_time = g.pre_time;
end
baseline_time = g.base_time/g.bin_time;
smoothed_zscore = zeros(size(psth_spx));
parfor ii = 1:size(psth_spx,1)   
    tempdata = psth_spx(ii,:);
    if any(g.artefactLength)
        tempdata(g.pre_time/g.bin_time+1:round((g.pre_time+g.artefactLength)/g.bin_time)) = NaN;
    end
    tempdata = smoothdata(tempdata,'movmean', g.smoothvalue);
    smoothed_zscore(ii,:) = (tempdata - mean(tempdata(1,1:baseline_time)))/std(tempdata(1,1:baseline_time));
end

%% Find responses
switch g.stim
    case {'shock', 'allsound'}     
        timewin1 = round((g.pre_time+g.test_win(1))/g.bin_time+1:(g.pre_time+g.test_win(2))/g.bin_time); % all response spikes
        timewin2 = 1:g.pre_time/g.bin_time; % all baseline spikes
        window = smoothed_zscore(:,timewin1);
        % Find exc lowFR neurons in 'timewin1'
        window_spx1 = psth_spx(:,timewin1);
        [idx_lowFR1,~] = find(sum(window_spx1,2)/num_ttl<0.4); % using 0.5 response spike/ttl instead of a discrete number. 
        idx_lowFR1 = unique(idx_lowFR1); 
        % % Find inh lowFR neurons in 'timewin2'    
        window_spx2 = psth_spx(:, timewin2);
        baseline_fr = sum(window_spx2,2)/(g.pre_time*num_ttl);
        [idx_lowFR2,~] = find(baseline_fr<1); % using 1 Hz baseline firing rate instead of a discrete number. 
        idx_lowFR2 = unique(idx_lowFR2); 
        % Find exc responses, remove lowFR from 'timewin1' 
        zResp_med = median(window,2);
        [idx_exc, ~] = find(zResp_med>=g.test_value);
        idx_exc = unique(idx_exc);
        idx_exc = setdiff(idx_exc, idx_lowFR1); 
        % Find inh responses, remove lowFR from 'timewin2'  
        [idx_inh,~] = find(zResp_med<=-g.test_value);
        idx_inh = unique(idx_inh);
        response_fr = sum(window_spx1,2)/((g.test_win(2)-g.test_win(1))*num_ttl);
        [idx_inh50,~] = find(response_fr./baseline_fr<0.5);
        idx_inh = setdiff(unique([idx_inh; idx_inh50]), idx_lowFR2);
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
% cut psth_spx
if any(g.psth_out)
    psth_spx = psth_spx(:, (g.pre_time+g.psth_out(1))/g.bin_time:(g.pre_time+g.psth_out(2))/g.bin_time);
end

