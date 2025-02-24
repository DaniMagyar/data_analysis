function [latencies] = BAfc_calc_latency(varargin)

% Default params
prs =  inputParser;
addRequired(prs, 'experiment',@ischar)
addRequired(prs,'ttl',@ischar)
addRequired(prs,'pre_time', @isnumeric)
addRequired(prs,'post_time', @isnumeric)
addRequired(prs,'bin_time', @isnumeric)
addParameter(prs,'abs',0,@isnumeric)
addParameter(prs,'TTLinclude',0,@isnumeric) % 0 includes all
addParameter(prs, 'TTLshift',0,@isnumeric) % shifting the PSTH compared to TTL, negative or positive number
parse(prs,varargin{:})
g = prs.Results;

% Firstspike calculation
% z-score, Smoothdata
[psth_spx, ~, ~] =  BAfc_psth_spx(g.experiment, g.ttl, g.pre_time, g.post_time, g.bin_time, 'TTLinclude', g.TTLinclude);
% inverting for inhibitory detection
if g.abs == 1
    psth_spx = abs(psth_spx);
end
timewin = g.pre_time/g.bin_time+1:(g.pre_time+g.post_time)/g.bin_time;
for ii = 1:size(psth_spx,1)   
    smoothed_spx(ii,:) = smoothdata(zscore(psth_spx(ii,:)),'movmean', 10); 
    PeakBin = find(smoothed_spx(ii,timewin)...
        == max(smoothed_spx(ii,timewin)), 1, 'first');
    if isempty(PeakBin) 
        PeakSpike(ii,1) = 0;
        Estimate(ii,1) = 0;
    else
        PeakSpike(ii,1) = PeakBin;
        Estimate(ii,1) = find(smoothed_spx(ii,timewin)...
            >= mean([max(smoothed_spx(ii,timewin))...
            min(smoothed_spx(ii,timewin))]), 1, 'first');
    end
end

latencies = Estimate; % in ms