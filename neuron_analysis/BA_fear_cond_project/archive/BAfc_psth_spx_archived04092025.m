function [psth_spx, num_ttl, bR] =  BAfc_psth_spx(varargin)

%v02: parfor loop , faster
% Reworked with parser
% OUTPUT:   -psth_spx: spikes
%           -num_ttl: number of stimulations
%           -bR: brain region

%% Default params
prs =  inputParser;
addRequired(prs,'recordings',@iscell)
addRequired(prs,'ttl',@ischar)
addRequired(prs,'pre_time', @isnumeric)
addRequired(prs,'post_time', @isnumeric)
addRequired(prs,'bin_time', @isnumeric)
addParameter(prs,'TTLinclude',0,@isnumeric) 
addParameter(prs, 'TTLshift',0,@isnumeric) % shifting the PSTH compared to TTL, negative or positive number (negative shifts backward, positive forward)
parse(prs,varargin{:})
g = prs.Results;

cell_metrics = BAfc_load_neurons('recordings', g.recordings, 'ttl', g.ttl);
bR = cell_metrics.brainRegion;

num_ttlAll = size(cell_metrics.general.(g.ttl){1},1);
if ~any(g.TTLinclude)
    g.TTLinclude = 1:num_ttlAll;
end
num_ttl = numel(g.TTLinclude);
%% Preallocate arrays for results
num_cells = numel(cell_metrics.cellID);
psth_spx = zeros(num_cells, round((g.pre_time + g.post_time) / g.bin_time)); % adjust size accordingly
parfor ii = 1:num_cells
    TTL = cell_metrics.general.(g.ttl){ii}(g.TTLinclude)+g.TTLshift;
    AP = cell_metrics.spikes.times{ii};
    % Preallocate cell arrays for spikes
    preAP = cell(numel(TTL), 1);
    postAP = cell(numel(TTL), 1);
    preAP_norm = cell(numel(TTL), 1);
    postAP_norm = cell(numel(TTL), 1);
    % Binning variables
    preAP_bin = zeros(g.pre_time/g.bin_time, numel(TTL));
    postAP_bin = zeros(g.post_time/g.bin_time, numel(TTL));
    % Loop over TTL events
    for jj = 1:numel(TTL)
        % Spikes before and after each TTL
        preAP{jj} = AP(AP >= (TTL(jj) - g.pre_time) & AP < TTL(jj));
        postAP{jj} = AP(AP > TTL(jj) & AP < (TTL(jj) + g.post_time));
        % Normalize spike times to TTL
        preAP_norm{jj} = preAP{jj} - TTL(jj);
        postAP_norm{jj} = postAP{jj} - TTL(jj);
    end
    % Bin spike times
    for tt = 1:numel(TTL)
        if g.pre_time == 0
            preAP_bin = []
        else
            preAP_bin(:, tt) = histcounts(preAP_norm{tt}, linspace(-g.pre_time, 0, g.pre_time/g.bin_time + 1));
        end
        postAP_bin(:, tt) = histcounts(postAP_norm{tt}, linspace(0, g.post_time, g.post_time/g.bin_time + 1));
    end
    % Sum bins across TTLs
    psth_spx(ii,:) = sum([preAP_bin; postAP_bin], 2);
end