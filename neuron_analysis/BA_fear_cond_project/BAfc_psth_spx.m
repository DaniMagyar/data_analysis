function [psth_spx, preAP_norm, postAP_norm, preAP_bin, postAP_bin] =  BAfc_psth_spx(varargin)

%v02: parfor loop , faster
% Reworked with parser
% OUTPUT:   -psth_spx: all spikes in all trialss. Each neuron is a row.
%           -preAP_norm, postAP_norm: normalized spike timestamps
%           -

%% Default params
prs =  inputParser;
addParameter(prs,'cell_metrics',[],@isstruct) % MUST CONTAIN: ttl in cell_metrics.general
addParameter(prs,'ttl',@ischar)
addParameter(prs,'pre_time', @isnumeric)
addParameter(prs,'post_time', @isnumeric)
addParameter(prs,'bin_time',[], @isnumeric)
addParameter(prs,'TTLinclude',0,@isnumeric) 
addParameter(prs, 'TTLshift',0,@isnumeric) % shifting the PSTH compared to TTL, negative or positive number (negative shifts backward, positive forward)
parse(prs,varargin{:})
g = prs.Results;

% num_ttlAll = size(g.cell_metrics.general.(g.ttl){1},1); % commented out bc of incompatibility with Gergo's data. It might be unused anyways.
% if ~any(g.TTLinclude)
%     g.TTLinclude = 1:num_ttlAll;
% end
% num_ttl = numel(g.TTLinclude);


%% Preallocate arrays for results
num_cells = numel(g.cell_metrics.cellID);
preAP_norm = cell(num_cells, 1);
postAP_norm = cell(num_cells, 1);
preAP_bin = cell(num_cells, 1);
postAP_bin = cell(num_cells, 1);
psth_spx = zeros(num_cells, round((g.pre_time + g.post_time) / g.bin_time)); % adjust size accordingly
for ii = 1:num_cells
    %TTL = g.cell_metrics.general.(g.ttl){ii}(g.TTLinclude)+g.TTLshift; % commented out bc of incompatibility with Gergo's data. It might be unused anyways.
    TTL = g.cell_metrics.general.(g.ttl){ii}+g.TTLshift;
    AP = g.cell_metrics.spikes.times{ii};
    % Preallocate cell arrays for spikes
    preAP = cell(numel(TTL), 1);
    postAP = cell(numel(TTL), 1);
    preAP_norm{ii} = cell(numel(TTL), 1);
    postAP_norm{ii} = cell(numel(TTL), 1);
    % Binning variables
    preAP_bin{ii} = zeros(g.pre_time/g.bin_time, numel(TTL));
    postAP_bin{ii} = zeros(g.post_time/g.bin_time, numel(TTL));
    % Loop over TTL events
    for jj = 1:numel(TTL)
        % Spikes before and after each TTL
        preAP{jj} = AP(AP >= (TTL(jj) - g.pre_time) & AP < TTL(jj));
        postAP{jj} = AP(AP > TTL(jj) & AP < (TTL(jj) + g.post_time));
        % Normalize spike times to TTL
        preAP_norm{ii}{jj} = preAP{jj} - TTL(jj);
        postAP_norm{ii}{jj} = postAP{jj} - TTL(jj);
        % Bin spike times
        if g.pre_time == 0
            preAP_bin{ii} = [];
        else
            preAP_bin{ii}(:, jj) = histcounts(preAP_norm{ii}{jj}, linspace(-g.pre_time, 0, g.pre_time/g.bin_time + 1));
        end
        if g.post_time == 0
            postAP_bin{ii} = [];
        else
            postAP_bin{ii}(:, jj) = histcounts(postAP_norm{ii}{jj}, linspace(0, g.post_time, g.post_time/g.bin_time + 1));
        end
    end
    % Sum bins across TTLs
    psth_spx(ii,:) = sum([preAP_bin{ii}; postAP_bin{ii}], 2);
end