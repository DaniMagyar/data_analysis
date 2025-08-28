function [psth_spx, preAP_bin_spx, postAP_bin_spx, preAP_norm_spx, postAP_norm_spx] =  BAfc_psth_spx(varargin)

%v02: parfor loop , faster
% Reworked with parser
% OUTPUT:   -psth_spx: spikes
%           -num_ttl: number of stimulations

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
psth_spx = zeros(num_cells, round((g.pre_time + g.post_time) / g.bin_time)); % adjust size accordingly
parfor ii = 1:num_cells
    %TTL = g.cell_metrics.general.(g.ttl){ii}(g.TTLinclude)+g.TTLshift; % commented out bc of incompatibility with Gergo's data. It might be unused anyways.
    TTL = g.cell_metrics.general.(g.ttl){ii}+g.TTLshift;
    AP = g.cell_metrics.spikes.times{ii};
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
            preAP_bin = [];
        else
            preAP_bin(:, tt) = histcounts(preAP_norm{tt}, linspace(-g.pre_time, 0, g.pre_time/g.bin_time + 1));
        end
        if g.post_time == 0
            postAP_bin = [];
        else
            postAP_bin(:, tt) = histcounts(postAP_norm{tt}, linspace(0, g.post_time, g.post_time/g.bin_time + 1));
        end
    end
    % Sum bins across TTLs
    psth_spx(ii,:) = sum([preAP_bin; postAP_bin], 2);
    preAP_bin_spx{ii} = sum(preAP_bin,1);
    postAP_bin_spx{ii} = sum(postAP_bin,1);
    preAP_norm_spx{ii} = cell2mat(preAP_norm);
    postAP_norm_spx{ii} = cell2mat(postAP_norm);
end