function params = BAfc_monosyn_params()
% BAfc_monosyn_params - Centralized parameters for monosynaptic response analysis
%
% SYNTAX:
%   params = BAfc_monosyn_params()
%
% OUTPUTS:
%   params - Structure containing all analysis parameters
%
% DESCRIPTION:
%   Returns default parameters for monosynaptic response detection.
%   Modify this function to change default behavior across all scripts.
%
% USAGE:
%   % Use defaults
%   params = BAfc_monosyn_params();
%
%   % Modify specific parameters
%   params = BAfc_monosyn_params();
%   params.zscore_threshold = 8;
%   params.brain_regions = {'LA', 'Astria'};
%
% SEE ALSO: BAfc_identify_responsive_neurons_v2
%
% AUTHOR: Claude Code
% DATE: 2025-10-28

%% ===== NEURON SELECTION =====

% Cell type filter: 'PN', 'IN', or 'All'
params.cellType = 'All';

% Brain regions to analyze: cell array of region names or 'All'
% Valid regions: 'LA', 'BA', 'Astria', 'CeA'
params.brain_regions = 'All';  % Will include all available regions

%% ===== TEMPORAL PARAMETERS =====

% PSTH time windows
params.pre_time = 5;           % Pre-stimulus time (s)
params.post_time = 5;          % Post-stimulus time (s)
params.bin_time = 0.001;       % Bin size (s) - 1ms resolution

% Monosynaptic response window
params.monosyn_window = 0.025; % Window for monosynaptic responses (s) - 0-25ms

% Stimulus artifact handling (for shock stimuli)
params.artifact_pre = 0.002;   % Pre-stimulus artifact window (s) - 2ms
params.artifact_post = 0.012;  % Post-stimulus artifact window (s) - 12ms

%% ===== RESPONSIVENESS CRITERIA (TWO-RULE SYSTEM) =====

% Rule 1: Standard criterion
params.zscore_threshold = 5;          % Minimum z-score
params.prob_threshold = 0.25;         % Minimum response probability (25%)

% Rule 2: Strict z-score with lenient probability
params.zscore_threshold_strict = 10;  % High z-score threshold
params.prob_threshold_lenient = 0.10; % Lenient probability (10%)

% A neuron is classified as responsive if it meets EITHER Rule 1 OR Rule 2

%% ===== Z-SCORE CALCULATION METHOD =====

% Method for calculating z-score in monosynaptic window
% 'peak': Maximum z-score in window (default)
% 'mean': Mean z-score in artifact-free portion of window
params.zscore_method = 'peak';

%% ===== LIGHT-INHIBITED NEURON DETECTION =====

% Enable detection and exclusion of light-inhibited neurons
params.detect_light_inhibited = true;
params.exclude_light_inhibited = true;

% Detection method: 'ttest' or 'zscore'
params.light_inhib_method = 'zscore';

% Time windows for comparison
params.light_inhib_window_recent = 0.5;     % Recent pre-stimulus window (s)
params.light_inhib_window_baseline = 4.5;   % Baseline window (s)

% Thresholds
params.light_inhib_zscore_threshold = -2;   % Z-score threshold (negative = inhibition)
params.light_inhib_p_threshold = 0.05;      % P-value threshold (for ttest method)

%% ===== DATA PROCESSING =====

% Smoothing for z-scored PSTHs
params.smoothvalue = 5;  % Savitzky-Golay filter window (bins)
                         % Set to 0 to disable smoothing

%% ===== VISUALIZATION PARAMETERS =====

% Heatmap parameters
params.heatmap_clim = [-5.5 13];    % Color limits for z-score heatmaps
params.heatmap_plotwin = 0.1;       % Time window to display (±s around stimulus)

% Font sizes
params.fontSize_title = 13;
params.fontSize_axis = 12;
params.fontSize_label = 11;

% Line widths
params.linewidth_stimulus = 2;
params.linewidth_trace = 1.5;

%% ===== ANALYSIS OPTIONS =====

% Print detailed output
params.verbose = true;

% Save intermediate results
params.save_results = false;
params.results_path = './results/';

end
