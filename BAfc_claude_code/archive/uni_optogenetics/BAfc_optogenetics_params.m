function params = BAfc_optogenetics_params()
% BAfc_optogenetics_params - Default parameters for optogenetics analysis
%
% SYNTAX:
%   params = BAfc_optogenetics_params()
%
% OUTPUTS:
%   params - Structure with default analysis parameters
%
% DESCRIPTION:
%   Returns default parameter structure used by BAfc_identify_responsive_neurons,
%   BAfc_monosyn_optogenetics, and BAfc_figure_optogenetics.
%
%   Modify parameters in this function to change defaults for all scripts.
%
% USAGE:
%   % Use defaults
%   params = BAfc_optogenetics_params();
%   results = BAfc_identify_responsive_neurons(cell_metrics, ttl, params);
%
%   % Modify specific parameters
%   params = BAfc_optogenetics_params();
%   params.zscore_threshold = 8;
%   params.bR = 'LA';
%   results = BAfc_identify_responsive_neurons(cell_metrics, ttl, params);
%
% AUTHOR: Claude Code
% DATE: 2025-10-19

%% ===== RESPONSIVENESS DETECTION PARAMETERS =====
% Used by BAfc_identify_responsive_neurons()

params.pre_time = 5;               % Baseline period (seconds)
params.post_time = 0.05;           % Post-stimulus period (seconds)
params.bin_time = 0.001;           % Bin size (1ms)
params.monosyn_window = 0.025;     % Monosynaptic response window (25ms)
params.zscore_method = 'peak';     % Z-score calculation method: 'peak' or 'mean'
% Two-rule system for responsiveness detection:
% Rule 1: zscore >= zscore_threshold AND prob >= prob_threshold
% Rule 2: zscore >= zscore_threshold_strict AND prob >= prob_threshold_lenient
params.zscore_threshold = 5;       % Z-score threshold for Rule 1
params.prob_threshold = 0.25;      % Probability threshold for Rule 1
params.zscore_threshold_strict = 10;   % Strict z-score threshold for Rule 2
params.prob_threshold_lenient = 0.1;   % Lenient probability threshold for Rule 2
params.bR = 'All';                 % Brain region filter ('LA', 'BA', 'All')
params.cellType = 'All';           % Cell type filter ('PN', 'IN', 'All')
params.smoothvalue = 5;            % Savitzky-Golay smoothing window (0 = no smoothing)

%% ===== FIGURE GENERATION PARAMETERS =====
% Additional parameters used by BAfc_figure_optogenetics()

params.pre_time_short = 0.1;           % 100ms baseline for raster plots
params.post_time_short = 0.1;          % 100ms post-stimulus for raster plots
params.short_latency_window = 0.025;   % 25ms short latency window
params.zscore_baseline_time = 5;       % 5s baseline for z-score calculation
params.pre_time_baseline = 5;          % 5s window for baseline calculation
params.change_threshold = 0.2;        % 20% change threshold for increased/decreased classification
params.pre_time_long = 2;              % 2s window for light inhibition verification
params.post_time_long = 2;             % 2s window

%% ===== ARTIFACT HANDLING PARAMETERS =====
% Parameters for shock stimulus artifact removal

params.artifact_pre = 0.002;       % 2ms before stimulus
params.artifact_post = 0.012;      % 12ms after stimulus

%% ===== LIGHT INHIBITION DETECTION PARAMETERS =====
% Parameters for detecting neurons directly inhibited by light in pre-stimulus period

params.light_inhib_window_recent = 0.5;     % Recent pre-stimulus window (s) - e.g., [-0.5, 0s]
params.light_inhib_window_baseline = 4.5;   % Baseline window duration (s) - e.g., [-5, -0.5s]
params.light_inhib_p_threshold = 0.05;      % P-value threshold for significance
params.light_inhib_method = 'zscore';        % 'ttest' or 'zscore'
params.light_inhib_zscore_threshold = -0.5;   % Z-score threshold (if using zscore method)
params.exclude_light_inhibited = true;      % Whether to exclude light-inhibited neurons

end
