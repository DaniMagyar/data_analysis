% example_monosyn_analysis.m
% Example script showing how to use the cleaned monosynaptic analysis workflow
%
% This script demonstrates:
% 1. Loading parameters
% 2. Identifying responsive neurons
% 3. Analyzing different brain regions
%
% AUTHOR: Claude Code
% DATE: 2025-10-28

clear all; close all;

%% ===== SETUP =====

% Define recordings to analyze
recordings = {...
    'MD292_002_kilosort',...
    'MD293_001_kilosort',...
    'MD294_001_kilosort',...
    'MD295_001_kilosort',...
    'MD296_001_kilosort',...
    'MD297_001_kilosort',...
    'MD298_001_kilosort',...
    'MD299_001_kilosort',...
    'MD300_001_kilosort',...
    'MD304_001_kilosort',...
    'MD305_001_kilosort',...
    'MD307_001_kilosort',...
    'MD309_001_kilosort',...
    'MD310_001_kilosort',...
    'MD311_002_kilosort',...
    'MD312_001_kilosort',...
    'MD313_001_kilosort',...
    'MD314_001_kilosort',...
    'MD315_001_kilosort',...
    'MD316_002_kilosort',...
    'MD317_001_kilosort',...
    'MD318_001_kilosort',...
    'MD318_002_kilosort',...
    'MD319_003_kilosort'};

% Define stimulus types
ttl = {'triptest_sound_only', 'triptest_shocks_only', 'triptest_both'};

%% ===== LOAD DATA =====

fprintf('\n===========================================\n');
fprintf('LOADING NEURAL DATA\n');
fprintf('===========================================\n');

cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
cell_metrics = BAfc_match_celltypes('cell_metrics', cell_metrics);

fprintf('\nLoaded %d neurons from %d recordings\n', ...
    length(cell_metrics.cellID), length(recordings));

%% ===== EXAMPLE 1: Analyze all brain regions =====

fprintf('\n===========================================\n');
fprintf('EXAMPLE 1: ANALYZE ALL BRAIN REGIONS\n');
fprintf('===========================================\n');

% Get default parameters
params = BAfc_monosyn_params();

% Modify if needed
% params.zscore_threshold = 8;  % Example: use different threshold

% Identify responsive neurons
results_all = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl, params);

% Check which regions were analyzed
fprintf('\nAnalyzed regions: %s\n', strjoin(results_all.analyzed_regions, ', '));

% Display summary for each stimulus
for stim_idx = 1:length(ttl)
    ttl_fieldname = matlab.lang.makeValidName(ttl{stim_idx});
    n_responsive = length(results_all.(ttl_fieldname).responsive_idx);
    fprintf('  %s: %d responsive neurons\n', ttl{stim_idx}, n_responsive);
end

%% ===== EXAMPLE 2: Analyze specific brain regions (LA and Astria) WITH VISUALIZATION =====

fprintf('\n===========================================\n');
fprintf('EXAMPLE 2: ANALYZE LA AND ASTRIA WITH HEATMAPS\n');
fprintf('===========================================\n');

% Get fresh parameters
params = BAfc_monosyn_params();

% Specify brain regions
params.brain_regions = {'LA', 'Astria'};

% Identify responsive neurons
results_LA_Astria = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl, params);

fprintf('\nAnalyzed regions: %s\n', strjoin(results_LA_Astria.analyzed_regions, ', '));

% Display summary
for stim_idx = 1:length(ttl)
    ttl_fieldname = matlab.lang.makeValidName(ttl{stim_idx});
    n_responsive = length(results_LA_Astria.(ttl_fieldname).responsive_idx);
    fprintf('  %s: %d responsive neurons\n', ttl{stim_idx}, n_responsive);
end

% Plot heatmaps
fprintf('\nGenerating heatmap visualization...\n');
BAfc_plot_monosynaptic_heatmaps(cell_metrics, results_LA_Astria, ttl, params);
