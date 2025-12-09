# Monosynaptic Optogenetic Analysis - Clean Supervised Version

**Author**: Claude Code
**Date**: 2025-10-28

## Overview

This folder contains the cleaned, supervised version of the monosynaptic response analysis workflow. Key improvements over the original version:

1. **Fixed Astria bug**: Original code hardcoded 'LA' and 'BA' as only valid brain regions
2. **Centralized parameters**: All parameters in `BAfc_monosyn_params.m`
3. **Cleaner API**: Simplified function interfaces
4. **Better documentation**: Comprehensive examples and documentation

## Files

### Core Functions

- **`BAfc_monosyn_params.m`**: Centralized parameter configuration
  - Returns default parameters for monosynaptic analysis
  - Modify this file to change defaults globally
  - All parameters documented with comments

- **`BAfc_identify_responsive_neurons_v2.m`**: Fixed responsive neuron detection
  - **FIXED**: Now properly handles all brain regions (LA, BA, Astria, CeA, etc.)
  - Two-rule responsiveness system (Rule 1 OR Rule 2)
  - Light-inhibited neuron detection and exclusion
  - Supports flexible brain region filtering

### Visualization Functions

- **`BAfc_plot_monosynaptic_heatmaps.m`**: Heatmap visualization
  - Plots monosynaptic responses for all brain regions/cell types
  - Displays CS, US, and CS+US responses side-by-side
  - Neurons sorted by category (CS only, US only, Both) and latency
  - Automatic layout adjustment based on number of groups

### Example Scripts

- **`example_monosyn_analysis.m`**: Comprehensive usage examples
  - Example 1: Analyze all brain regions
  - Example 2: Analyze specific regions (LA + Astria) WITH HEATMAPS
  - Example 3: Analyze cell-type specific (e.g., PNs in LA)
  - Example 4: Access detailed neuron data
  - Example 5: Parameter comparison

## Bug Fix Details

### Original Bug (in `uni_optogenetics/BAfc_identify_responsive_neurons.m`)

```matlab
% Lines 92-97 - HARDCODED LA and BA only
if strcmp(params.bR, 'All')
    idx_region = strcmp(cell_metrics.brainRegion, 'LA') | ...
                 strcmp(cell_metrics.brainRegion, 'BA');  % BUG: Excludes Astria!
else
    idx_region = strcmp(cell_metrics.brainRegion, params.bR);
end
```

**Problem**: When `params.bR = 'All'`, it only includes LA and BA, completely ignoring Astria, CeA, or any other brain regions in the dataset.

### Fixed Version (in `BAfc_identify_responsive_neurons_v2.m`)

```matlab
% Lines 91-113 - FIXED: Properly handles all regions
if ischar(params.brain_regions) && strcmp(params.brain_regions, 'All')
    % Get all unique brain regions present in data
    unique_regions = unique(cell_metrics.brainRegion);
    fprintf('Brain regions found in data: %s\n', strjoin(unique_regions, ', '));

    % Include all neurons with valid brain region labels
    idx_region = false(size(cell_metrics.brainRegion));
    for i = 1:length(unique_regions)
        idx_region = idx_region | strcmp(cell_metrics.brainRegion, unique_regions{i});
    end

    analyzed_regions = unique_regions;

elseif iscell(params.brain_regions)
    % Specific list of regions provided
    idx_region = false(size(cell_metrics.brainRegion));
    for i = 1:length(params.brain_regions)
        idx_region = idx_region | strcmp(cell_metrics.brainRegion, params.brain_regions{i});
    end

    analyzed_regions = params.brain_regions;
```

**Fix**:
1. Automatically detects all unique brain regions in the dataset
2. Supports flexible specification: `'All'`, `'LA'`, or `{'LA', 'Astria'}`
3. Returns `analyzed_regions` in results for transparency

## Usage

### Basic Usage

```matlab
% Load data
cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);
cell_metrics = BAfc_match_celltypes('cell_metrics', cell_metrics);

% Get default parameters
params = BAfc_monosyn_params();

% Identify responsive neurons
results = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl, params);
```

### Analyze LA and Astria with Visualization

```matlab
% Get parameters
params = BAfc_monosyn_params();

% Specify regions
params.brain_regions = {'LA', 'Astria'};

% Run analysis
results = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl, params);

% Check analyzed regions
fprintf('Analyzed: %s\n', strjoin(results.analyzed_regions, ', '));

% Plot heatmaps
BAfc_plot_monosynaptic_heatmaps(cell_metrics, results, ttl, params);
```

### Analyze Only PNs in BA

```matlab
% Get parameters
params = BAfc_monosyn_params();

% Specify cell type and region
params.cellType = 'PN';
params.brain_regions = 'BA';

% Run analysis
results = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl, params);
```

### Modify Detection Thresholds

```matlab
% Get parameters
params = BAfc_monosyn_params();

% Modify responsiveness criteria
params.zscore_threshold = 8;              % Rule 1: more stringent
params.zscore_threshold_strict = 12;      % Rule 2: more stringent
params.prob_threshold = 0.30;             % Rule 1: higher probability

% Run analysis
results = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl, params);
```

## Parameter Reference

### Neuron Selection
- `cellType`: `'PN'`, `'IN'`, or `'All'`
- `brain_regions`: `'All'`, `'LA'`, `{'LA', 'Astria'}`, etc.

### Temporal Parameters
- `pre_time`: Pre-stimulus window (default: 5s)
- `post_time`: Post-stimulus window (default: 5s)
- `bin_time`: Bin size (default: 0.001s = 1ms)
- `monosyn_window`: Monosynaptic window (default: 0.025s = 25ms)

### Responsiveness Criteria (Two-Rule System)
**Rule 1** (Standard):
- `zscore_threshold`: 5
- `prob_threshold`: 0.25

**Rule 2** (Strict z-score, lenient probability):
- `zscore_threshold_strict`: 10
- `prob_threshold_lenient`: 0.10

A neuron is classified as responsive if it meets **EITHER** rule.

### Z-Score Calculation
- `zscore_method`: `'peak'` (max in window) or `'mean'` (mean in artifact-free portion)

### Light-Inhibited Neuron Detection
- `detect_light_inhibited`: Enable/disable detection
- `exclude_light_inhibited`: Exclude from responsive detection
- `light_inhib_method`: `'ttest'` or `'zscore'`
- `light_inhib_zscore_threshold`: -2 (negative = inhibition)

## Results Structure

```matlab
results
  .params                  % Copy of parameters used
  .ttl_types               % List of TTL types analyzed
  .num_analyzed            % Total neurons analyzed
  .analyzed_regions        % Brain regions included
  .light_inhibited_neurons % Indices of light-inhibited neurons
  .[ttl_name]              % For each stimulus type:
    .responsive_idx        % Indices of responsive neurons
    .all_peak_zscore       % Z-scores for all neurons
    .all_response_probability  % Probabilities for all neurons
    .all_mean_latency      % Latencies for all neurons (ms)
    .neuron_idx            % Indices (responsive only)
    .peak_zscore           % Z-scores (responsive only)
    .response_probability  % Probabilities (responsive only)
    .mean_latency          % Latencies (responsive only, ms)
```

## Compatibility

### Works With
- `BAfc_figure_monosynaptic.m` (change function call)
- `BAfc_figure_optogenetics.m` (change function call)
- Any script using monosynaptic response detection

### Migration from Original

Replace:
```matlab
results = BAfc_identify_responsive_neurons(cell_metrics, ttl, params);
```

With:
```matlab
params = BAfc_monosyn_params();  % Get clean params
% Modify params as needed
results = BAfc_identify_responsive_neurons_v2(cell_metrics, ttl, params);
```

**Note**: The original `BAfc_optogenetics_params()` and the new `BAfc_monosyn_params()` have similar structures but the new version is cleaner and uses `brain_regions` instead of `bR`.

## See Also

- Original version: `uni_optogenetics/BAfc_identify_responsive_neurons.m`
- Original params: `uni_optogenetics/BAfc_optogenetics_params.m`
- Figure scripts: `uni_monosynaptic/BAfc_figure_monosynaptic.m`
- Optogenetics figures: `uni_optogenetics/BAfc_figure_optogenetics.m`
