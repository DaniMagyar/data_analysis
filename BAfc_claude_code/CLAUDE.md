# CLAUDE.md

## Communication Style
**IMPORTANT**: Minimal text unless requested. Save tokens. Code changes and brief confirmations only.

## Project Overview
MATLAB scripts for basal amygdala (BA) fear conditioning neuronal data with optogenetic manipulation. K-means clustering on PSTH data to identify response patterns.

## Directory Structure
- **BAfc_claude_code/**: Clustering analysis scripts
- **neuron_analysis/BA_fear_cond_project/**: Core BAfc helper functions
- **Data**: `C:\Users\dmagyar\Desktop\BA_fear_cond\[recording_name]\kilosort25preprocess\`

## Core Functions (../neuron_analysis/BA_fear_cond_project/)

**BAfc_load_neurons()**: Loads cell_metrics from sessions
- Inputs: recordings, ttl, mainFolder
- Animal ID format: `MD###_###` from `MD###_###_kilosort`

**BAfc_psth_spx()**: Calculates PSTH
- Returns: `psth_spx` (neurons × bins), `preAP_norm/postAP_norm` (trial timestamps), `preAP_bin/postAP_bin` (binned per trial)
- Enables raster plots and single-trial analyses

**BAfc_putative_cellTypes()**: Excitatory/inhibitory classification

**BAfc_colors()**: Standard color scheme

### TTL Event Types
- `triptest_sound_only` / `triptest_sound_only_light`: CS ± light
- `triptest_shocks_only` / `triptest_shocks_only_light`: US ± light
- `triptest_both` / `triptest_both_light`: CS+US ± light

## Main Clustering Scripts

### K_means_clustering_light_claude_002_05vs05_zscore.m
Clusters stimulus-excited neurons with optogenetics.

**4 neuron groups**:
1. K-means clusters (1 to g.clustnum): Excited neurons
2. Light-inhibited (g.clustnum+1): Inhibited during pre-stim light
3. Stimulus-inhibited (g.clustnum+2): Inhibited by stimulus
4. Non-responsive (g.clustnum+3): No response

**Key params**: `g.clustnum`, `g.roi_pca` (−0.5 to +0.5s), `g.bR` (brain region), `g.bin_time` (0.01s), `g.smoothvalue` (5)

**Outputs**: Side-by-side heatmaps, z-scored/Hz traces, delta maps, scatter plots

### comparison_claude_003_supervised.m
Compares CS, US, CS+US without optogenetics. Includes elbow/silhouette for optimal k.

## Population Decoding Scripts

### BAfc_stimulus_decoder_v9.m
LDA/logistic decoder for binary stimulus discrimination. Cross-validated, temporal generalization matrices, permutation testing.

**Key params**: `g.cell_type_filter`, `g.brain_regions`, `ttl` (stimulus pair), `g.bin_time` (0.005s)

### BAfc_single_neuron_discrimination_v1.m
ROC/AUC and d-prime for single neurons. Sliding window (25-50ms), mixed-effects models.

**Key params**: `g.cellType`, `g.window_size` (0.05s), `g.step_size` (0.01s)

### BAfc_correlation_CS_US_CSandUS.m
Correlations between CS, US, CS+US responses (z-scored, 0-1s window).

### BAfc_trajectory.m
PCA trajectories for LA/BA, PN/IN. Calculates trajectory length, Mahalanobis distance, displacement angles.

### BAfc_return_to_baseline.m
Response offset timing analysis by cell type and region.

### BAfc_figure_optogenetics.m
Publication figure: PV silencing effects on short latency responses (0-25ms).

**Figure 1**: Example PN/IN rasters (US ± light)
**Figure 2**: Population short latency analysis (CS/US responsive neurons, 4×4 grid)
- Peak z-score, response probability, latency
- Color-coded by change (red >20%, blue <−20%, gray no change)
- Pie charts, stats on changed neurons only

**Key params**:
- `idx_example_IN/PN`: Example neurons
- `idx_CS/US_responsive`: Response groups
- `g.monosyn_window` (0.025s), `g.zscore_baseline_time` (5s), `g.change_threshold` (0.20)

## Waveform Classification (uni_classifier/)

**waveform_classifier_main.m**: Entry point (fresh or load session)
**calculate_waveform_features.m**: Trough-to-peak and half-width
**waveform_curator_ui.m**: Interactive GUI

- 10x supersampling for precision
- Draggable markers (trough/peak/half-width)
- Brain region filtering (LA/BA/CeA/Astria)
- Outputs: `trough_to_peak_corrected_ms`, `half_width_corrected_ms`

## Monosynaptic Response Analysis (uni_optogenetics/)

### BAfc_monosyn_optogenetics.m
Identifies monosynaptically activated neurons.

**Criteria**: Peak z-score ≥ threshold, response probability ≥ threshold, latency calculated
**Key params**: `g.monosyn_window` (0.025s), `g.zscore_threshold` (10), `g.prob_threshold` (0.20), `g.bR`, `g.cellType`
**Outputs**: `monosyn_results` struct, launches raster UI

### BAfc_monosyn_raster_ui.m
Interactive raster viewer (4 TTL types, −50 to +50ms). Shows responsive neurons, z-score, probability, latency.

### BAfc_match_celltypes()
Matches `cell_metrics` to `putative_celltypes.mat` by animal/brainRegion/cellID/cluID (not batchIDs).

## Modular Optogenetics System (2025-10-20 update)

**BAfc_optogenetics_params.m**: Centralized parameters
- **Two-rule responsiveness**: Rule 1 (z≥5, p≥0.25) OR Rule 2 (z≥10, p≥0.1)
- **Light-inhibited neuron detection**: Detects pre-stimulus inhibition during light trials
  - Method: 'ttest' or 'zscore'
  - Windows: recent (−0.5 to 0s) vs baseline (−5 to −0.5s)
  - Thresholds: p<0.05 (ttest) or z<−2 (zscore)
  - Auto-exclusion from responsive detection
- Z-score method: 'peak' or 'mean'
- Artifact handling for shocks (2-12ms exclusion)

**BAfc_identify_responsive_neurons.m**: Classification function
- Inputs: cell_metrics, ttl_list, params
- Returns: results struct with `all_*` fields, `responsive_idx`, `light_inhibited_neurons`, `light_inhibited_zscores`
- Two-rule system prevents missing high z-score/low probability neurons
- Light-inhibited neurons stored separately and excluded

**BAfc_figure_optogenetics.m**: Publication figures
- **Figure 3 (new)**: Response change visualization (2×3 grid)
  - Col 1: Paired scatter (no-light vs light, unity line)
  - Col 2: Bar chart (increased/decreased/no change proportions)
  - Col 3: Violin plot (percent change distribution)
  - Separates increased/decreased for proper stats

**BAfc_monosyn_raster_ui.m**: Interactive raster viewer
- Two-rule responsive detection (Rule 1 OR Rule 2)
- Red title and latency line for responsive neurons
- 4 TTL types, −50 to +50ms window
- Shows z-score, spike count, probability, latency

**Usage**:
```matlab
params = BAfc_optogenetics_params();
params.zscore_threshold_strict = 8;  % Adjust Rule 2
params.light_inhib_method = 'zscore';  % or 'ttest'
results = BAfc_identify_responsive_neurons(cell_metrics, ttl, params);
```

## Monosynaptic Response Visualization (uni_monosynaptic/)

### BAfc_figure_monosynaptic.m
Visualizes responsive neurons across stimulus types.

**Layout**: 3×3 grid (rows = stimuli, columns = metrics)
- **Row 1**: CS (Sound) responses
- **Row 2**: US (Shock) responses
- **Row 3**: CS+US (Both) responses

**Columns**: Response magnitude (z-score), probability, latency
**Groups**: LA-PN, LA-IN, BA-PN, BA-IN

**Plot style**: Individual neurons (scatter) + mean ± SD (errorbar)
- Uses SD (not SEM) per neuroscience convention when showing individual points
- Jittered x-positions for visibility
- Labels show group and sample size: "LA-PN (n=24)"

## Common Modifications

**Clusters**: `g.clustnum = 5;` (line 42)
**Region**: `g.bR = 'LA';` (line 45)
**Stimulus**: Edit `ttl` variable (lines 20-22)
**Time window**: `g.roi_pca = round((g.pre_time-0.5)/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);`

## Implementation Notes

- K-means: 20 replicates, squared Euclidean, MaxIter 500
- Z-score heatmaps: [−5.5, 13], delta: [−3, 3]
- Neuron ordering: centroid distance (excited), activity (inhibited), FR (non-responsive)
- Scripts are standalone (not functions), clear workspace on run
