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
Comprehensive monosynaptic response analysis with heatmaps and quantitative metrics.

**Figure 1: Heatmap Visualization** (4×3 grid)
- **Groups**: LA-PN, LA-IN, BA-PN, BA-IN (rows)
- **Stimuli**: CS (Sound), US (Shock), CS+US (Both) (columns)
- **Time window**: -0.1 to 0.1s around stimulus onset
- **Neuron inclusion**: Union of all neurons responsive to ANY stimulus within each group
- **Consistent ordering**: Same neuron list and order across all 3 stimuli within each group
- **Performance**: Pre-calculates all PSTHs once (3 total vs 12), ~4x faster
- **Color scale**: Z-score [-5.5, 13], standard BAfc heatmap colormap

**Sorting logic** (within each group, top to bottom):
1. **Category 1 (CS only)**: Neurons responsive to CS but not US, sorted by CS latency (earliest first)
2. **Category 2 (US only)**: Neurons responsive to US but not CS, sorted by US latency (earliest first)
3. **Category 3 (Both CS and US)**: Neurons responsive to both, sorted by CS+US latency (earliest first, CS latency fallback)

**Figure 2: Population Latency Analysis** (2×3 grid)
- Row 1: Latency histograms for CS, US, CS+US (0-50ms, 1ms bins)
- Row 2: CS vs US latency scatter, latency difference scatter plot

**Figure 3: Convergence, Integration, and Selectivity** (3×3 grid)
- **Convergence proportions**: Stacked bars (CS only/US only/Both %)
- **Integration type**: Stacked bars (Supralinear/Linear/Sublinear %)
- **Selectivity index**: Stacked bars (CS-selective/Non-selective/US-selective %)
- **Response magnitude**: Scatter+errorbar comparing CS vs US vs CS+US peak z-scores
- **Integration index scatter**: Individual neurons with thresholds (±0.2)
- **Selectivity index scatter**: Individual neurons with thresholds (±0.2)

**Quantitative Analyses** (printed to console):

1. **CS-US Convergence Metrics**
   - % neurons: CS only, US only, Both
   - Convergence index per group

2. **Integration Type Analysis** (dual-responsive neurons)
   - Supralinear: CS+US > (CS+US) by >20%
   - Linear: CS+US ≈ (CS+US) within ±20%
   - Sublinear: CS+US < (CS+US) by >20%
   - Integration index = (CS+US_actual - CS+US_predicted) / CS+US_predicted

3. **Temporal Precision & Reliability**
   - Mean latency ± SEM per stimulus/group
   - Population latency SD (variability across neurons)
   - Response probability (mean across neurons)

4. **Latency Comparisons** (dual-responsive neurons)
   - CS vs US latency difference (paired t-test)
   - % neurons with faster CS vs US response

5. **Response Magnitude Analysis** (dual-responsive neurons)
   - Peak z-scores: CS, US, CS+US (0-50ms window)
   - Statistical comparisons (paired t-tests)
   - Summation test: CS+US_actual vs CS+US_predicted

6. **Selectivity Index** (dual-responsive neurons)
   - Selectivity = (CS - US) / (CS + US)
   - Positive (>0.2) = CS-selective
   - Negative (<-0.2) = US-selective
   - Near-zero (±0.2) = Non-selective

7. **Cell Type/Region Comparisons**
   - LA vs BA convergence rates
   - PN vs IN convergence rates
   - Latency comparisons per group
   - Integration index comparisons

**Console output**: Category breakdown, statistics with significance stars (*, **, ***)

## Publication Figures (figures/figure_1/)

### BAfc_figure_1.m
Final publication figure with modular panel organization. A4-optimized layout (1500×1500px).

**Structure**: 4×4 tile layout, each panel in separate function

**Row 1**: Experimental setup (2 cols), Example traces (duplicated, 1 col each)
**Row 2**: PN ISI, IN ISI, Spike features, Normalized waveforms
**Row 3**: Firing rate distributions (LA, BA, Astria, CeA)
**Row 4**: Burst index distributions (LA, BA, Astria, CeA)

**Key features**:
- PN/IN plotted together with color coding
- ACG/waveform insets on ISI plots (top right corner)
- Log-scaled axes for FR/burst/ISI (3 ticklabels: min, mid, max)
- Firing rate: `g.cell_metrics.firingRate`
- Burst index: `g.cell_metrics.burstIndex_Royer2012`
- Waveforms: LineWidth 3, symbols MarkerSize 6

**Spike features panel** (Row 2, col 3):
- X-axis: Trough-to-peak (ms)
- Y-axis: Firing rate (Hz) from `g.cell_metrics.firingRate`
- Dashed lines: L-shape (horizontal from left to 0.4 at y=10, vertical from 10 to top at x=0.4)
- Separates fast-spiking high-FR neurons (upper left) from regular-spiking low-FR (lower right)

**Fonts**: Title 14pt, Axis 12pt

### BAfc_figure_2.m (figures/figure_2/)
Multi-region comparison figure with heatmaps, cluster PSTHs, and proportion analysis.

**Layout**: 4 rows (LA, BA, Astria, CeA) × 5 columns
- Columns 1-2: Nested CS+US heatmaps (z-score, side-by-side, spanning 2 cols)
- Column 3: Stacked bar chart (half width)
- Columns 4-5: Separate CS and US lineplots (firing rate Hz, spanning 2 cols)

**Brain region filters**:
- LA: all neurons
- BA: all neurons
- Astria: all neurons
- CeA: all neurons

**Clustering**: Same 5-cluster manual system as comparison_heatmap_005
- Cluster 1: CS-selective (red: 0.8 0.2 0.2)
- Cluster 2: US-selective (blue: 0.2 0.4 0.8)
- Cluster 3: Multisensory (purple: 0.6 0.2 0.6)
- Cluster 4: Non-responsive (gray: 0.6 0.6 0.6)
- Cluster 5: Inhibited (green: 0.2 0.6 0.3)
- Rank score sorting within clusters
- Inhibition: 50% FR drop threshold
- Excitation: z-score ≥2

**Heatmaps** (columns 1-2):
- CS and US heatmaps nested horizontally with compact spacing
- Titles: "CS" and "US" (bold)
- Shared global color limits (99th percentile)
- Black cluster boundary lines (1pt)
- Brain region labels (vertical, left side, FontSize 16, bold)
- Y-ticks: first and last neuron only
- Y-label: "Neuron #" (LA only)
- X-ticks: [-1 0 1] (bottom row only)
- Manual colorbar (left side, width 0.01, "Z-score" label)

**Stacked bar chart** (column 3):
- Vertical stacked bar (BarWidth: 0.5, positioned left at x=0.5, xlim [0 2])
- Order: bottom-to-top matches heatmap (cluster 5→1)
- Percentage labels on segments (white, bold)
- Colored text legend (right side, matching lineplot rows, FontSize 12, bold)
- Only shows clusters present in that brain region
- Total n at bottom center
- No axes visible

**Cluster lineplots** (columns 4-5):
- Separate CS (col 4) and US (col 5) panels with titles (bold)
- **Firing rate in Hz** (converted: psth_spx / (num_trials * bin_time))
- 5 stacked panels per column (1 per cluster)
- Mean firing rate trace per cluster (LineWidth 2.5, cluster color)
- **Separate y-limits**: clusters 1-4 (responsive) vs cluster 5 (inhibited, doubled: 2.2× multiplier)
- Dashed vertical line at t=0 (LineWidth 2)
- Y-axis hidden, no labels
- X-axis: [-1 0 1] on bottom panels only (CeA, cluster 5)
- Empty clusters: axis kept (for continuous dashed line), no trace
- **Scalebars**: 20 Hz on cluster 4 (non-responsive), 5 Hz on cluster 5 (inhibited), positioned near top

**Key features**:
- Figure size: 1500×3000px (tall)
- Global color limits (99th percentile) across all heatmaps
- Firing rate properly calculated (Hz): `psth_spx / (num_trials * bin_time)`
- Smoothing: Savitzky-Golay (201 points)
- Time window: ±2s display

**Fonts**: Title 14pt, Axis 12pt, Lineplot axis 12pt, Scalebar label 10pt

### BAfc_figure_3.m (figures/figure_3/)
CS+US convergence figure showing responsive neurons (CS-selective, US-selective, Multisensory) with quantitative metrics.

**TTL**: `triptest_sound_only`, `triptest_shocks_only`, `triptest_both`

**Layout**: 7×2 grid (1500×3000px)
- **Rows 1-2**: LA heatmaps (col 1, span 2 rows) | LA lineplots (col 2, span 2 rows)
- **Rows 3-4**: Astria heatmaps (col 1, span 2 rows) | Astria lineplots (col 2, span 2 rows)
- **Rows 5-7**: Metrics bar charts (span both columns, nested 3×3 grid)

**Brain regions**: LA, Astria only

**Clustering**: Identical to figure 2 (all 5 clusters computed)
- Classification based on CS and US responses only (not CS+US)
- Same rank score sorting as figure 2
- **Only clusters 1-3 plotted** (CS-selective, US-selective, Multisensory)
- Same neuron ordering as figure 2 for first 3 clusters

**Heatmaps** (column 1):
- Shows CS, US, and CS+US responses for same neurons (3 nested heatmaps horizontally)
- Same neurons in same order across all 3 stimuli
- Global color limits (99th percentile) across all 6 heatmaps (LA and Astria combined)
- Black cluster boundary lines (1pt)
- Brain region labels (vertical, left side, first heatmap only)
- Y-ticks: first and last neuron only
- Y-label: "Neuron #" (LA first heatmap only)
- X-label: "Time (s)" (Astria only)
- X-ticks: [-1 0 1] (Astria only)
- Manual colorbar (left side of Astria heatmap, cb_left=0.025, cb_bottom=0.48, width 0.01, height 0.15, "Z-score" label)

**Lineplots** (column 2):
- 3 stimuli (CS, US, CS+US) × 3 clusters (stacked)
- Titles: "CS", "US", "CS+US" (LA only, bold)
- **Firing rate in Hz** (converted: psth_spx / (num_trials * bin_time))
- Mean firing rate trace per cluster (LineWidth 2.5, cluster color)
- Shared y-limits across all 3 stimuli (CS, US, CS+US)
- 20 Hz scalebar on cluster 3 (Multisensory) of CS+US panel, positioned near top

**Metrics bar charts** (rows 5-7, nested 3×3 grid):
- **Rows**: CS-sel, US-sel, Multisensory (cluster names on left column only)
- **Columns**: Response AUC, Response Peak, Peak Time (titles on top row, bold)
- **Bar arrangement**: LA bars [2 3 4], gap, Astria bars [7 8 9]
  - Position 2/7: CS response
  - Position 3/8: US response
  - Position 4/9: CS+US response
  - Bar colors: Cluster-specific (CS-sel red, US-sel blue, Multi purple)
    - LA: Darker version (cluster_color × 0.7)
    - Astria: Lighter version (cluster_color + (1-cluster_color) × 0.5)
  - Bar width: 0.6 (narrower)
  - xlim: [0.5 10.5] (margins on both sides)
  - Gap between LA and Astria: 2 positions
- **X-tick labels**: Only on bottom row (Multisensory), "CS/US/CS+US" repeated twice
- **Y-axis**:
  - Left column (AUC): ylabel = "CS-sel/US-sel/Multi\nZ-score"
  - Middle column (Peak): ylabel = "Z-score"
  - Right column (Time): ylabel = "Time (ms)"
  - Y-ticks: Bottom, middle, top (rounded to 1 decimal) on all columns
- **Calculations**:
  - Response AUC: sum of positive z-scores in ROI
  - Response Peak: max z-score in ROI
  - Peak Time: latency to max z-score in ROI (ms)
- **Statistics**: Mean ± SEM with error bars (LineWidth 1.5, CapSize 4)
- **LA/Astria labels**: Text inside top row panels
  - Column 1 (AUC): 5% from top
  - Columns 2-3 (Peak, Time): 10% from top (lower position)
  - Top row ylim extended by 15% to accommodate labels

**Key features**:
- Figure size: 1500×3000px (adjustable height)
- **Figure visibility**: Use `'Visible', 'off'` to create figures larger than screen resolution
- **Position order**: Must set 'Units' before 'Position': `figure('Units', 'pixels', 'Position', [...])`
- Same clustering and sorting as figure 2
- Shows CS+US integration for responsive neurons
- Quantitative metrics for cross-region comparison
- Smoothing: Savitzky-Golay (201 points)
- Time window: ±2s display

**Fonts**: Title 14pt, Axis 12pt

## Clustering Analysis (BAfc_claude_code/)

### comparison_heatmap_004.m
Clean heatmap-only comparison of CS vs US responses. PNs and INs clustered separately, plotted together (PNs top, INs bottom).

**Two clustering methods**:

**1. K-means clustering** (data-driven)
- Clusters: 5 (adjustable via `g.clustnum`)
- ROI: ±0.5s around stimulus onset
- Sorting: onset latency → offset latency → distance to centroid

**2. Manual clustering** (hypothesis-driven)
- Cluster 1: CS-selective (CS excited, US not)
- Cluster 2: US-selective (US excited, CS not)
- Cluster 3: Multisensory (both CS and US excited)
- Cluster 4: Inhibited (not excited, but inhibited by ≥1 stimulus)
- Cluster 5: Non-responsive (neither excited nor inhibited)
- Sorting: onset latency → offset latency → peak response

**Key parameters**:
- `g.excitation_threshold = 2` (z-score)
- `g.inhibition_threshold = -2` (z-score)
- `g.use_percentile = true` (colormap limits)
- `g.clim_percentile = 99` (2.5th to 97.5th percentile)
- `g.smoothvalue = 101` (Savitzky-Golay filter)

**Z-scoring**: Baseline-only (pre-stimulus period) mean/SD used

**Performance**: PSTHs calculated once at start, reused throughout

**Latency detection**:
- Onset: First sustained crossing above threshold (≥20ms)
- Offset: First sustained return below threshold after onset
- Both used for sorting within clusters

**Output**: Two figures (k-means and manual), side-by-side CS/US heatmaps, thick black line separating PNs/INs, colored lines separating clusters

### comparison_heatmap_005.m
Rank score-sorted manual clustering with improved inhibition detection.

**Manual clustering**:
- Cluster 1: CS-selective (CS excited, US not)
- Cluster 2: US-selective (US excited, CS not)
- Cluster 3: Multisensory (both CS and US excited)
- Cluster 4: Non-responsive (neither excited nor inhibited)
- Cluster 5: Inhibited (firing rate drop ≥30%)

**Cluster ordering**: 1 → 2 → 3 → 4 → 5 (inhibited at bottom, just above INs)

**Sorting criteria** (rank score method):
- **Rank score**: `onset_latency + α * duration`, where `duration = offset_latency - onset_latency`
- **α = 0.8** (adjustable weight)
- Cluster 1: CS rank score (ascending)
- Cluster 2: US rank score (ascending)
- Cluster 3: Mean of CS and US rank scores (ascending)
- Cluster 4/5: Mean z-score during test time (descending)

**Inhibition detection** (firing rate-based):
- Compares baseline (−5 to 0s) vs test time (0 to +0.5s) in **Hz** (not z-score)
- `FR_drop = (baseline_FR - test_FR) / baseline_FR`
- Threshold: `g.inhibition_fr_drop = 0.30` (30% drop)
- Stored in `psthHz_full` (smoothed raw firing rates)

**INs classification**:
- Separated into responsive vs non-responsive
- Responsive: sorted by mean rank score (CS+US)
- Non-responsive: sorted by firing rate (descending)
- Gray line separates groups

**Key parameters**:
- `g.excitation_threshold = 2` (z-score for excitation)
- `g.inhibition_fr_drop = 0.30` (30% FR drop)
- `g.alpha = 0.8` (duration weight)
- `g.onset_threshold = 2` (z-score)
- `g.min_consec_bins = 20` (20ms sustained response)
- `g.plotwin = [1 1]` (±1s display)
- `g.smoothvalue = 101` (Savitzky-Golay)
- `g.clim_percentile = 95`

**Visual styling**:
- Cluster colors: Softer palette (0.8 0.2 0.2 for CS-sel, 0.2 0.4 0.8 for US-sel, etc.)
- Cluster borders: 2pt, PN/IN separator: 4pt black, IN separator: 2pt dark gray
- Colormap: `parula` or `g.colors.Heatmap`

**Output**: Side-by-side CS/US heatmaps, PNs top/INs bottom, colored cluster boundaries, rank score-sorted within clusters

## Common Modifications

**Clusters**: `g.clustnum = 5;` (line 42)
**Region**: `g.bR = 'LA';` (line 45)
**Stimulus**: Edit `ttl` variable (lines 20-22)
**Time window**: `g.roi_pca = round((g.pre_time-0.5)/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);`

## Implementation Notes

- K-means: 20 replicates, squared Euclidean, MaxIter 500
- Z-score heatmaps: baseline-referenced (pre-stimulus period only)
- Neuron ordering: onset latency (primary) → offset latency (secondary) → response magnitude/distance to centroid
- Scripts are standalone (not functions), clear workspace on run
- Baseline period: all time before stimulus (t < 0)
