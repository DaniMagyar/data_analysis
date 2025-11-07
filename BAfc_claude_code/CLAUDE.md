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
Monosynaptic analysis with heatmaps and quantitative metrics.

**Fig 1 (4×3)**: LA-PN/IN, BA-PN/IN × CS/US/CS+US heatmaps. Sorting: CS-only (CS latency), US-only (US latency), Both (CS+US latency). Z-score [-5.5, 13], -0.1 to 0.1s window.

**Fig 2 (2×3)**: Latency histograms (0-50ms) and CS vs US scatter plots.

**Fig 3 (3×3)**: Convergence proportions, integration type (supralinear/linear/sublinear), selectivity index, response magnitude comparisons.

**Console Stats**: Convergence metrics, integration analysis (±20% thresholds), temporal precision, latency comparisons (paired t-tests), selectivity index ((CS-US)/(CS+US), ±0.2 thresholds), cell type/region comparisons.

## Publication Figures (figures/figure_1/)

### BAfc_figure_1.m (updated 2025-11-06)
Final publication figure with modular panel organization. A4-optimized layout (1000×1000px).

**Structure**: 4×4 tile layout, each panel in separate function

**Row 1**: Experimental setup (2 cols), Example CS response (tile 3), Example US response (tile 4)
**Row 2**: PN ISI, IN ISI, Spike features, Normalized waveforms
**Row 3**: Firing rate distributions (LA, BA, Astria, CeA)
**Row 4**: Burst index distributions (LA, BA, Astria, CeA)

**Example Traces** (Row 1, tiles 3-4):
- **Example CS response**: MD312_001, channel 122, `triptest_sound_only([2 3 4 5 6])`, time window [-40, +60]ms
  - Y-limits: [-0.12, 0.1], yticks: [-0.1, 0, 0.1]
  - Red speaker symbol (🔊) at x=7ms, top of plot
- **Example US response**: MD292_002, channel 32, `shocks([5 7 36 41 43 47])`, time window [-40, +60]ms
  - Y-limits: [-0.5, 0.3], yticks: [-0.5, -0.1, 0.3]
  - Gray shaded area (0-10ms, artifact period)
  - Red lightning bolt (⚡) at x=5ms, top of plot, FontSize 20
- **Both examples**: Dashed black line at x=0, X-axis label "Time (ms)", title "Example CS/US response"

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

**Fonts**: All text 10pt (titles, axes, labels, legends, annotations)

### BAfc_figure_2.m (figures/figure_2/, updated 2025-11-06)
Multi-region comparison: 4 rows (LA/BA/Astria/CeA) × 5 cols (1000×800px).

**Layout**: Cols 1-2 nested CS+US heatmaps (z-score, 99th percentile limits), col 3 stacked bar (cluster proportions), cols 4-5 separate CS/US lineplots (Hz).

**Clustering**: 5 clusters (CS-sel red, US-sel blue, Multi purple, Non-resp gray, Inhibited green). Rank score sorting. Excitation z≥2. Smoothing: Savitzky-Golay (201).

**Inhibition Detection (Two-Rule System)**:
- **Rule 1**: ≥50% FR drop in 0-0.2s window (short, transient inhibition)
- **Rule 2**: ≥50% FR drop in 0-1.0s window (long, sustained inhibition)
- Neuron classified as inhibited if **EITHER** rule satisfied (OR logic)
- Parameters: `g.inhibition_fr_drop` (0.50), `g.inhibition_window_short` (0.2s), `g.inhibition_window_long` (1.0s)
- Console output shows breakdown: Rule1-only, Rule2-only, Both rules

**Within-Region Cluster Merging**:
- Adaptive minimum cluster size: 2% of neurons in each brain region (`g.min_cluster_percent = 2`)
- Merges small clusters (<2% threshold) into nearest cluster by centroid distance (Euclidean on concatenated CS+US PSTHs)
- Region-specific thresholds prevent overly strict/lenient requirements
- Fallback: If all clusters small, merge into cluster 4 (Non-responsive)
- Console output shows cluster sizes (absolute + percentage) and merge operations

**Statistical Analysis**:
- **Chi-square test**: Contingency table (4 regions × 5 clusters) with permutation testing (10,000 iterations)
- **Cramér's V**: Global and pairwise effect sizes
- **FDR correction**: Benjamini-Hochberg for pairwise comparisons
- **All 4 regions included**: LA, BA, Astria, CeA

**Separate Figures**:

1. **US Metrics Figure** (`fig_US`, 1000×250px, 1×4 layout):
   - **Column 1**: Chi-square p-values matrix (4×4, all regions)
     - Binary color: White (p≤0.001), Light grey (p>0.001)
     - No colorbar, p-values and significance stars overlaid
   - **Columns 2-4**: US response metrics bar charts (3 regions only: LA, BA, Astria)
     - ΔnSpikes, ΔPeak FR (Hz), Response Length (ms)
     - US-selective and Multisensory neurons only (clusters 2 & 3)
     - Simple grey bars [0.6 0.6 0.6], mean ± SEM
     - Wilcoxon rank-sum tests (3 pairwise comparisons: LA-BA, LA-Astria, BA-Astria)
     - No y-labels (titles indicate metrics)

2. **Supplementary Statistical Figure** (`fig_stats`, 1000×300px, 1×4 layout):
   - **Panel 1**: Permutation test distribution
   - **Panel 2**: Contingency table (4 regions × 5 clusters, shows χ² statistic)
   - **Panel 3**: Cramér's V effect size matrix (4×4, `parula` colormap, clim [0, 0.5])
   - **Panel 4**: FDR correction scatter plot (uncorrected vs corrected p-values)

**Heatmaps**: Black cluster lines (1pt), manual colorbar (left, 0.01 width). **Bars**: Vertical stacked (BarWidth 0.5), percentages on segments. **Lineplots**: Hz (psth_spx/trials/bin_time), separate y-limits (clusters 1-4 vs 5: 2.2× multiplier), scalebars 20/5 Hz.

### BAfc_figure_3.m (figures/figure_3/, updated 2025-11-07)
CS+US convergence for LA/Astria responsive neurons (clusters 1-3 only). 6×2 grid (1000×1000px).

**Layout**:
- Rows 1-3: Heatmaps and lineplots (2×6 tiledlayout with tight spacing, 50% of figure)
  - Row 1: LA (6 panels)
  - Row 2: Astria (6 panels)
  - Columns 1-3: CS, US, CS+US heatmaps
  - Columns 4-6: CS, US, CS+US lineplots (3 clusters stacked in nested 3×1 layouts)
- Rows 4-6: Metrics bar charts (3×3 nested, 50% of figure)

**Heatmaps**: Same neurons across CS/US/CS+US, global 99th percentile limits, onset/offset dots (MarkerSize 4), manual colorbar (Astria left, [0.025, 0.52], 0.01×0.15).

**Lineplots**: 3 stimuli × 3 clusters stacked (nested tiledlayouts), Hz (psth_spx/trials/bin_time), shared y-limits, 20 Hz scalebar on cluster 3 of CS+US panel.

**Bar Charts** (3 metrics): ΔnSpikes, ΔPeak FR, Response length (onset-to-offset window, baseline-subtracted).
- **Positions**: LA [2 3 4], Astria [7 8 9]
- **Colors**: LA darker (0.7×), Astria lighter (+0.5 blend)
- **Error bars**: SEM (appropriate for group mean comparisons with statistical tests)
- **Y-limits**: ΔnSpikes [0 40], ΔPeak FR [0 120], Response length [0 500] (adjustable - positioning is dynamic)
- **Stats**: Wilcoxon signed-rank (CS vs US, CS vs CS+US, US vs CS+US), separate per region
- **Significance lines** (2025-11-05 optimization):
  - Level 1 (neighbor comparisons): 8% of y-range above data max, narrower lines [x±0.4]
  - Level 2 (two-side comparisons): Fixed spacing above Level 1 (ΔnSpikes: 4 units, ΔPeak: 12 units, Response length: 50 units)
  - Stars positioned 5% of y-range above lines with `VerticalAlignment: 'middle'` for optimal spacing
  - Dynamic positioning: scales automatically with any ylim changes via `y_range = curr_ylim(2) - curr_ylim(1)`
- **Fonts**: All 10pt (titles, axes, labels)

### BAfc_figure_3_supp_PN_IN.m (figures/figure_3/, 2025-11-07)
Supplementary figure: LA cell type comparison (PN vs IN). Same structure as BAfc_figure_3 but separates PNs and INs.

**Layout**: 6×2 grid (1000×1000px), identical to main Figure 3
- Rows 1-3: Heatmaps and lineplots (50% of figure)
  - Row 1: LA-PN (6 panels)
  - Row 2: LA-IN (6 panels)
  - Columns 1-3: CS, US, CS+US heatmaps
  - Columns 4-6: CS, US, CS+US lineplots (3 clusters stacked)
- Rows 4-6: Bar charts (2 rows × 3 metrics, 50% of figure)
  - **US-selective and Multisensory only** (CS-selective excluded)
  - Row 1: US-selective
  - Row 2: Multisensory

**Clustering**: Same as main Figure 3 (z≥2 excitation, 50% FR drop inhibition, clusters 1-3 responsive)

**Bar Charts**: ΔnSpikes, ΔPeak FR, Response length
- **Positions**: PN [2 3 4], IN [7 8 9]
- **Colors**: PN darker (0.7×), IN lighter (+0.5 blend)
- **Y-limits**: ΔnSpikes [0 100], ΔPeak FR [0 300], Response length [0 700]
- **Stats**: Wilcoxon signed-rank within each cell type (CS vs US, CS vs CS+US, US vs CS+US)
- **Significance lines**: Same 2-level system as main Figure 3

## Clustering Analysis (BAfc_claude_code/)

### comparison_heatmap_004.m & 005.m
CS vs US heatmaps with PNs/INs clustered separately (PNs top, INs bottom).

**004**: Two methods: (1) K-means (5 clusters, ±0.5s ROI, sorted by onset→offset→distance), (2) Manual (CS-sel, US-sel, Multi, Inhibited, Non-resp). Excitation z≥2, Inhibition z≤-2. 99th percentile limits, smoothing 101, latency detection ≥20ms sustained.

**005**: Rank score sorting (`onset + α*duration`, α=0.8). Inhibition: FR-based (30% drop, Hz not z-score). Clusters ordered 1→5 (inhibited at bottom). INs: responsive (mean rank score) vs non-responsive (FR descending), gray separator. Excitation z≥2, clim 95th percentile, softer colors (CS-sel 0.8 0.2 0.2, US-sel 0.2 0.4 0.8, Multi 0.6 0.2 0.6).

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

## Monosynaptic Response Analysis - Clean Version (uni_monosyn_opto/, 2025-10-28)

**Bug Fix**: Original `BAfc_identify_responsive_neurons.m` hardcoded LA/BA only, excluded Astria/CeA when `params.bR='All'`.

**New Files**: `BAfc_monosyn_params.m` (centralized params, `brain_regions` not `bR`, supports 'All'/string/cell array), `BAfc_identify_responsive_neurons_v2.m` (fixed region filtering, returns `analyzed_regions`), `BAfc_plot_monosynaptic_heatmaps.m` (plots regions×cell types×stimuli, 3 categories: CS-only/US-only/Both sorted by latency, z-score [-5.5,13], -0.1 to +0.1s), `example_monosyn_analysis.m` (5 examples), `README.md`.

**Usage**: Load data → `params=BAfc_monosyn_params()` → `results=BAfc_identify_responsive_neurons_v2(cell_metrics,ttl,params)` → `BAfc_plot_monosynaptic_heatmaps()`.

**BAfc_plot_monosynaptic_heatmaps**: 8 analysis sections (convergence, integration ±20%, temporal precision, latency comparisons, response magnitude, selectivity (CS-US)/(CS+US) ±0.2, latency distributions, convergence/integration grid). 3 figures + console stats. Dynamic colors for any regions.

## Publication Figures (figures/figure_4/)

### BAfc_figure_4.m (2025-10-29)

Comprehensive monosynaptic response figure with heatmaps, bar charts, and latency analysis.

**Layout**: 4×5 grid (1500×1200px)
- **Heatmaps**: Rows 1-4, Columns 1-3
  - LA: Row 1 (spanning 2 rows), tiles 1, 2, 3
  - Astria: Row 3 (spanning 2 rows), tiles 11, 12, 13
  - Shows PNs and INs separately (PNs top, INs bottom with white separator line)
- **Bar charts**: Rows 1-3, Columns 4-5 (ΔMean FR, ΔPeak FR)
- **Latency boxplots**: Row 4, Columns 4-5 (spanning 2 columns)

**Clustering & Sorting**:
- **Categories**: CS-selective, US-selective, Multisensory
- **Sorting by rank score**: `onset_latency + g.alpha * duration`
  - CS-selective: sorted by CS rank score
  - US-selective: sorted by US rank score
  - Multisensory: sorted by CS+US rank score
- **Offset handling**: If no offset detected, uses end of response window (`g.monosyn_window`)
- **Separate PN/IN clustering**: LA shows both cell types on heatmaps (white LineWidth 1 separator)

**Artifact Indication**:
- **Gray shaded overlay**: 0-10ms on all heatmaps (30% transparency)
- **Lightning bolt symbol** (⚡): Yellow, FontSize 16, only on US and CS+US heatmaps at top

**Bar Charts** (LA-PN and Astria, no INs):
- **Metrics**: ΔMean FR (ylim [0 200]), ΔPeak FR (ylim [0 300])
- **Window**: 12ms to `g.monosyn_window` (adjustable, default 25ms)
- **Positions**: LA-PN [2 3 4], Astria [7 8 9]
- **Colors**: LA-PN 0.7× cluster, Astria +0.5 blend
- **Stats**: Wilcoxon signed-rank (CS vs US, US vs CS+US, CS vs CS+US)
- **Labels**: Cluster names only (left column)

**Latency Boxplots** (ylim [0 50]):
- **4 pairs**: LA CS, Astria CS, LA US, Astria US (Selective vs Multisensory)
- **Scatter**: Jittered (gray=selective, purple=multisensory)
- **Stats**: Ranksum per pair, shows stars or "n.s."

**Colorbar**: Left of Astria heatmap, manual axes (0.010×0.15, [0.025, 0.51]), `g.colors.Heatmap`

**Parameters**: `g.monosyn_window=0.025` (adjustable), `g.bin_time=0.001`, `g.smoothvalue=7`, `g.plotwin=[0.05 0.05]`, `g.onset_threshold=5`, `g.min_consec_bins=1`, `g.alpha=0.0`, `g.use_two_rule=false`, `g.zscore_threshold_one_rule=5`

**Helper**: `compute_onset_offset_latency()` - if no offset, uses `g.monosyn_window` end

## Publication Figures (figures/figure_5/)

### BAfc_figure_5.m (2025-10-30)
Optogenetic manipulation effects on monosynaptic responses. Compares CS ± light and US ± light.

**TTLs**: `triptest_sound_only`, `triptest_sound_only_light`, `triptest_shocks_only`, `triptest_shocks_only_light`

**Detection**: Identifies monosynaptically responsive neurons in each condition separately, then compares light vs no-light effects.

**Statistics**: Wilcoxon signed-rank test on trial-by-trial spike counts in response window (12ms to `g.monosyn_window`).

**Output**: Compatible with `BAfc_monosyn_raster_ui` for visual inspection. Detailed summary with cellIDs for each category (Increased/Decreased/Unchanged).

### BAfc_figure_5_v2.m (2025-10-30, updated 2025-10-30)
Optogenetic effects on CS+US responses with fine-grained multi-window statistical testing.

**TTLs**: `triptest_both`, `triptest_both_light` (or `triptest_sound_only`, `triptest_shocks_only` variants)

**Key Feature - Fine-Grained Multi-Window Testing**: Tests for optogenetic effects in 28 time windows with 1ms increments:
- **Windows**: 12-13ms, 12-14ms, 12-15ms, ..., 12-40ms (28 windows total)
- **Artifact exclusion**: 12ms (hardcoded for shock stimuli)
- **Window generation**: `window_ends = 0.013:0.001:0.040` (1ms steps)

**Detection Strategy**:
1. Identify neurons responsive to EITHER condition (no-light OR light) using `g.monosyn_window`
2. For each responsive neuron, test all 28 windows separately
3. Classify as **Increased** if ANY window shows p<0.05 with positive change
4. Classify as **Decreased** if ANY window shows p<0.05 with negative change (no increased windows)
5. Classify as **Unchanged** if no windows reach significance

**Rationale**: Fine-grained temporal resolution captures precise timing of optogenetic effects. Different neurons may show effects at different latencies (e.g., early 12-20ms vs late 25-40ms). Testing multiple overlapping windows ensures no effects are missed due to specific window choice.

**Statistics**: Wilcoxon signed-rank test (paired, non-parametric) on trial-by-trial spike counts for each window independently.

**Visualization**: Pie charts and bar charts showing proportions (Increased/Decreased/Unchanged) for LA and Astria regions. Colors: red (increased), blue (decreased), gray (unchanged).

**Output Format**: For each neuron, reports minimum p-value and corresponding window:
```
neuron_idx=123, cellID=45, animal=MD298_001, min_p=0.0023 (12-18ms)
```

**UI Compatibility**: Outputs `g_ui` and `monosyn_results` structures for `BAfc_monosyn_raster_ui(g_ui, monosyn_results, ttl)`. Display window adjustable via `g_ui.params.pre_time_short` and `g_ui.params.post_time_short` (default 0.1s = ±100ms).

**Key Parameters**:
- `g.monosyn_window`: Detection window (e.g., 0.05s = 0-50ms)
- `g.zscore_threshold_rule1`: 3 (Rule 1 z-score threshold, lowered from 5 for CS+US)
- `g.use_two_rule`: true (two-rule detection: Rule 1 OR Rule 2)
- Artifact exclusion: 12ms (hardcoded)
- Test windows: 12-13ms to 12-40ms in 1ms steps (28 windows total)

**Output Summary**: Console printout with neuron counts and cellIDs for each category per region. Example:
```
LA: 45 responsive neurons (in either condition)
  Increased (p<0.05): 12 (26.7%)
  Decreased (p<0.05): 3 (6.7%)
  Unchanged (p>=0.05): 30 (66.7%)
```

### BAfc_figure_5_supp_light_inhibited.m (2025-11-07)
Supplementary figure showing light-inhibited neurons' responses to CS and US (no-light vs light).

**Purpose**: Visualize how neurons inhibited by pre-stimulus light respond to CS and US stimuli with and without optogenetic manipulation.

**TTLs**:
- No-light: `triptest_sound_only`, `triptest_shocks_only`
- Light: `triptest_sound_only_light`, `triptest_shocks_only_light`

**Light-Inhibited Detection** (matching BAfc_figure_2.m):
- **Method**: Wilcoxon ranksum test + FR drop threshold
- **Criteria**: p<0.05 (ranksum, right-tailed) AND FR drop ≥50%
- **Windows**: Recent (−0.5 to 0s) vs Baseline (−5 to −0.5s)
- Tests ALL LA and Astria neurons across all light conditions
- Neuron marked as inhibited if detected in ANY light condition

**Layout**: 4×4 grid (1200×1400px)
- **Row 1** (LA-CS): No-light heatmap | Light heatmap | No-light lineplot | Light lineplot
- **Row 2** (LA-US): No-light heatmap | Light heatmap | No-light lineplot | Light lineplot
- **Row 3** (Astria-CS): No-light heatmap | Light heatmap | No-light lineplot | Light lineplot
- **Row 4** (Astria-US): No-light heatmap | Light heatmap | No-light lineplot | Light lineplot

**Clustering & Sorting**:
- 5 clusters based on no-light responses: CS-sel, US-sel, Multi, Non-resp, Inhibited
- Same classification and sorting as BAfc_figure_2/3 (z≥2 excitation, 50% FR drop inhibition)
- Sorted by rank score (onset + 0.5×duration) within clusters

**Heatmaps**:
- Z-scored, 99th percentile color limits
- Black cluster lines (LineWidth 1)
- Onset/offset markers (black dots, MarkerSize 4) for both conditions
- Time window: ±2s, smoothing: Savitzky-Golay (201)

**Lineplots**:
- Population mean (Hz) for all light-inhibited neurons in that region
- Green line (no-light), Red line (light)
- Same time window as heatmaps

**Row Labels**: Brain region + "Neuron #" (left ylabel on first stimulus per region)

**Colorbar**: Manual axes, left side ([0.025, 0.40], 0.010×0.20), Z-score

**Key Parameters**:
- `g.light_inhib_window_recent = 0.5` (0-0.5s before stimulus)
- `g.light_inhib_window_baseline = 4.5` (4.5s baseline)
- `g.light_inhib_p_threshold = 0.05`
- `g.light_inhib_fr_drop = 0.5` (50% threshold)
- `g.excitation_threshold = 2`, `g.inhibition_fr_drop = 0.50`

**File Location**: `figures/figure_5/BAfc_figure_5_supp_light_inhibited.m`

### BAfc_figure_5.m - Current Production Version (2025-10-31)
Comprehensive optogenetic manipulation figure analyzing CS and US separately. Processes all 4 TTL conditions in single run.

**Layout**: 4×4 grid (1500×800px)
- **Top-left (rows 1-2, cols 1-2)**: Injection site image (`drawed_injection.png`)
- **Bottom-left (rows 3-4, cols 1-2)**: Example neuron rasters (3×4 nested)
- **Right side (rows 1-4, cols 3-4)**: Population analysis (3×4 nested, tight spacing)

**TTLs Analyzed**: All 4 conditions in single run
- CS: `triptest_sound_only` vs `triptest_sound_only_light`
- US: `triptest_shocks_only` vs `triptest_shocks_only_light`

**Results Storage**: `results_all{region, stimulus}` where stimulus: 1=CS, 2=US

**Multi-Window Statistical Testing**: 28 windows tested (12-13ms to 12-40ms, 1ms increments)
- Wilcoxon signed-rank test on trial-by-trial spike counts per window
- Enhanced: ANY window p<0.05 with positive change
- Decreased: ANY window p<0.05 with negative change (no increased windows)
- Unchanged: NO windows p<0.05

**Example Neurons** (Bottom-left rasters, 3 rows × 4 columns):
- Row 1: No-light rasters (4 neurons: MD309_001/20 CS, MD307_001/13 US, MD318_001/46 CS, MD317_001/43 US)
- Row 2: Light rasters (same neurons, light conditions)
- Row 3: FR line plots (black=no-light, red=light, ylim [-10 300])

Raster parameters:
- Time window: ±50ms, xlim [-50 50]
- MarkerSize 6, red vertical line (LineWidth 2) at time 0
- Y-labels: First column only (Trial/Trial/FR(Hz))
- X-labels: Bottom row only (Time (ms))

**Population Analysis** (Right side, 3 rows × 4 columns, tight spacing):

Row 1 - Pie charts (LA-CS, LA-US, Astria-CS, Astria-US):
- Enhanced (red 0.8 0.2 0.2): Increased with light
- Non-enhanced (gray 0.7 0.7 0.7): Decreased + Unchanged combined
- Data stored separately (n_increased, n_decreased, n_unchanged) but visualized as Enhanced vs Non-enhanced

Row 2 - Slope graphs (ΔPeak FR trajectories):
- Individual enhanced neurons: Red lines (alpha 0.3) connecting no-light to light
- Gray scatter (position 1): No-light ΔPeak FR
- Red scatter (position 2): Light ΔPeak FR
- Black thick line (LineWidth 3): Population mean trajectory
- Y-label: Only first panel (LA-CS) shows "\DeltaPeak FR (Hz)"
- ΔPeak FR calculated in neuron-specific optimal window (12ms to their most significant window_end)

Row 3 - Population FR plots (mean ± SEM for all enhanced neurons):
- Black line + gray shading (alpha 0.3): No-light condition
- Red line + red shading (alpha 0.3): Light condition
- Red vertical line: Stimulus onset
- Y-label: Only first panel (LA-CS) shows "FR (Hz)"
- ylim [-10 300], xlim [-50 50]

**Detection Parameters**:
- Two-rule system: Rule 1 (z≥3, p≥0.25) OR Rule 2 (z≥10, p≥0.1)
- `g.monosyn_window`: 0.05s (0-50ms)
- `g.zscore_threshold_rule1`: 3 (lowered from 5 for better sensitivity)
- `g.use_two_rule`: true
- Artifact exclusion: 12ms (hardcoded)
- Smoothing: Savitzky-Golay (window 7)
- Bin time: 1ms

**Spacing Optimization** (2025-10-31):
- Right nested layout: `'TileSpacing', 'tight'` (changed from 'compact')
- Y-axis labels minimized: Only first panel in each row shows ylabel
- Slope graphs: `YTickLabel` cleared for panels after first
- Population FR plots: `YTickLabel` cleared for panels after first

**Console Output**: Detailed summary per region/stimulus with neuron indices, cellIDs, animals, minimum p-value, and corresponding window for each category (Increased/Decreased/Unchanged).

## Figure Legends (2025-10-31)

Comprehensive figure legends created for all 5 publication figures in `figure_legends.txt`:

**Figure 1**: Neuron characterization across 4 brain regions (LA/BA/Astria/CeA)
- Row 1: Experimental setup, example traces
- Row 2: ISI distributions (PN/IN), spike features (trough-to-peak vs FR with L-shape), normalized waveforms
- Row 3-4: Firing rate and burst index distributions per region
- 4×4 layout, 1500×1500px, A4-optimized

**Figure 2**: Regional CS/US response comparison (4 regions × 5 columns, 1500×3000px)
- Cols 1-2: Nested z-score heatmaps (CS/US, 99th percentile)
- Col 3: Stacked bar (5 cluster proportions)
- Cols 4-5: Firing rate line plots (Hz) per cluster
- 5 clusters: CS-sel, US-sel, Multi, Non-resp, Inhibited
- Savitzky-Golay smoothing (201)

**Figure 3**: CS+US convergence in LA/Astria responsive neurons (7×2 grid, 1500×3000px)
- Rows 1-2: LA (3 nested heatmaps + lineplots)
- Rows 3-4: Astria (3 nested heatmaps + lineplots)
- Rows 5-7: ΔMean FR and ΔPeak FR bar charts (3×2 nested)
- Wilcoxon signed-rank stats, mean±SEM

**Figure 4**: Monosynaptic responses (4×5 grid, 1500×1200px)
- Cols 1-3: Heatmaps (LA rows 1-2, Astria rows 3-4) showing PNs/INs, CS/US/CS+US
- Cols 4-5 rows 1-3: ΔMean/ΔPeak FR bar charts
- Cols 4-5 row 4: Latency boxplots (Selective vs Multisensory)
- Artifact overlay 0-10ms, lightning bolt on US/CS+US
- Z-score [-5.5, 13], window ±50ms

**Figure 5**: Optogenetic manipulation (4×4 grid, 1500×800px)
- Top-left: Injection site
- Bottom-left: 4 example neurons (3 rows: no-light rasters, light rasters, FR plots)
- Right: 3×4 nested (pie charts, slope graphs, population FR)
- Multi-window testing (28 windows, 1ms increments)
- Enhanced vs Non-enhanced categories

All legends include panel descriptions, technical parameters, statistical methods, color schemes, dimensions, and key analysis details. File location: `C:\Users\dmagyar\Documents\data_analysis\BAfc_claude_code\figure_legends.txt`

## Supplementary Figures (figures/supplementary_stability/)

### BAfc_supp_stability.m (2025-11-04, updated 2025-11-06)

Response stability analysis showing neuronal responses remained stable throughout experiment.

**Purpose**: Demonstrate that neuronal responses did not change over time during the recording session.

**Method**:
- Identifies responsive neurons (clusters 1-3: CS-selective, US-selective, Multisensory) using same criteria as BAfc_figure_3
- Divides trials into 10 blocks per neuron (adaptive block size based on total trial count)
- Calculates ΔFR (response FR - baseline FR) for each block
  - Baseline: -5 to 0s window (all trials)
  - Response: 0 to 1s window (per block)
- Mean ± SEM across neurons for each block

**Layout**: 1×3 grid (1500×500px), 'tight' TileSpacing
- Panel 1: CS responses
- Panel 2: US responses
- Panel 3: CS+US responses

**Brain Regions**:
- LA (red line): All 3 stimuli
- BA (green line): US and CS+US only (no CS)
- Astria (blue line): All 3 stimuli

**Plot Details**:
- X-axis: Block # (1-10)
- Y-axis: ΔFiring Rate (Hz), ylim [0 30]
- Line: Mean ΔFR (LineWidth 3)
- Shaded area: SEM (FaceAlpha 0.3, excluded from legend via 'HandleVisibility', 'off')
- Legend: CS+US panel only, shows region and neuron count

**Key Parameters**:
- `g.excitation_threshold = 2` (z-score)
- `g.inhibition_fr_drop = 0.50` (50% FR drop)
- `g.test_time = 1` (response window 0-1s)
- `g.pre_time = 5` (baseline window -5 to 0s)
- `g.bin_time = 0.001` (1ms bins)

**File Location**: `figures/supplementary_stability/BAfc_supp_stability.m`
