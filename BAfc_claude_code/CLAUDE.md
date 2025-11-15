# CLAUDE.md

## Communication Style
**IMPORTANT**: Minimal text unless requested. Save tokens. Code changes and brief confirmations only.

## Project Overview
MATLAB scripts for basal amygdala (BA) fear conditioning neuronal data with optogenetics. K-means clustering on PSTH data.

## Directory Structure
- **BAfc_claude_code/**: Clustering analysis scripts
- **neuron_analysis/BA_fear_cond_project/**: Core BAfc helper functions
- **Data**: `C:\Users\dmagyar\Desktop\BA_fear_cond\[recording_name]\kilosort25preprocess\`

## Core Functions (../neuron_analysis/BA_fear_cond_project/)

**BAfc_load_neurons()**: Loads cell_metrics from sessions (recordings, ttl, mainFolder). Animal ID: `MD###_###` from `MD###_###_kilosort`

**BAfc_psth_spx()**: Returns `psth_spx` (neurons × bins), `preAP_norm/postAP_norm` (trial timestamps), `preAP_bin/postAP_bin` (binned per trial)

**BAfc_putative_cellTypes()**: Excitatory/inhibitory classification

**BAfc_colors()**: Standard color scheme

### TTL Event Types
- `triptest_sound_only` / `triptest_sound_only_light`: CS ± light
- `triptest_shocks_only` / `triptest_shocks_only_light`: US ± light
- `triptest_both` / `triptest_both_light`: CS+US ± light

## Analysis Methods

### Clustering
- 5 clusters: CS-sel (red), US-sel (blue), Multi (purple), Non-resp (gray), Inhibited (green)
- Excitation: z≥2, Inhibition: 50% FR drop (two-rule: 0-0.2s OR 0-1.0s)
- Rank score sorting: `onset + α*duration`
- K-means: 20 replicates, squared Euclidean, MaxIter 500

### Responsiveness Detection (Two-Rule System)
Rule 1 (z≥3, p≥0.25) OR Rule 2 (z≥10, p≥0.1). Monosynaptic: 0-25ms window, artifact exclusion 0-12ms.

### Light-Inhibited Detection
Wilcoxon ranksum + 50% FR drop (recent −0.5-0s vs baseline −5 to −0.5s, p<0.05)

### Savitzky-Golay Filter Correction (2025-11-12)
Corrects temporal delay by shifting forward `filter_delay = floor(smoothvalue/2)` bins. Applied to all figures.

## Publication Figures

### Figure Labeling Style
- Labels: 14pt bold, positioned outside panels (annotation textbox)
- Main figure panels: A, B, C, D...
- Supplementary panels: A, B, C...
- Position via normalized coordinates, avoid title overlap

### BAfc_figure_1.m (1000×1000px, 4×4 grid)
**A**: Experimental setup (empty panel, manual insertion), **B**: Example traces (CS/US, 2 cols), **C**: ISI (PN/IN, 2 cols), **D**: Spike features + waveforms (2 cols), **E**: FR distributions (nested 1×4), **F**: Burst indices (nested 1×4)
Rows 3-4 use nested tiledlayouts with tight spacing, ylabel 'Count' on left only

### BAfc_figure_2_example_rasters.m (1000×300px, 2×2 grid)
**A**: 4 example neurons (CS-sel, US-sel, Multi, Inhibited), each with CS/US rasters (nested 1×2)

### BAfc_figure_2.m (1000×700px, 4 regions × 5 cols)
**B**: Main heatmaps + lineplots. Cols 1-2: CS+US heatmaps (z-score, 99th percentile), Col 3: Stacked bars, Cols 4-5: CS/US lineplots (Hz)
**C**: Chi-square p-values (1×4 layout, tile 1)
**D**: US metrics bar graphs (nested 1×3, tiles 2-4, title "Comparison of US-evoked excitatory responses")
**Supp A-C**: Permutation test, Contingency table (spans 2 tiles), Cramér's V (spans 2 tiles). 1×5 layout.

### BAfc_figure_3.m (1000×1000px, 5×7 grid)
Rows 1-4: LA/AStria heatmaps (cols 1-3: CS/US/CS+US, "LA neurons", "AStria neurons"), lineplots (cols 4-6: 3 clusters stacked), ΔFR bars (col 7)
Row 5: Pie charts (left: LA/AStria proportions), Region comparison bars (right: CS/US/CS+US, grey circles with jitter, y-axis to 95th percentile)
**Kruskal-Wallis gating**: Main figure bar charts only show signrank stars if KW p<0.05
**Supp figure** (6×3, 1000×1200px): Chi-square (row 1, 1×3), Metrics bars (rows 2-4, 3 clusters × 3 metrics), KW tests (rows 5-6, LA/AStria × 3 clusters). Post-hoc only shown if KW p<0.05
**Supp PN_IN**: LA PN vs IN comparison (6×2 grid: rows 1-3 heatmaps/lineplots, rows 4-6 bar charts)

### BAfc_figure_4.m (1200×800px, 4×6 grid)
Cols 1-3: Monosynaptic heatmaps (LA rows 1-2, AStria rows 3-4, CS/US/CS+US, "LA neurons", "AStria neurons")
Col 4: ΔPeak FR bars (6×1 nested, 3 clusters each region)
Cols 5-6: Pie charts (row 1, FontSize 14), Region comparison (row 2, grey circles with jitter, y-axis to 95th percentile), Latency comparison (row 3)
α=0.0 for rank score sorting, `g.monosyn_window=0.025`
**Kruskal-Wallis gating**: Main figure bar charts only show signrank stars if KW p<0.05
**Supp figure** (3×1, 1000×900px): Selective vs multisensory latency (row 1), Chi-square test (row 2, 1×3), KW tests (row 3, 2×3 nested, LA/AStria × 3 clusters). Post-hoc only shown if KW p<0.05

### BAfc_figure_5.m (1000×500px, 2×4 grid) - Updated 2025-11-14
**A**: Example rasters (rows 1-2), **B**: FR lineplots (row 3), **C**: Pie charts (row 1 right), **D**: Spaghetti plots (row 2 right)
Left: Example neurons (3×4: no-light rasters, light rasters, FR plots)
Right: Population (2×4: pie charts, spaghetti plots with mean spike count)
**Testing method** (configurable via `g.testing_method`):
- `'single_window'` (DEFAULT, reviewer-safe): Test predefined window (12-50ms), no multiple comparisons issue
- `'multi_window'` (exploratory): Test 12ms to `g.monosyn_window` in 1ms steps, ANY window p<0.05 with positive change (WARNING: 86% false positive rate)
Window parameters: `g.test_window_start = 0.012`, `g.test_window_end = 0.050`
Example neurons: MD309_001 cell 20 (LA CS/US), MD318_001 cell 46 (AStria CS), MD317_001 cell 43 (AStria US)
**Supp light_inhibited** (5×4 grid, 1000×1000px): **A**: Example rasters (row 1), **B**: LA heatmaps (row 2), **C**: LA lineplots (row 3), **D**: AStria heatmaps (row 4), **E**: AStria lineplots (row 5)

### BAfc_supp_stability.m (1500×500px, 1×3 grid)
Response stability across 10 trial blocks. ΔFR = response FR (0-1s) - baseline FR (−5-0s). LA/BA/Astria, mean±SEM.

## Key Parameters
- `g.bin_time = 0.001` (1ms bins)
- `g.smoothvalue = 201` (figures 2-3) or 5 (figures 4-5)
- `g.excitation_threshold = 2`, `g.inhibition_fr_drop = 0.50`
- `g.monosyn_window = 0.025` (monosynaptic: 0-25ms)
- `g.zscore_threshold_rule1 = 3`, `g.zscore_threshold_rule2 = 10`
- `g.prob_threshold_rule1 = 0.25`, `g.prob_threshold_rule2 = 0.1`
- Artifact exclusion: 0-12ms (shocks)
- Fonts: 10pt throughout (section titles 12pt)

## Statistical Methods (2025-11-14)
### Kruskal-Wallis Gating (Figures 3-4)
- Main figure bar charts: Perform KW test first, only show signrank stars if KW p<0.05
- Supplementary KW panels: Only display post-hoc p-values if KW p<0.05
- Prevents multiple comparisons inflation in bar chart visualizations

### Figure 5 Optogenetic Testing
- Default: Single pre-defined window (12-50ms), statistically sound, reviewer-safe
- Alternative: Multi-window exploratory testing (requires justification, high false positive rate)

## Figure Legends
All figures have accompanying `.docx` legends created via Python scripts:
- `create_figure_X_legends.py` in each figure directory
- Format: Heading 1 titles, Normal paragraphs with bold panel labels (A), 10pt font
- Regenerate: `cd figures/figure_X && python create_figure_X_legends.py`

## Common Modifications
**Clusters**: `g.clustnum = 5;`
**Region**: `g.bR = 'LA';` (or 'BA', 'Astria', 'CeA')
**Stimulus**: Edit `ttl` variable
**Time window**: `g.roi_pca = round((g.pre_time-0.5)/g.bin_time+1:(g.pre_time+0.5)/g.bin_time);`
