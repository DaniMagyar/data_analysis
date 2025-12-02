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
**C**: Chi-square p-values (1×3 layout, tile 1)
**D**: US metrics bar graphs (nested 1×2, tiles 2-3, title "Comparison of US-evoked excitatory responses"). Shows ΔFR and Response length only (ΔnSpikes removed 2025-11-26)
**Supp A-C**: Permutation test, Contingency table (spans 2 tiles), Cramér's V (spans 2 tiles). 1×5 layout. Shows ΔFR and Response length (ΔnSpikes removed 2025-11-26). K-W post-hoc only shown if K-W p<0.05.

### BAfc_figure_3.m (1000×1000px, 5×7 grid)
Rows 1-4: LA/AStria heatmaps (cols 1-3: CS/US/CS+US, "LA neurons", "AStria neurons"), lineplots (cols 4-6: 3 clusters stacked), ΔFR bars (col 7)
Row 5: Pie charts (left: LA/AStria proportions), Region comparison bars (right: CS/US/CS+US, grey circles with jitter, y-axis to 95th percentile)
**Kruskal-Wallis gating**: Main figure bar charts only show signrank stars if KW p<0.05
**ΔPeak FR calculation** (2025-11-26): Uses fixed 0-1s response window for ALL neurons regardless of latency detection (ensures non-zero values for sub-threshold responses)
**Supp figure** (6×3, 1000×1200px): Chi-square (row 1, 1×3), Metrics bars (rows 2-4, 3 clusters × 2 metrics - ΔFR and Response length only, ΔnSpikes removed 2025-11-26), KW tests (rows 5-6, LA/AStria × 3 clusters). Post-hoc only shown if KW p<0.05
**Supp PN_IN**: LA PN vs IN comparison (6×2 grid: rows 1-3 heatmaps/lineplots, rows 4-6 bar charts)

### BAfc_figure_4.m (1200×800px, 4×6 grid)
Cols 1-3: Monosynaptic heatmaps (LA rows 1-2, AStria rows 3-4, CS/US/CS+US, "LA neurons", "AStria neurons")
Col 4: ΔPeak FR bars (6×1 nested, 3 clusters each region)
Cols 5-6: Pie charts (row 1, FontSize 14), Region comparison (row 2, grey circles with jitter, y-axis to 95th percentile), **G**: Latency comparison (row 3, Wilcoxon rank-sum LA vs AStria)
α=0.0 for rank score sorting, `g.monosyn_window=0.025`
**Kruskal-Wallis gating**: Main figure bar charts only show signrank stars if KW p<0.05
**Supp figure removed** (2025-12-01): All supplementary data now in Excel export only
**Excel export** (2025-12-01): `figure_4_data.xlsx` with sheets:
  - `PanelsBD_DeltaFR`: Bar chart data (monosynaptic 12-25ms window) with Friedman tests
  - `PanelE_PieCharts`: Cluster proportions with Chi-square tests
  - `PanelF_RegionComparison`: LA vs AStria comparison (monosynaptic 12-25ms)
  - `PanelG_LatencyComparison`: Onset latency statistics with Wilcoxon rank-sum tests
  - `RawData_DeltaFR_OnsetLat`: Individual neuron values (Local #, Global Index, Animal ID, ΔPeak FR CS/US/CS+US, Onset Lat CS/US/CS+US)
**NOTE**: Low smoothvalue (5) with 1ms bins can cause extremely high ΔPeak FR values (800+ Hz) in 12-25ms window. Consider increasing smoothvalue to 15-25 for monosynaptic analysis.

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
**Stats export** (2025-11-27): `figure_5_supp_stats.txt` includes PV silencing effect (pooled CS+US, -0.5-0s baseline comparison, Wilcoxon signed-rank), firing rates by stimulus, latencies, cluster distributions

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

## Data Export to Excel (2025-12-01)
All figures export data to `figure_X_data.xlsx` for verification and supplementary tables:

### Implementation Pattern (see BAfc_figure_3.m as reference)
At end of figure script, call: `export_figureX_to_excel_simple(...)`

### Export Function Structure
```matlab
function export_figureX_to_excel_simple(results_all, kw_data_storage, kw_results, ...)
    if exist(output_filename, 'file')
        delete(output_filename);
    end

    % Create sheets with writecell()
    % ...

    writecell(sheet_data, output_filename, 'Sheet', 'SheetName');
end
```

### Required Sheets (Figure 3 example)
1. **Summary sheets** (3-4 sheets):
   - `PanelsCF_DeltaFR_RespLen`: Bar chart data with ΔPeak FR and Response Length side-by-side
     - Columns 1-5: ΔPeak FR (Mean, SEM, Median, SD in Hz)
     - Column 6: Empty separator
     - Columns 7-10: Response Length (Mean, SEM, Median, SD in ms)
     - Statistical tests displayed horizontally (Friedman + post-hoc for both metrics side-by-side)
   - `PanelG_PieCharts`: Cluster proportions with Chi-square tests
   - `PanelH_RegionComparison`: LA vs AStria comparison with Wilcoxon rank-sum tests
   - All sheets: descriptive stats + p-values with significance markers (*, **, ***, n.s.)

2. **Raw data sheets** (2-3 sheets):
   - `RawData_DeltaFR_RespLen`: Individual neuron values for both ΔPeak FR and Response Length
     - Columns 1-6: `Local #`, `Global Index`, `Animal ID`, ΔPeak FR CS/US/CS+US (Hz)
     - Column 7: Empty separator
     - Columns 8-10: Response Length CS/US/CS+US (ms)
     - Include Mean and SEM at bottom of each section for verification
   - `RawData_RegionComparison`: Individual neuron values for Panel H
   - `Lineplots_NeuronIndices`: Global indices + animal IDs only (NOT full time-series data)
   - Columns: `Local #`, `Global Index`, `Animal ID`, data values

### Figure 4 Excel Export Pattern
Similar to Figure 3 but for monosynaptic responses:
- `PanelsBD_DeltaFR`: ΔPeak FR only (no response length - monosynaptic)
- `PanelG_LatencyComparison`: Onset latency with Mean, SEM, Median, SD, Range
- `RawData_DeltaFR_OnsetLat`: Combined ΔFR + Onset Latency columns with Global Index and Animal ID

### Helper Function
```matlab
function sig_str = format_significance(p_val)
    if p_val < 0.001
        sig_str = '***';
    elseif p_val < 0.01
        sig_str = '**';
    elseif p_val < 0.05
        sig_str = '*';
    else
        sig_str = 'n.s.';
    end
end
```

### File Size Guidelines
- Keep file small (~40-100 KB) by exporting summaries + raw neuron values
- DO NOT export full time-series data (creates 30+ MB files)
- For lineplots: export only neuron indices, not 4000 time points per neuron

### Console Output
- Minimal messages only
- Remove verbose fprintf statements
- No "Exported..." or "EXPORTING..." messages

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
