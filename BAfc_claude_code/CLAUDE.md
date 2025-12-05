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

### Clustering (Figures 2-3)
- 5 clusters: CS-sel (red), US-sel (blue), Multi (purple), Non-resp (gray), Inhibited (green)
- Excitation: z≥1.5 (peak in 0-1s window), Inhibition: 50% FR drop (0-0.2s OR 0-1.0s)
- Rank score sorting: `onset + α*duration`

### Onset/Offset Latency Detection
- Threshold: z≥1.5 (same as excitation threshold)
- Requires 20ms consecutive bins above/below threshold
- Onset: First sustained excursion above threshold
- Offset: First sustained return below threshold after onset
- If no offset found: uses end of 1s window
- If no onset found: onset/offset = NaN

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
**NaN handling**: Burst index statistics use `'omitnan'` option for mean/std/median, and exclude NaNs before ranksum tests
**Excel export** (2025-12-03): `figure_1_data.xlsx` with sheets:
  - `Summary_SampleSizes`: Animal counts, neuron counts by region
  - `PanelD_SpikeFeatures`: TTP and firing rate statistics (LA+BA)
  - `PanelE_FiringRate`: Firing rate distributions by region
  - `PanelF_BurstIndex`: Burst index distributions by region (handles NaN values)
  - `CrossRegion_FiringRate`: Kruskal-Wallis test + post-hoc comparisons

### BAfc_figure_2_example_rasters.m (1000×300px, 2×2 grid)
**A**: 4 example neurons (CS-sel, US-sel, Multi, Inhibited), each with CS/US rasters (nested 1×2)

### BAfc_figure_2.m (1000×700px, 4 regions × 5 cols)
**B**: Main heatmaps + lineplots. Cols 1-2: CS+US heatmaps (z-score, 99th percentile), Col 3: Stacked bars, Cols 4-5: CS/US lineplots (Hz)
**C**: Chi-square p-values (1×3 layout, tile 1)
**D**: US metrics bar graphs (nested 1×2, tiles 2-3, title "Comparison of US-evoked excitatory responses"). Shows Peak FR (Hz) and Response length (ms)
**Peak FR calculation** (2025-12-03): Uses maximum firing rate in 0-1s response window (`max(psth_Hz(:, g.roi))`). No baseline subtraction - absolute peak values.
**Response length**: Onset-to-offset duration. If no offset detected (sustained response), uses end of 1s window. NaN if no onset detected.
**Supplementary figure removed** (2025-12-03): Permutation test and Cramér's V plots removed from script. Data retained in Excel export.
**Excel export** (2025-12-03): `figure_2_data.xlsx` with sheets:
  - `Summary_SampleSizes`: Animal counts and neuron counts by region
  - `PanelB_Clusters`: Contingency table and cluster proportions
  - `PanelC_ChiSquare`: Chi-square test with FDR-corrected pairwise comparisons
  - `PanelD_USMetrics`: Peak FR and Response Length statistics with Kruskal-Wallis tests
  - `ClusterLatencies`: Cluster-specific onset/offset latencies for CS and US
  - `RawData_AllNeurons`: Individual neuron values (all clusters) with Global Index, Animal ID, Cluster, CS Peak (z-score & Hz), US Peak (z-score & Hz), Response Length (ms)

### BAfc_figure_3.m (1000×1000px, 5×7 grid)
Rows 1-4: LA/AStria heatmaps (cols 1-3: CS/US/CS+US, "LA neurons", "AStria neurons"), lineplots (cols 4-6: 3 clusters stacked), FR bars (col 7)
Row 5: Pie charts (left: LA/AStria proportions), Region comparison bars (right: CS/US/CS+US, grey circles with jitter, y-axis to 95th percentile)
**Panel labels**: A (LA heatmaps), B (LA lineplots), C (LA FR bars), D (AStria heatmaps), E (AStria lineplots), F (AStria FR bars), G (pie charts), H (region comparison)
**Friedman test gating**: Main figure bar charts only show signrank stars if Friedman p<0.05
**Panel H (Region comparison)**: Uses only responsive neurons (CS-sel, US-sel, Multi clusters 1-3) for LA vs AStria comparison
**Peak FR calculation** (2025-12-03): Uses maximum firing rate in 0-1s response window (`max(psth_Hz(:, g.roi))`). No baseline subtraction - absolute peak values. Applied to CS, US, and CS+US responses.
**Response length**: Computed from onset/offset latencies for ALL neurons (not filtered). Excel export retains NaN values when no onset/offset detected.
**Heatmaps**: Display only responsive neurons (clusters 1-3); all 5 clusters used for pie charts and statistics
**Excel export** (2025-12-03): `figure_3_data.xlsx` with sheets:
  - `PanelsCF_PeakFR_RespLen`: Bar chart data with Peak FR and Response Length (Mean, SEM, Median, SD) with Friedman tests
  - `PanelG_PieCharts`: Cluster proportions with Chi-square tests
  - `PanelH_RegionComparison`: LA vs AStria comparison with Wilcoxon rank-sum tests
  - `RawData_PeakFR_RespLen`: Individual neuron values with Local #, Global Index, Animal ID, Peak FR CS/US/CS+US (Hz), Response Length CS/US/CS+US (ms)
**Supp figure** (6×3, 1000×1200px): Chi-square (row 1, 1×3), Metrics bars (rows 2-4, 3 clusters × 2 metrics - Peak FR and Response length), KW tests (rows 5-6, LA/AStria × 3 clusters). Post-hoc only shown if KW p<0.05

### BAfc_figure_4.m (1200×800px, 4×6 grid)
Cols 1-3: Monosynaptic heatmaps (LA rows 1-2, AStria rows 3-4, CS/US/CS+US, "LA neurons", "AStria neurons")
Col 4: ΔFR bars (6×1 nested, 3 clusters each region)
Cols 5-6: Pie charts (row 1, FontSize 14), Region comparison (row 2, grey circles with jitter, y-axis to 95th percentile), **G**: Latency comparison (row 3, Wilcoxon rank-sum LA vs AStria)
α=0.0 for rank score sorting, `g.monosyn_window=0.025`
**Kruskal-Wallis gating**: Main figure bar charts only show signrank stars if KW p<0.05
**Panel F (Region comparison)**: CS responses from CS-sel + Multi (clusters 1,3), US responses from US-sel + Multi (clusters 2,3), CS+US from all responsive (clusters 1,2,3). Uses absolute FR, no baseline subtraction.
**ΔFR calculation** (2025-12-02): `ΔFR = (nSpikes/trial) / window_duration`, where window_duration = 13ms (0.013s). This avoids smoothing artifacts from high-frequency binning while maintaining consistent Hz units. CS+US uses actual `triptest_both` trials, NOT pooled CS+US data.
**Supp figure removed** (2025-12-01): All supplementary data now in Excel export only
**Excel export** (2025-12-01, updated 2025-12-05): `figure_4_data.xlsx` with sheets:
  - `PanelsBD_DeltaFR`: Bar chart data (monosynaptic 12-25ms window) with Friedman tests. ΔFR calculated from spike counts per trial.
  - `PanelE_PieCharts`: Cluster proportions with Chi-square tests
  - `PanelF_RegionComparison`: LA vs AStria comparison (monosynaptic 12-25ms). CS: clusters 1,3; US: clusters 2,3; CS+US: clusters 1,2,3
  - `PanelG_LatencyComparison`: Onset latency statistics with Wilcoxon rank-sum tests
  - `RawData_DeltaFR_OnsetLat`: Individual neuron values (Local #, Global Index, Animal ID, ΔFR CS/US/CS+US, Onset Lat CS/US/CS+US)

### BAfc_figure_5.m (1000×500px, 2×4 grid) - Updated 2025-12-05
**A**: Example rasters (rows 1-2), **B**: FR lineplots (row 3), **C**: Pie charts (row 1 right), **D**: Spaghetti plots (row 2 right)
Left: Example neurons (3×4: no-light rasters, light rasters, FR plots)
Right: Population (2×4: pie charts, spaghetti plots with mean spike count)
**Testing method**: Single pre-defined window (10-50ms), statistically sound, no multiple comparisons issue
Window parameters: `g.artifact_end = 0.010`, `g.monosyn_window = 0.050`
Example neurons: MD309_001 cell 20 (LA US ↑), MD317_001 cell 43 (AStria US ↑), MD307_001 cell 16 (LA US ↓), MD315_001 cell 15 (AStria US ↓)
**Pie chart label overlap fix** (2025-12-03): Manual position adjustments applied to prevent overlap between adjacent pie charts. LA US 34% moved down (-0.1), 54% moved down (-0.3). AStria CS 16% moved right-up (0.05, 0.1), 66% moved down (-0.15).
**Excel export**: `figure_5_data.xlsx` with sheets:
  - `PanelE_PieCharts`: Proportions of enhanced/decreased/unchanged neurons with test window
  - `PanelF_SpaghetttiPlots`: FR statistics (Mean, SEM, Median, SD) for all three classifications (increased, decreased, unchanged) with paired t-tests
  - `RawData_Classifications`: Individual neuron classifications with p-values and FR values (no-light and with-light)
  - `ChiSquare_RegionComparisons`: LA vs AStria chi-square tests for CS and US

### Supplementary Figure 2 (supplementary_figure_2/)
Response stability across trial blocks (1000×500px, 1×2 grid: CS/US only). Uses first 50 trials (10 blocks, 5 trials/block). Shows simple FR (Hz) in 0-1s response window, not ΔFR. Cluster assignments loaded from `figure_2_data.xlsx` (RawData_AllNeurons sheet). CS panel: CS-sel + Multi neurons. US panel: US-sel + Multi neurons. LA/BA/Astria, mean±SEM.
**Script**: `BAfc_supplementary_figure_2_v2.m`
**Excel export**: `supplementary_figure_2_data.xlsx` with sheets:
  - `Summary`: Parameters (5 trials/block, 50 max trials, 10 blocks), neuron counts by region
  - `Block_Statistics`: Mean ± SEM FR for each block
  - `Neurons_Per_Block`: Region-Stimulus | Block | Global Index (column format)
  - `Neuron_List`: All neurons with Global Index, Animal ID, Cell ID, Region, Cell Type

### Supplementary Figure 3 (supplementary_figure_3/)
LA PN vs IN comparison (6×2 grid: rows 1-3 heatmaps/lineplots, rows 4-6 bar charts)

### Supplementary Figure 4 (supplementary_figure_4/)
Light-inhibited neurons (5×4 grid, 1000×1000px): **A**: Example rasters (row 1), **B**: LA heatmaps (row 2), **C**: LA lineplots (row 3), **D**: AStria heatmaps (row 4), **E**: AStria lineplots (row 5)
**Excel export**: `supplementary_figure_4_data.xlsx` with sheets:
  - `Summary`: Sample sizes, detection criteria, cluster distributions
  - `PV_Silencing_Effect`: Pooled CS+US, -0.5-0s baseline comparison (Mean, SEM, Median) with Wilcoxon signed-rank test
  - `Firing_Rate_Analysis`: Mean±SEM±Median±SD FR in three time windows (-1 to -0.5s baseline, -0.5 to 0s illumination, 0 to 0.5s response) for CS/US ± light (4 conditions × 3 windows). Statistical comparisons: (1) Baseline vs illumination effect (no-light trials pooled), (2) CS response no-light vs light, (3) US response no-light vs light
  - `Complete_Neuron_Lists`: Individual neuron details (position, animal, cell ID, region, cell type, cluster, p-value, FR drop) + all calculated FR values for 4 conditions × 3 time windows (12 columns). Includes Mean and SEM rows at bottom for verification

## Key Parameters
- `g.bin_time = 0.001` (1ms bins)
- `g.smoothvalue = 201` (figures 2-3) or 5 (figures 4-5)
- `g.excitation_threshold = 1.5`, `g.inhibition_fr_drop = 0.50`
- `g.onset_threshold = g.excitation_threshold` (for latency detection)
- `g.min_consec_bins = 20` (20ms sustained response required for onset/offset detection)
- `g.monosyn_window = 0.025` (monosynaptic: 0-25ms)
- `g.roi = 0-1s` response window for peak FR and latency detection
- Artifact exclusion: 0-12ms (shocks)
- Fonts: 10pt throughout (section titles 12pt)

## Statistical Methods (2025-11-14)
### Kruskal-Wallis Gating (Figures 3-4)
- Main figure bar charts: Perform KW test first, only show signrank stars if KW p<0.05
- Supplementary KW panels: Only display post-hoc p-values if KW p<0.05
- Prevents multiple comparisons inflation in bar chart visualizations

### Figure 5 Optogenetic Testing (Updated 2025-12-05)
- Single pre-defined window (10-50ms): Statistically sound, no multiple comparisons issue
- Wilcoxon signed-rank test: Compares no-light vs light conditions for each neuron
- Classification: Increased (p<0.05, positive change), Decreased (p<0.05, negative change), Unchanged (p≥0.05 or no change)

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

### Required Sheets (Figure 2 and 3 examples)
1. **Summary sheets** (3-4 sheets):
   - `PanelsCF_PeakFR_RespLen` (Fig 3) or `PanelD_USMetrics` (Fig 2): Bar chart data with Peak FR and Response Length side-by-side
     - Columns 1-5: Peak FR (Mean, SEM, Median, SD in Hz)
     - Column 6: Empty separator
     - Columns 7-10: Response Length (Mean, SEM, Median, SD in ms)
     - Statistical tests displayed horizontally (Friedman/Kruskal-Wallis + post-hoc for both metrics)
   - `PanelG_PieCharts`: Cluster proportions with Chi-square tests
   - `PanelH_RegionComparison`: LA vs AStria comparison with Wilcoxon rank-sum tests
   - All sheets: descriptive stats + p-values with significance markers (*, **, ***, n.s.)

2. **Raw data sheets** (2-3 sheets):
   - `RawData_PeakFR_RespLen` (Fig 3) or `RawData_AllNeurons` (Fig 2): Individual neuron values
     - Figure 2: CS Peak (z-score & Hz), US Peak (z-score & Hz), Response Length (ms)
     - Figure 3: Peak FR CS/US/CS+US (Hz), Response Length CS/US/CS+US (ms)
     - Include `Local #`, `Global Index`, `Animal ID` for each neuron
     - Include Mean and SEM at bottom for verification
   - `RawData_RegionComparison`: Individual neuron values for Panel H
   - `Lineplots_NeuronIndices`: Global indices + animal IDs only (NOT full time-series data)

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
