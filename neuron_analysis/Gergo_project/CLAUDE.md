# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview
MATLAB analysis scripts for amygdala projection neuron recordings in response to unconditioned stimuli (US/shock). Focuses on comparing responses between PFC-projecting and DMS-projecting neurons (Figures 1, 4, 5) as well as LA vs BA neurons (Figures 2, 3) across different experimental conditions (anesthetized, awake, optogenetically tagged).

## Data Organization

### Data Sources
Main data folder: `C:\Users\dmagyar\Desktop\Gergo`

**Projection neuron data** (Figure 1):
- PFC projecting: `Gergo\Opto 1 Gergo mPFCbe vetito\Opto 1\`
- DMS projecting: `Gergo\Opto 2 Gergo DMSbe vetito\Opto 2\`
- PFC tagging: `Gergo\NagyGergo_Opto1_PrL\`
- DMS tagging: `Gergo\NagyGergo_Opto2_DMS\`

**Brain region data** (Figures 2-5):
- LA awake: `Gergo\LA_awake\`
- BA awake: `Gergo\BA_awake\`
- LA anesthetized: `Gergo\LA_urethane Gergo\LA_urethane\`
- BA anesthetized: `Gergo\BA_urethane Gergo\BA_urethane\`

**Probe recordings** (Figures 4-5):
- Awake/anesthetized data loaded via `BAfc_load_neurons()` with `MD313_001_kilosort` recording
- Batch data: `Gergo\AMY_shok_rest_run_Dani adatok\cell_metrics_batch.mat`
- TTL files: `TTLsKS_Batch[1-4].mat` with `shock_motor_rest` timestamps

### Data File Structure
Each recording folder contains:
- Spike files: `GR*.mat` (MClust format with `TS` field, divide by 10000 for seconds)
- TTL file: `shocKTTL.mat` with `shockTTL` field (shock timestamps)
- Events file: `*events*` (Open Ephys format for optogenetic stimulation)

## Cell Metrics Structure
Scripts build a `cell_metrics` structure with:
- `spikes.times`: Cell array of spike timestamps (seconds)
- `projection`: 'PFC' or 'DMS' (Figure 1)
- `brainRegion`: 'LA' or 'BA' (Figures 2-5)
- `general.shockTTL`: Shock timestamps per neuron
- `general.optoTTL`: Optogenetic stimulation timestamps (Figure 1 only)
- `cellID`: Neuron identifiers

## Core Helper Functions
Located in `../BA_fear_cond_project/`:

**BAfc_psth_spx()**: Generates peri-stimulus time histograms
- Usage: `psth_spx = BAfc_psth_spx('cell_metrics', g.cell_metrics, 'ttl', 'shockTTL', 'pre_time', 1, 'post_time', 1, 'bin_time', 0.005)`
- Returns: Matrix of spike counts (neurons × time bins)

**BAfc_putative_cellTypes()**: Classifies excitatory/inhibitory neurons
- Used in Figures 4-5 probe data

**BAfc_colors()**: Standard color scheme
- Returns struct with `Heatmap` colormap

**BAfc_load_neurons()**: Loads Kilosort data
- Usage: `cell_metrics = BAfc_load_neurons('mainFolder', g.mainFolder, 'recordings', {'MD313_001_kilosort'}, 'ttl', {'triptest_shocks_only'})`

**load_open_ephys_data()**: Parses Open Ephys event files
- Returns data, timestamps, info

## Analysis Scripts

### Gergo_figure_1.m - Projection Neuron Analysis
Compares PFC-projecting vs DMS-projecting neurons responding to US (shock).

**Two analysis sections**:
1. **US response** (lines 38-220): Shock-evoked activity
   - Time window: -1 to +1s, 5ms bins, smoothing window 21
   - Test window: 0-600ms post-stimulus
2. **Optotagging** (lines 221-274): Validates projection identity
   - Time window: -20 to +20ms, 1ms bins, smoothing window 3

**Outputs**: 4×2 figure layout
- Rows 1-2: US heatmaps (PFC left, DMS right) + zscore traces + firing rate traces
- Rows 3-4: Optotagging heatmaps (same layout)
- Neurons sorted by response direction → cluster → onset → offset timing

### Gergo_figure_2.m - Anesthetized LA/BA Response
Basic comparison of LA vs BA neurons under urethane anesthesia.

**Parameters**:
- Time window: -1 to +1s, 5ms bins, smoothing window 20
- Test window: 0-600ms post-shock

**Outputs**: 3×2 figure layout
- Row 1: Schematic (`mouse_drawing_anesthetized.png`)
- Row 2: Heatmaps (LA left, BA right)
- Row 3: Firing rate histograms

### Gergo_figure_3.m - Awake LA/BA Response
Same structure as Figure 2 but for awake recordings. Combines Kilosort batch data with MClust data.

### Gergo_figure_4.m - Awake LA/BA with Examples
Extended version of Figure 3 with raster plot examples.

**Layout**: 12×8 tile grid
- Tiles 1-4: Schematic (`probe_shank_neurons.png`)
- Tiles 5,21: Example rasters (Cell 530 LA, Cell 468 BA from MD313_001)
- Tiles 33-40: LA/BA heatmaps (only excited neurons)
- Tiles 65-72: LA/BA spike count bar plots

**Key feature**: Only includes excited neurons (clusters 1-2 from `Gergo_psth_sorter`)
- `n_LA_excited = sum(clusters(idx_LA) == 1) + sum(clusters(idx_LA) == 2)`
- Neurons sorted by response direction → onset only

**Parameters**:
- PSTH: 200ms window, 1ms bins, smoothing 11
- Bar plots: 2ms bins for higher resolution

### Gergo_figure_5.m - Anesthetized LA/BA with Examples
Same structure as Figure 4 but for urethane-anesthetized data. Example neuron plots commented out.

**Special handling**:
- Manually removes artifact spikes at stimulus onset: `psth_spx(:, g.pre_time/g.bin_time+1:g.pre_time/g.bin_time+2) = 0`
- Must use 5ms bins (enforced via console message)

### Gergo_psth_sorter.m - Response Classifier
Categorizes neurons by response profile using z-scored PSTH data.

**Inputs**:
- `g.roi`: Region of interest bins (typically 0-600ms post-stimulus)
- `g.exctestvalue`: Z-score threshold for excitation (typically 1)
- `g.inhtestvalue`: Z-score threshold for inhibition (typically 1)
- `g.bin_time`: Bin size for converting indices to milliseconds

**Response classification logic**:
- **Excited (respDir = -1)**: Peak z-score ≥ exctestvalue in ROI
- **Inhibited (respDir = 1)**: Peak inverted z-score ≥ inhtestvalue
- **Non-responsive (respDir = 0)**: Neither criterion met

**Temporal features**: Calculates onset/offset indices with quality control
- Response mean must be ≥ half the threshold to prevent spurious detections
- Iteratively removes trailing low-quality bins

**Cluster assignments** (based on onset/offset timing):
- **Cluster 1**: Fast excited (onset 0-50ms, duration 0-50ms)
- **Cluster 2**: Sustained/late excited (all other excited)
- **Cluster 3**: Inhibited (all timing patterns)
- **Cluster 5**: Non-responsive

**Returns**:
- `results.respDir`: Response direction per neuron
- `results.onsetIdx`: Response onset bin
- `results.offsetIdx`: Response offset bin
- `clusters`: Cluster assignment (1, 2, 3, or 5)

### Gergo_raster_plot_gui.m - Interactive Raster Viewer
GUI for exploring single-neuron rasters aligned to shock TTL.

### Gergo_firstspike_latency.m - Latency Calculator
Calculates first spike latency after stimulus onset.

## Common Parameters

**Standard PSTH settings**:
- `g.pre_time = 1` / `0.2` (shock response / probe examples)
- `g.post_time = 1` / `0.2`
- `g.bin_time = 0.005` / `0.01` / `0.001` (coarse / medium / fine)
- `g.smoothvalue = 21` / `11` / `7` (Savitzky-Golay window)

**Analysis windows**:
- `g.test_time = 0.6` / `0.1` / `0.2` (shock response duration)
- `g.roi = g.pre_time/g.bin_time+1:(g.pre_time+g.test_time)/g.bin_time` (bins)

**Response thresholds** (z-score):
- `g.exctestvalue = 1` / `3` (excitation, anesthetized / awake)
- `g.inhtestvalue = 1` (inhibition)

**Figure aesthetics**:
- `g.fontSize1 = 15` / `13` (titles)
- `g.fontSize2 = 12` (axes)
- `g.xlinewidth = 2` (stimulus marker)
- `g.clim = [min(psth_spx_zscore(:)) max(psth_spx_zscore(:))]` (auto heatmap scale)

## Typical Workflow

1. **Configure data paths**: Set `g.mainFolder` and subfolders
2. **Build cell_metrics**: Loop through folders, load spike/TTL files
3. **Generate PSTH**: Call `BAfc_psth_spx()` with desired parameters
4. **Z-score and smooth**: `zscore(psth_spx,0,2)` then `smoothdata(...,'sgolay',g.smoothvalue)`
5. **Classify responses**: Call `Gergo_psth_sorter()` to get clusters
6. **Sort neurons**: Use `sortrows()` on response direction, cluster, timing
7. **Plot results**: Heatmaps (`imagesc`), traces (`plot` with shaded SEM), rasters (`scatter`)

## Figure Generation
All scripts generate publication-ready figures with `tiledlayout` for multi-panel layouts. Use `exportgraphics(gcf, 'filename.tiff', 'Resolution', 300)` to save.

## Notes on Data Conversion
- **MClust timestamps**: Divide by 10000 to convert to seconds (`spikes.TS/10000`)
- **Open Ephys events**: Keep odd-indexed timestamps for onsets (`timestamps(1:2:end)`)
- **Opto TTL separation**: Use `setdiff(timestamps, ttl.shockTTL)` to isolate light-only trials

## Dependencies
Requires functions from `../BA_fear_cond_project/` (BAfc toolbox) and Open Ephys analysis tools.

## Gergo_figure_1_claude.m - Publication-Ready Figure Configuration

Located in `Gergo_claude/Gergo_figure_1_claude.m`, this is an optimized version of `Gergo_figure_1.m` for publication export.

### Figure Specifications

**Figure 1 (Main)**: 60×25 cm at 300 DPI
- 2×4 tile layout comparing BA→PFC vs BA→DMS projection neurons
- Panels: A-H (raster examples, heatmaps, population activity, responsive percentages)

**Figure 2 (Supplementary)**: 45×25 cm at 300 DPI
- 2×3 tile layout for optotagging validation
- Panels: A-F (rasters, light response heatmaps, z-scored activity)

### Font Sizing for Export
All fonts scaled for 60×25 cm figure that will be resized in Word:
- `g.fontSize2 = 22`: Axis labels (scales to ~6-7pt)
- `g.fontSize1 = 25`: Titles (scales to ~8pt)
- `g.xlinewidth = 3`: Stimulus marker line
- `g.axisLineWidth = 2`: Axis lines
- `g.markerSize = 8`: Raster dots

### Figure Layout Details

**Panel A/E**: Example rasters (BA→PFC / BA→DMS)
- Black dashed line at stimulus onset (time 0)
- Red lightning symbol (⚡) at x=-0.05, positioned slightly left of dashed line
- Limited to first 80 trials for both panels
- Marker size: 8 (scales appropriately when exported)
- Y-label: Two lines with "BA→PFC"/"BA→DMS" (bold, fontSize2+8) on top, "Trials" (normal, fontSize2) below

**Panel B/F**: US response heatmaps
- Black dashed vertical line at time 0
- Red lightning symbol (⚡) at x=-0.05, near top of panel (y = size-2)
- Colorbar on right side showing z-score scale
- Titles: "US response" (non-bold, normal weight)
- Neurons sorted by onset latency (earliest first)

**Panel C/G**: Population activity traces
- Firing rate (Hz) with SEM shading
- Red lightning symbol (⚡) at x=-0.05, y=60 (top of plot)
- Line width: 3 for visibility after scaling
- Titles: "Population activity" (non-bold, normal weight)

**Panel D/H**: Stacked bar charts (responsive percentages)
- **Stacking order**: Non-responsive (gray) at bottom, responsive (red) on top
- Legend direction reversed to match visual order
- Shows total neuron count: `(n=X)` positioned at (0.85, 0.55) in normalized coordinates
- Only displays percentage labels if they round to ≥1% (avoids showing "0%")
- Titles: "Responsiveness to US" (non-bold, normal weight)

### Optotagging (Figure 2)

**Figure 2 (Supplementary)**: 35×25 cm at 300 DPI
- 2×2 tile layout (z-score panels removed)
- Panels: A-D (rasters and light response heatmaps only)

**Panels A/C**: Optotagging rasters (BA→PFC / BA→DMS)
- Y-axis limited to actual number of trials plotted
- Black dashed line at light stimulus onset
- Light blue shaded area from 0-10ms indicating blue light pulse duration
- Y-label: Two lines with "BA→PFC"/"BA→DMS" (bold, fontSize2+8) on top, "Trials" (normal, fontSize2) below
- Titles: "Optotagging" (non-bold, normal weight)

**Panels B/D**: Light response heatmaps
- Black dashed vertical line at time 0
- Colorbar showing z-score scale
- Titles: "Light response" (non-bold, normal weight)

### Key Parameter Settings
```matlab
g.pre_time = 1;
g.post_time = 1;
g.bin_time = 0.001;
g.smoothvalue = 201;
g.test_time = 0.6;
g.exctestvalue = 2;
g.inhtestvalue = 2;
```

### Export Instructions
1. Run script to generate figures
2. Use `exportgraphics(gcf, 'figure_name.tiff', 'Resolution', 300)`
3. Insert into Word/PowerPoint
4. Resize to desired final dimensions (~18-20 cm width)
5. Fonts will scale appropriately to publication standards

### Design Decisions
- **Black dashed lines** at stimulus onset (time 0) for all plots - more visible than gray
- **Baseline z-scoring**: Z-score calculated from baseline period only (pre-stimulus) for US response
- **Color scale**: Uses percentile clipping (0.5-99.5%) to avoid extreme value dominance
- **Neuron sorting**: By onset latency (earliest responders first) for clearer temporal patterns
- **Consistent styling**: All panels use same line widths, fonts, and marker sizes for uniformity
- **Lightning symbols**: Red lightning symbols (⚡) rendered as figure-level annotations using `annotation('textbox',...)` to ensure they appear on top of all plot elements including dashed lines. Positioned at x=-0.05 (slightly left of stimulus line) using normalized figure coordinates calculated from axes positions.

### Sample Sizes
- BA→PFC neurons: n=15
- BA→DMS neurons: n=17

### Lightning Symbol Implementation
Lightning symbols are placed on panels A, B, C, E, F, G using figure-level annotations after all axes elements are drawn. This ensures they appear on top of everything (raster dots, heatmap data, traces, and dashed lines). Example positioning:
```matlab
annotation('textbox', [0.14, 0.925, 0.01, 0.01], ...
    'String', '\bf\fontsize{35}\color{red}⚡', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FitBoxToText', 'off');
```
Normalized figure coordinates are manually adjusted for each panel to position at shock onset.

## Publication Materials

### Figure_Legends_and_Methods.txt
Located in `Gergo_claude/Figure_Legends_and_Methods.txt`, contains:

**Complete figure legends** for both figures:
- Figure 1 panels A-H: Detailed descriptions of raster plots, heatmaps, population traces, and responsive neuron proportions
- Figure 2 panels A-D: Supplementary optotagging validation descriptions

**Methods section** following standard neuroscience format (Journal of Neuroscience style):
- Data Analysis section with complete technical details:
  - 1-ms binning for PSTH construction (-1000 to +1000 ms window)
  - Baseline-specific z-score normalization (mean and SD from -1000 to 0 ms)
  - Savitzky-Golay smoothing (201 ms frame for shock, 3 ms for opto)
  - Response classification criteria (|z| ≥ 2 within 0-600 ms post-stimulus)
  - Onset latency calculation (first bin exceeding threshold)
  - Neuron sorting by onset latency in heatmaps
  - Firing rate conversion and SEM calculation for population traces

**Updated with actual sample sizes**:
- BA→PFC: n=15 neurons
- BA→DMS: n=17 neurons

The methods text is formatted for direct insertion into manuscript with proper technical detail on all analysis parameters matching the actual script implementation.
