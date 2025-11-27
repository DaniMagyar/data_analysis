# Project Agent: LA-AStr Functional Unit Analysis

## Scope & Identity
I am a specialized local agent dedicated exclusively to the research project located in `C:\Users\dmagyar\Documents\data_analysis\BAfc_claude_code`. My primary mandate is to analyze, interpret, and write manuscript content for the "BAfc" (Basal Amygdala fear conditioning) dataset.

## Scientific Core Hypothesis
**The Lateral Amygdala (LA) and Amygdalostriatal Transition Area (AStria) operate as a unified, parallel processing module during fear conditioning.** 
*   **Novelty:** This challenges the canonical "Basolateral Amygdala" (BLA) model, which posits the LA and Basal Amygdala (BA) as the primary functional unit.
*   **Evidence:** 
    *   LA and AStria exhibit matched firing dynamics, response proportions, and multisensory integration.
    *   BA and Central Amygdala (CeA) are largely non-responsive to simple CS/US stimuli in this paradigm.
    *   Both LA and AStria receive parallel monosynaptic sensory inputs (matched latencies).
    *   Both regions share identical inhibitory control mechanisms (PV+ interneurons).

## Project Context
*   **Subjects:** Awake, head-fixed mice (N=24 total; 14 PV-Cre, 10 WT).
*   **Methods:** In vivo high-density electrophysiology (silicon probes) and optogenetics (PV-Cre/ArchT).
*   **Regions Recorded:** Lateral Amygdala (LA), Basal Amygdala (BA), Amygdalostriatal Transition Area (AStria), Central Amygdala (CeA).
*   **Stimuli:** Auditory Tone (CS) and Tail Shock (US).
*   **Analysis Environment:** MATLAB R2024b, Kilosort 2.5.

## Operational Guidelines
1.  **Local Focus:** I rely strictly on the files provided within this specific project directory (MATLAB scripts, text drafts, figure data).
2.  **Consistency:** All generated text (Results, Discussion, etc.) must align with the core hypothesis that LA and AStria are the functional unit.
*   Precise terminology: differentiating "AStria" from "Striatum".
*   **Decimal Standard:** 2 decimal places for percentages, firing rates, latencies, and p-values (except p < 0.001). 3 decimal places for spike width.

## Manuscript Figures & Statistical Sources

This section details the visual data structure and where to find the corresponding statistical outputs.

### Main Figures

**Figure 1: Characterization of neurons in the basal amygdala**
*   **A:** Experimental setup (Schematic).
*   **B:** Example traces (Raw data).
*   **C:** Example ISI and ACG (PN vs IN).
*   **D:** Spike feature classification (Waveform width vs Firing Rate).
*   **E:** Firing rate distributions (PN/IN/Region).
*   **F:** Burst index distributions.
*   **Source Stats:** `figures/figure_1/figure_1_stats.txt` (PN vs IN waveform/FR stats, Burst Index).

**Figure 2: Regional differences in neuronal response profiles**
*   **A:** Example rasters (CS-sel, US-sel, Multi, Inhibited).
*   **B:** Population heatmaps (LA, BA, AStria, CeA) & Bar charts (Cluster proportions).
*   **C:** Chi-square significance matrix (Region vs Region comparisons).
*   **D:** US-evoked metrics (Delta FR magnitude, Response Length).
*   **Source Stats:** `figures/figure_2/figure_2_stats.txt` (Chi-square tests, Cramér's V, Kruskal-Wallis for metrics).

**Figure 3: Similar integration in LA and AStria**
*   **A:** LA Heatmaps (CS, US, CS+US responses).
*   **B:** LA Lineplots (Firing rates per cluster).
*   **C:** LA Delta FR bar charts (CS vs US vs CS+US integration dynamics).
*   **D:** AStria Heatmaps.
*   **E:** AStria Lineplots.
*   **F:** AStria Delta FR bar charts.
*   **G:** Pie charts (Responsive proportions LA vs AStria).
*   **H:** Signal gain comparison (Delta FR across regions for all stimuli).
*   **Source Stats:** `figures/figure_3/figure_3_stats.txt` (Wilcoxon signed-rank for integration, Wilcoxon rank-sum for gain).

**Figure 4: Monosynaptic responses**
*   **A:** LA Heatmaps (Monosynaptic window 0-25ms).
*   **B:** LA Delta FR bar charts (Monosynaptic magnitude).
*   **C:** AStria Heatmaps.
*   **D:** AStria Delta FR bar charts.
*   **E:** Pie charts (Monosynaptic proportions).
*   **F:** Region comparison (Delta FR magnitude).
*   **G:** Onset latency comparison.
*   **Source Stats:** `figures/figure_4/figure_4_stats.txt` (Latencies, Magnitudes, Chi-square).

**Figure 5: Optogenetic modulation (PV silencing)**
*   **A:** Viral strategy schematic.
*   **B:** Histology (ArchT expression images).
*   **C:** Example rasters (Single neurons + Light).
*   **D:** Example lineplots (Single neurons + Light).
*   **E:** Pie charts (Proportion of enhanced neurons).
*   **F:** Population gain (Spaghetti plots showing FR change with light).
*   **Source Stats:** `figures/figure_5/figure_5_stats.txt` (Enhancement proportions, gain stats).

### Supplementary Figures

**Figure S1: Electrode positions**
*   **A-E:** Maps, tracks, and schematics of recording sites.
*   **Source:** Histological reconstruction (No text stats).

**Figure S2: Statistical analysis of cluster distributions (Related to Fig 2)**
*   **A:** Permutation test histogram.
*   **B:** Contingency table (Observed counts).
*   **C:** Cramér's V effect size matrix.
*   **D:** Kruskal-Wallis results for metrics.
*   **Source Stats:** `figures/figure_2/figure_2_stats.txt`.

**Figure S3: Statistical analysis of regional cluster distributions (Related to Fig 3)**
*   **A:** Permutation test (LA vs AStria).
*   **B:** Contingency table.
*   **C:** Stats summary.
*   **D:** Response metrics (Delta Spikes, Delta FR, Length).
*   **E:** Detailed Kruskal-Wallis tests.
*   **Source Stats:** `figures/figure_3/figure_3_stats.txt`.

**Figure S4: PN vs IN comparison in LA**
*   **A-E:** Heatmaps and metrics comparing Principal Neurons vs Interneurons.
*   **Source Stats:** Derived from `figures/figure_3/BAfc_figure_3_supp_PN_IN.m`.

**Figure S5: Monosynaptic latency analysis (Related to Fig 4)**
*   **A:** Boxplots (Latency comparison across subtypes).
*   **B:** Permutation distribution.
*   **C:** Contingency table.
*   **D:** Chi-square stats.
*   **E:** Kruskal-Wallis for monosynaptic magnitude.
*   **Source Stats:** `figures/figure_4/figure_4_stats.txt` (Latency section).

**Figure S6: Light-inhibited PV interneurons (Related to Fig 5)**
*   **A:** Example rasters (PV inhibition).
*   **B:** Heatmaps (LA PVs).
*   **C:** Lineplots (LA PV firing drops).
*   **D:** Heatmaps (AStria PVs).
*   **E:** Lineplots (AStria PV firing drops).
*   **Source Stats:** `figures/figure_5/figure_5_supp_stats.txt` (Manipulation check stats: baseline FR drops).