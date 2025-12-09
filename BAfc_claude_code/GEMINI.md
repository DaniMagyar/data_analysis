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
*   **Decimal Standard:** 2 decimal places for percentages, firing rates, and latencies. 3 decimal places for spike width. For p-values: use "p < 0.001" if the value is less than 0.001; otherwise, report the exact p-value with 3 decimal places (e.g., p = 0.041). **It is critical to ensure these p-values are precisely extracted from the statistical sources and accurately reflected in the text.**
4.  **DOCX Formatting:** When generating or modifying DOCX files, always use "Times New Roman" font and 12pt font size for body text.

## Tables

*   **Table 1:** Baseline characterization of firing rates and bursting properties
*   **Table 2:** Block-by-block firing rate evolution across regions and stimuli
*   **Table 3:** Comparison of US-evoked excitatory responses across amygdala regions
*   **Table 4:** Firing rate dynamics and signal gain in LA and AStria
*   **Table 5:** Comparison of signal gain (firing rate) between LA and AStria across stimulus conditions
*   **Table 6:** Monosynaptic firing rate responses (0-25ms) in LA and AStria
*   **Table 7:** Comparison of monosynaptic response magnitude and onset latency between LA and AStria
*   **Table 8:** Comparison of Principal Neuron (PN) and Interneuron (IN) responses in LA
*   **Table 9:** Modulation of Firing Rates by PV Interneuron Inhibition

## Manuscript Figures & Statistical Sources

This section details the visual data structure and where to find the corresponding statistical outputs.

### Data Source Map (Excel Files)
*   **Figure 1:** `figures/figure_1/figure_1_data.xlsx`
*   **Figure 2:** `figures/figure_2/figure_2_data.xlsx`
*   **Figure 3:** `figures/figure_3/figure_3_data.xlsx`
*   **Figure 4:** `figures/figure_4/figure_4_data.xlsx`
*   **Figure 5:** `figures/figure_5/figure_5_data.xlsx`
*   **Figure S2:** `figures/supplementary_figure_2/supplementary_figure_2_data.xlsx`
*   **Figure S3:** `figures/supplementary_figure_3/supplementary_figure_3_data.xlsx`
*   **Figure S4:** `figures/supplementary_figure_4/supplementary_figure_4_data.xlsx`

### Main Figures

**Figure 1: Characterization of neurons in the basal amygdala**
*   **A:** Experimental setup (Schematic).
*   **B:** Example traces (Raw data).
*   **C:** Example ISI and ACG (PN vs IN).
*   **D:** Spike feature classification (Waveform width vs Firing Rate).
*   **E:** Firing rate distributions (PN/IN/Region).
*   **F:** Burst index distributions.
*   **Source Stats:** `figures/figure_1/figure_1_data.xlsx` (PN vs IN waveform/FR stats, Burst Index).

**Figure 2: Regional differences in neuronal response profiles**
*   **A:** Example rasters (CS-sel, US-sel, Multi, Inhibited).
*   **B:** Population heatmaps (LA, BA, AStria, CeA) & Bar charts (Cluster proportions).
*   **C:** Chi-square significance matrix (Region vs Region comparisons).
*   **D:** US-evoked metrics (Delta FR magnitude, Response Length).
*   **Source Stats:** `figures/figure_2/figure_2_data.xlsx` (Chi-square tests, Cramér's V, Kruskal-Wallis for metrics).

**Figure 3: Similar integration in LA and AStria**
*   **A:** LA Heatmaps (CS, US, CS+US responses).
*   **B:** LA Lineplots (Firing rates per cluster).
*   **C:** LA Delta FR bar charts (CS vs US vs CS+US integration dynamics).
*   **D:** AStria Heatmaps.
*   **E:** AStria Lineplots.
*   **F:** AStria Delta FR bar charts.
*   **G:** Pie charts (Responsive proportions LA vs AStria).
*   **H:** Signal gain comparison (Delta FR across regions for all stimuli).
*   **Source Stats:** `figures/figure_3/figure_3_data.xlsx` (Wilcoxon signed-rank for integration, Wilcoxon rank-sum for gain).

**Figure 4: Monosynaptic responses**
*   **A:** LA Heatmaps (Monosynaptic window 0-25ms).
*   **B:** LA Delta FR bar charts (Monosynaptic magnitude).
*   **C:** AStria Heatmaps.
*   **D:** AStria Delta FR bar charts.
*   **E:** Pie charts (Monosynaptic proportions).
*   **F:** Region comparison (Delta FR magnitude).
*   **G:** Onset latency comparison.
*   **Source Stats:** `figures/figure_4/figure_4_data.xlsx` (Latencies, Magnitudes, Chi-square).

**Figure 5: Optogenetic modulation (PV silencing)**
*   **A:** Viral strategy schematic.
*   **B:** Histology (ArchT expression images).
*   **C:** Example rasters (Single neurons + Light).
*   **D:** Example lineplots (Single neurons + Light).
*   **E:** Pie charts (Proportion of enhanced neurons).
*   **F:** Population gain (Spaghetti plots showing FR change with light).
*   **Source Stats:** `figures/figure_5/figure_5_data.xlsx` (Enhancement proportions, gain stats).

### Supplementary Figures

**Figure S1: Electrode positions**
*   **A-E:** Maps, tracks, and schematics of recording sites.
*   **Source:** Histological reconstruction (No text stats).

**Figure S2: Stability of neuronal responses across trial blocks.**
*   **Source Stats:** `figures/supplementary_figure_2/supplementary_figure_2_data.xlsx`.

**Figure S3: PN vs IN: Cell type comparison of response profiles in the lateral amygdala.**
*   **Source Stats:** `figures/supplementary_figure_3/supplementary_figure_3_data.xlsx`.

**Figure S4: Light-inhibited PV interneurons during optogenetic manipulation.**
*   **Source Stats:** `figures/figure_5/figure_5_data.xlsx` (Manipulation check stats: baseline FR drops).

**Figure S5: Monosynaptic latency analysis (Related to Fig 4)**
*   **A:** Boxplots (Latency comparison across subtypes).
*   **B:** Permutation distribution.
*   **C:** Contingency table.
*   **D:** Chi-square stats.
*   **E:** Kruskal-Wallis for monosynaptic magnitude.
*   **Source Stats:** `figures/figure_4/figure_4_data.xlsx` (Latency section).