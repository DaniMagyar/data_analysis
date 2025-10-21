========================================
BAfc_figure_optogenetics.m - AUTOMATIC RESPONSIVENESS DETECTION
========================================

OVERVIEW:
The script now supports both MANUAL and AUTOMATIC neuron selection modes.

========================================
HOW TO USE:
========================================

At the top of the script (line 23), change this variable:

NEURON_SELECTION_METHOD = 'automatic';  % For automatic detection
   OR
NEURON_SELECTION_METHOD = 'manual';     % For manual selection


========================================
AUTOMATIC MODE PARAMETERS (lines 26-28):
========================================

g.responsiveness_zscore_threshold = 3;    % Peak z-score must exceed this value
g.responsiveness_prob_threshold = 0.15;   % At least 15% of trials must have spikes
g.responsiveness_window = 0.05;           % Time window for detection (50ms after stimulus)


========================================
RESPONSIVENESS CRITERIA:
========================================

A neuron is classified as RESPONSIVE if it meets BOTH criteria:

1. Peak Z-score ≥ threshold in response window
   - Calculated from z-scored PSTH (baseline = 5s pre-stimulus)
   - Response window = 0 to 50ms after stimulus onset

2. Response Probability ≥ threshold
   - Fraction of trials with at least 1 spike in response window
   - Default: at least 15% of trials must respond


========================================
NEURON CATEGORIES:
========================================

The algorithm detects 4 stimulus conditions:
1. CS no-light (triptest_sound_only)
2. CS light (triptest_sound_only_light)
3. US no-light (triptest_shocks_only)
4. US light (triptest_shocks_only_light)

Then classifies neurons as:

CS RESPONSIVE:
  - Responsive to CS during NO-LIGHT **OR** LIGHT condition
  - Logical: responsive_CS_nolight | responsive_CS_light

US RESPONSIVE:
  - Responsive to US during NO-LIGHT **OR** LIGHT condition
  - Logical: responsive_US_nolight | responsive_US_light

CSUS RESPONSIVE:
  - Responsive to BOTH CS and US (in either light condition)
  - Logical: CS_responsive & US_responsive

CS-ONLY:
  - Responsive to CS but NOT US
  - Logical: CS_responsive & ~US_responsive

US-ONLY:
  - Responsive to US but NOT CS
  - Logical: US_responsive & ~CS_responsive


========================================
EXAMPLE NEURONS:
========================================

Automatically selected from responsive populations:
- Example PN: First pyramidal neuron in US-responsive group
- Example IN: First interneuron in US-responsive group

Fallback if no responsive neurons found:
- Example PN: First PN in entire dataset
- Example IN: First IN in entire dataset


========================================
OUTPUT:
========================================

When running in automatic mode, you'll see:

===== AUTOMATIC RESPONSIVENESS DETECTION =====

Analyzing CS no-light...
  Found X responsive neurons

Analyzing CS light...
  Found X responsive neurons

Analyzing US no-light...
  Found X responsive neurons

Analyzing US light...
  Found X responsive neurons

===== RESPONSIVENESS SUMMARY =====
CS responsive (CS no-light OR light): X neurons
US responsive (US no-light OR light): X neurons
CSUS responsive (both CS and US): X neurons
CS-only responsive: X neurons
US-only responsive: X neurons

Example neurons:
  PN: neuron XXX (LA, PN)
  IN: neuron XXX (BA, IN)


========================================
MANUAL MODE:
========================================

If you set NEURON_SELECTION_METHOD = 'manual', the script uses:

idx_example_PN_manual = 423;
idx_example_IN_manual = 97;
idx_CS_responsive_manual = [133, 238, 94, 105, 117, 128, 147];
idx_US_responsive_manual = [57, 82, 74, 91, 96, 97, 127, 409, 418, 423, 446, 448, 94, 105, 117, 128, 147];

These are defined at lines 51-62 in the script.


========================================
TUNING RECOMMENDATIONS:
========================================

If you get TOO MANY responsive neurons:
  - Increase g.responsiveness_zscore_threshold (try 4 or 5)
  - Increase g.responsiveness_prob_threshold (try 0.20 or 0.25)
  - Decrease g.responsiveness_window (try 0.025 = 25ms)

If you get TOO FEW responsive neurons:
  - Decrease g.responsiveness_zscore_threshold (try 2 or 2.5)
  - Decrease g.responsiveness_prob_threshold (try 0.10)
  - Increase g.responsiveness_window (try 0.10 = 100ms)


========================================
BACKUP FILES:
========================================

The original script is saved as:
  BAfc_figure_optogenetics_backup.m
  BAfc_figure_optogenetics_old.m

You can restore the original by:
  copy BAfc_figure_optogenetics_backup.m BAfc_figure_optogenetics.m

========================================

========================================
BRAIN REGION FILTERING (NEW):
========================================

Line 26:
g.brain_regions_filter = {'LA', 'BA'};

The script now only analyzes neurons from specified brain regions.
Default: Only LA (lateral amygdala) and BA (basal amygdala) neurons.

To include other regions, modify this parameter:
g.brain_regions_filter = {'LA', 'BA', 'CeA'};  % Add central amygdala
g.brain_regions_filter = {'LA'};               % Only lateral amygdala
g.brain_regions_filter = {};                   % All regions (no filter)

The filter applies to BOTH manual and automatic modes.

Console output will show:
Total neurons: XXX
Neurons in LA/BA: XXX

========================================
