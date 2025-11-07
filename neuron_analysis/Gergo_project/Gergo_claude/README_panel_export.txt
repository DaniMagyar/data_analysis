INSTRUCTIONS FOR GENERATING PUBLICATION-READY FIGURE PANELS
===========================================================

Modified by Claude Code to:
1. Change all fonts to Arial with 10pt size
2. Remove panel labels (A, B, C, D, etc.)
3. Save each panel as a separate high-resolution PNG file

HOW TO USE:
-----------

Step 1: Run the main figure script
   >> Gergo_figure_1_claude

   This will:
   - Generate the combined figure (Figure 1 with 8 panels)
   - Generate the supplementary figure (Figure 2 with 4 panels)
   - Display them in MATLAB figure windows
   - All fonts will be Arial 10pt
   - No panel labels will be shown

Step 2: Save individual panels
   >> save_panels_separately

   This will:
   - Create 8 separate PNG files for the main figure panels:
     * Panel_1_PFC_example_raster.png
     * Panel_2_PFC_heatmap.png
     * Panel_3_PFC_firing_rate.png
     * Panel_4_PFC_responsiveness.png
     * Panel_5_DMS_example_raster.png
     * Panel_6_DMS_heatmap.png
     * Panel_7_DMS_firing_rate.png
     * Panel_8_DMS_responsiveness.png

   - Create 4 separate PNG files for the supplementary figure:
     * Panel_S1_PFC_opto_raster.png
     * Panel_S2_PFC_opto_heatmap.png
     * Panel_S3_DMS_opto_raster.png
     * Panel_S4_DMS_opto_heatmap.png

   All files are saved to:
   C:\Users\dmagyar\Documents\data_analysis\neuron_analysis\Gergo_project\Gergo_claude

PANEL DIMENSIONS:
-----------------
- Standard panels: 8 cm × 6 cm
- Heatmap panels: 10 cm × 6 cm (extra width for colorbar)
- Bar chart panels: 4.8 cm × 6 cm (narrower for single bar)
- Resolution: 300 DPI (publication quality)
- Font: Arial 10pt (final size after insertion into Word)

USING THE PANELS IN WORD/POWERPOINT:
-------------------------------------
1. Insert the PNG files into your document using Insert > Pictures
2. The panels are already sized correctly - just drag and drop
3. Fonts will be 10pt Arial as specified
4. No need to resize - the dimensions are publication-ready
5. You can add panel labels (A, B, C, etc.) in Word if needed

CUSTOMIZING PANEL SIZE:
-----------------------
If you need different dimensions, edit the variables in save_panels_separately.m:
   panel_width = 8;   % cm - change as needed
   panel_height = 6;  % cm - change as needed

Then re-run: save_panels_separately

NOTES:
------
- The script uses exportgraphics() which requires MATLAB R2020a or later
- If you encounter any errors, ensure all variables from Gergo_figure_1_claude
  are still in the workspace before running save_panels_separately
- The original combined figures are still generated and can be saved manually
  using File > Save As in the figure window if needed
