function waveform_classifier_main()
% WAVEFORM_CLASSIFIER_MAIN Main entry point for waveform classification
%
% This function asks the user whether to start a fresh session or load
% a saved curation session, then launches the curator UI

% Add path to functions
addpath(fileparts(mfilename('fullpath')));

% Ask user for mode
choice = questdlg('Do you want to start a fresh calculation or load a saved session?', ...
    'Waveform Classifier', ...
    'Fresh Calculation', 'Load Saved Session', 'Cancel', 'Fresh Calculation');

switch choice
    case 'Fresh Calculation'
        % Run fresh calculation
        run_fresh_calculation();

    case 'Load Saved Session'
        % Load saved session - user will select file
        load_saved_session();

    case 'Cancel'
        fprintf('Classification cancelled by user.\n');
        return;
end

end

function run_fresh_calculation()
% RUN_FRESH_CALCULATION Perform fresh waveform feature calculation

% Define all available recordings
all_recordings = {...
    'MD292_002_kilosort',...
    'MD293_001_kilosort',...
    'MD294_001_kilosort',...
    'MD295_001_kilosort',...
    'MD296_001_kilosort',...
    'MD297_001_kilosort',...
    'MD298_001_kilosort',...
    'MD299_001_kilosort',...
    'MD300_001_kilosort',...
    'MD304_001_kilosort',...
    'MD307_001_kilosort',...
    'MD309_001_kilosort',...
    'MD311_002_kilosort',...
    'MD312_001_kilosort',...
    'MD313_001_kilosort',...
    'MD314_001_kilosort',...
    'MD315_001_kilosort',...
    'MD316_002_kilosort',...
    'MD317_001_kilosort',...
    'MD318_001_kilosort',...
    'MD318_002_kilosort',...
    'MD319_003_kilosort'};

% Let user select recordings
[selected_idx, ok] = listdlg('ListString', all_recordings, ...
    'SelectionMode', 'multiple', ...
    'InitialValue', 1:length(all_recordings), ...
    'PromptString', 'Select recordings to analyze:', ...
    'ListSize', [300, 400]);

if ~ok
    fprintf('Recording selection cancelled by user.\n');
    return;
end

recordings = all_recordings(selected_idx);
fprintf('Selected %d recordings to analyze\n', length(recordings));

% Load neurons
g.cell_metrics = BAfc_load_neurons('recordings', recordings);

timeaxis = g.cell_metrics.waveforms.time{1};
waveforms_raw_orig = cell2mat(g.cell_metrics.waveforms.raw(:));
firing_rates = g.cell_metrics.firingRate;
brain_regions = g.cell_metrics.brainRegion;

% Extract metadata for each neuron
animals = g.cell_metrics.animal;
batch_ids = g.cell_metrics.batchIDs;
cell_ids = g.cell_metrics.cellID;
clu_ids = g.cell_metrics.cluID;

% Keep truly original waveforms for raw display (unchanged)
waveforms_raw_display = waveforms_raw_orig;

% Create flipped version for calculations
waveforms_raw_flipped = waveforms_raw_orig;
waveforms_raw_flipped(g.cell_metrics.polarity>0,:) = -waveforms_raw_flipped(g.cell_metrics.polarity>0,:);

% Use flipped waveforms for calculations (will be supersampled)
waveforms = waveforms_raw_flipped;

% Normalize for peak-to-peak amplitude
for i = 1:size(waveforms, 1)
    ptp_amp = max(waveforms(i, :)) - min(waveforms(i, :));
    waveforms(i, :) = waveforms(i, :) / ptp_amp;
end

% Parameters
use_supersampling = true;
supersample_factor = 10;
% Baseline is mean of first 25% of samples
baseline_samples = round(size(waveforms, 2) * 0.25);

% Calculate features
fprintf('Calculating waveform features...\n');
[trough_to_peak, half_width, waveforms_super, timeaxis_super] = ...
    calculate_waveform_features(waveforms, timeaxis, use_supersampling, supersample_factor, baseline_samples);

fprintf('Calculation complete. Launching curator UI...\n');

% Launch curator UI
manual_corrections = struct();
waveform_curator_ui(waveforms_super, timeaxis_super, trough_to_peak, half_width, baseline_samples, manual_corrections, waveforms_raw_display, timeaxis, firing_rates, timeaxis, brain_regions, animals, batch_ids, cell_ids, clu_ids);

end

function load_saved_session()
% LOAD_SAVED_SESSION Load a previously saved curation session

% Open file selection dialog
[filename, pathname] = uigetfile('*.mat', 'Select a saved session file');

% Check if user cancelled
if isequal(filename, 0)
    fprintf('Load cancelled by user.\n');
    return;
end

% Construct full file path
fullpath = fullfile(pathname, filename);
fprintf('Loading saved session from %s...\n', fullpath);

try
    loaded = load(fullpath);
catch ME
    warndlg(sprintf('Error loading file: %s', ME.message), 'Error');
    return;
end

% Check required fields
if ~isfield(loaded, 'waveforms_super') || ~isfield(loaded, 'timeaxis_super')
    warndlg('Saved session is missing required data. Starting fresh calculation.', 'Warning');
    run_fresh_calculation();
    return;
end

% Extract data
waveforms_super = loaded.waveforms_super;
timeaxis_super = loaded.timeaxis_super;
trough_to_peak = loaded.trough_to_peak;
half_width = loaded.half_width;
baseline_samples = loaded.baseline_samples;

if isfield(loaded, 'manual_corrections')
    manual_corrections = loaded.manual_corrections;
    fprintf('Loaded %d manual corrections\n', length(fieldnames(manual_corrections)));
else
    manual_corrections = struct();
end

% Load raw waveforms if available
if isfield(loaded, 'waveforms_raw') && isfield(loaded, 'timeaxis_raw')
    waveforms_raw = loaded.waveforms_raw;
    timeaxis_raw = loaded.timeaxis_raw;
else
    waveforms_raw = [];
    timeaxis_raw = [];
    fprintf('Warning: Raw waveforms not found in saved session.\n');
end

% Load firing rates if available
if isfield(loaded, 'firing_rates')
    firing_rates = loaded.firing_rates;
else
    firing_rates = [];
    fprintf('Warning: Firing rates not found in saved session.\n');
end

% Load timeaxis if available
if isfield(loaded, 'timeaxis')
    timeaxis = loaded.timeaxis;
else
    timeaxis = [];
    fprintf('Warning: Time axis not found in saved session.\n');
end

% Load brain regions if available
if isfield(loaded, 'brain_regions')
    brain_regions = loaded.brain_regions;
else
    brain_regions = {};
    fprintf('Warning: Brain regions not found in saved session.\n');
end

% Load metadata if available
if isfield(loaded, 'animals')
    animals = loaded.animals;
else
    animals = {};
    fprintf('Warning: Animals not found in saved session.\n');
end

if isfield(loaded, 'batch_ids')
    batch_ids = loaded.batch_ids;
else
    batch_ids = [];
    fprintf('Warning: Batch IDs not found in saved session.\n');
end

if isfield(loaded, 'cell_ids')
    cell_ids = loaded.cell_ids;
else
    cell_ids = {};
    fprintf('Warning: Cell IDs not found in saved session.\n');
end

if isfield(loaded, 'clu_ids')
    clu_ids = loaded.clu_ids;
else
    clu_ids = [];
    fprintf('Warning: Cluster IDs not found in saved session.\n');
end

fprintf('Session loaded. Launching curator UI...\n');

% Launch curator UI
waveform_curator_ui(waveforms_super, timeaxis_super, trough_to_peak, half_width, baseline_samples, manual_corrections, waveforms_raw, timeaxis_raw, firing_rates, timeaxis, brain_regions, animals, batch_ids, cell_ids, clu_ids);

end
