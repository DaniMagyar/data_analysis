function waveform_curator_ui(waveforms_super, timeaxis_super, trough_to_peak, half_width, baseline_samples, manual_corrections, waveforms_raw, timeaxis_raw, firing_rates, timeaxis, brain_regions, animals, batch_ids, cell_ids, clu_ids)
% WAVEFORM_CURATOR_UI Interactive UI for manual curation of waveform features
%
% Inputs:
%   waveforms_super - Supersampled waveforms matrix
%   timeaxis_super - Supersampled time axis
%   trough_to_peak - Calculated trough-to-peak values
%   half_width - Calculated half-width values
%   baseline_samples - Number of samples before trough for baseline
%   manual_corrections - Struct containing manual corrections (optional)
%   waveforms_raw - Raw waveforms matrix (optional)
%   timeaxis_raw - Raw time axis (optional)
%   firing_rates - Mean firing rate for each neuron (optional)
%   timeaxis - Original time axis in milliseconds (optional)
%   brain_regions - Cell array of brain region labels (optional)
%   animals - Cell array of animal IDs (optional)
%   batch_ids - Array of batch IDs (optional)
%   cell_ids - Cell array of cell IDs (optional)
%   clu_ids - Array of cluster IDs (optional)

if nargin < 6
    manual_corrections = struct();
end

if nargin < 8
    waveforms_raw = [];
    timeaxis_raw = [];
end

if nargin < 9
    firing_rates = [];
end

if nargin < 10
    timeaxis = [];
end

if nargin < 11
    brain_regions = {};
end

if nargin < 12
    animals = {};
end

if nargin < 13
    batch_ids = [];
end

if nargin < 14
    cell_ids = {};
end

if nargin < 15
    clu_ids = [];
end

num_neurons = size(waveforms_super, 1);

% Get unique brain regions for filter
if ~isempty(brain_regions)
    unique_regions = unique(brain_regions);
    unique_regions = ['All', unique_regions(:)'];  % Add 'All' option
else
    unique_regions = {'All'};
end

% Get unique animals for filter
if ~isempty(animals)
    if iscell(animals)
        unique_animals = unique(animals);
        unique_animals = ['All', unique_animals(:)'];  % Add 'All' option
    else
        unique_animals_nums = unique(animals);
        unique_animals = cell(1, length(unique_animals_nums) + 1);
        unique_animals{1} = 'All';
        for i = 1:length(unique_animals_nums)
            unique_animals{i+1} = num2str(unique_animals_nums(i));
        end
    end
else
    unique_animals = {'All'};
end

% Interactive visualization
fig_interactive = figure('Name', 'Interactive Waveform Curator', 'Position', [50 50 1400 700]);

% Create UI data structure
ui_data.current_neuron = 1;
ui_data.num_neurons = num_neurons;
ui_data.timeaxis_super = timeaxis_super;
ui_data.waveforms_super = waveforms_super;
ui_data.trough_to_peak = trough_to_peak;
ui_data.half_width = half_width;
ui_data.baseline_samples = baseline_samples;
ui_data.manual_corrections = manual_corrections;
ui_data.edit_mode = false;
ui_data.dragging_point = '';
ui_data.waveforms_raw = waveforms_raw;
ui_data.timeaxis_raw = timeaxis_raw;
ui_data.firing_rates = firing_rates;
ui_data.timeaxis = timeaxis;
ui_data.panning_scatter = false;
ui_data.pan_start_point = [];
ui_data.brain_regions = brain_regions;
ui_data.selected_region = 'All';
ui_data.animals = animals;
ui_data.selected_animal = 'All';
ui_data.batch_ids = batch_ids;
ui_data.cell_ids = cell_ids;
ui_data.clu_ids = clu_ids;
ui_data.filtered_indices = 1:num_neurons;

% Create UI controls
ui_data.prev_btn = uicontrol('Style', 'pushbutton', 'String', 'Previous', ...
    'Position', [20 20 80 30], 'Callback', @prevNeuron);
ui_data.next_btn = uicontrol('Style', 'pushbutton', 'String', 'Next', ...
    'Position', [110 20 80 30], 'Callback', @nextNeuron);
ui_data.neuron_text = uicontrol('Style', 'text', 'String', sprintf('Neuron: %d / %d', ui_data.current_neuron, ui_data.num_neurons), ...
    'Position', [200 20 150 30], 'FontSize', 10);
ui_data.edit_btn = uicontrol('Style', 'togglebutton', 'String', 'Edit Mode OFF', ...
    'Position', [360 20 100 30], 'Callback', @toggleEditMode);
ui_data.save_btn = uicontrol('Style', 'pushbutton', 'String', 'Save All', ...
    'Position', [470 20 80 30], 'Callback', @saveCorrections);
ui_data.reset_btn = uicontrol('Style', 'pushbutton', 'String', 'Reset Current', ...
    'Position', [560 20 100 30], 'Callback', @resetCurrent);
ui_data.invert_btn = uicontrol('Style', 'pushbutton', 'String', 'Invert Polarity', ...
    'Position', [670 20 100 30], 'Callback', @invertPolarity);
ui_data.goto_text = uicontrol('Style', 'edit', 'String', '1', ...
    'Position', [780 20 60 30], 'FontSize', 10);
ui_data.goto_btn = uicontrol('Style', 'pushbutton', 'String', 'Go To', ...
    'Position', [850 20 60 30], 'Callback', @gotoNeuron);
ui_data.add_halfwidth_btn = uicontrol('Style', 'pushbutton', 'String', 'Add Half-Width', ...
    'Position', [920 20 100 30], 'Callback', @addHalfWidth);
ui_data.region_label = uicontrol('Style', 'text', 'String', 'Region:', ...
    'Position', [1030 20 50 30], 'FontSize', 10);
ui_data.region_popup = uicontrol('Style', 'popupmenu', 'String', unique_regions, ...
    'Position', [1085 20 100 30], 'Callback', @changeRegion);
ui_data.animal_label = uicontrol('Style', 'text', 'String', 'Animal:', ...
    'Position', [1195 20 50 30], 'FontSize', 10);
ui_data.animal_popup = uicontrol('Style', 'popupmenu', 'String', unique_animals, ...
    'Position', [1250 20 100 30], 'Callback', @changeAnimal);

% Plot areas - filtered waveform top left, raw waveform bottom left, scatter on right
ui_data.ax = axes('Position', [0.05 0.55 0.5 0.4]);
ui_data.ax_raw = axes('Position', [0.05 0.15 0.5 0.35]);
ui_data.ax_scatter = axes('Position', [0.6 0.2 0.35 0.7]);

% Draggable point handles
ui_data.point_handles = struct('trough', [], 'peak', [], 'desc', [], 'asc', []);

% Scatter plot handle for current neuron marker
ui_data.current_marker = [];

% Store data in figure
guidata(fig_interactive, ui_data);

% Set up mouse callbacks on figure
set(fig_interactive, 'WindowButtonMotionFcn', @mouseMove);
set(fig_interactive, 'WindowButtonUpFcn', @mouseUp);
set(fig_interactive, 'WindowScrollWheelFcn', @mouseScroll);

updatePlot(fig_interactive);

    function prevNeuron(src, ~)
        ui_data = guidata(src);
        % Find previous neuron in filtered list
        current_idx = find(ui_data.filtered_indices == ui_data.current_neuron, 1);
        if ~isempty(current_idx) && current_idx > 1
            ui_data.current_neuron = ui_data.filtered_indices(current_idx - 1);
        elseif isempty(current_idx) && ~isempty(ui_data.filtered_indices)
            ui_data.current_neuron = ui_data.filtered_indices(1);
        end
        set(ui_data.goto_text, 'String', num2str(ui_data.current_neuron));
        guidata(src, ui_data);
        updatePlot(gcbf);
    end

    function nextNeuron(src, ~)
        ui_data = guidata(src);
        % Find next neuron in filtered list
        current_idx = find(ui_data.filtered_indices == ui_data.current_neuron, 1);
        if ~isempty(current_idx) && current_idx < length(ui_data.filtered_indices)
            ui_data.current_neuron = ui_data.filtered_indices(current_idx + 1);
        elseif isempty(current_idx) && ~isempty(ui_data.filtered_indices)
            ui_data.current_neuron = ui_data.filtered_indices(1);
        end
        set(ui_data.goto_text, 'String', num2str(ui_data.current_neuron));
        guidata(src, ui_data);
        updatePlot(gcbf);
    end

    function changeRegion(src, ~)
        ui_data = guidata(src);
        region_idx = get(src, 'Value');
        regions_list = get(src, 'String');
        ui_data.selected_region = regions_list{region_idx};

        % Save updated data before calling updateFilteredIndices
        guidata(src, ui_data);

        % Update filtered indices based on both region and animal
        updateFilteredIndices(src);
    end

    function changeAnimal(src, ~)
        ui_data = guidata(src);
        animal_idx = get(src, 'Value');
        animals_list = get(src, 'String');
        ui_data.selected_animal = animals_list{animal_idx};

        % Save updated data before calling updateFilteredIndices
        guidata(src, ui_data);

        % Update filtered indices based on both region and animal
        updateFilteredIndices(src);
    end

    function updateFilteredIndices(src)
        ui_data = guidata(src);

        % Start with all neurons
        mask = true(ui_data.num_neurons, 1);

        % Debug: Show animal data type and sample values
        if ~isempty(ui_data.animals)
            fprintf('Animals data type: %s\n', class(ui_data.animals));
            if iscell(ui_data.animals)
                fprintf('Sample animals: %s, %s, %s\n', ui_data.animals{1}, ui_data.animals{min(2,length(ui_data.animals))}, ui_data.animals{min(3,length(ui_data.animals))});
            else
                fprintf('Sample animals: %s, %s, %s\n', num2str(ui_data.animals(1)), num2str(ui_data.animals(min(2,length(ui_data.animals)))), num2str(ui_data.animals(min(3,length(ui_data.animals)))));
            end
        end
        fprintf('Selected animal from dropdown: "%s"\n', ui_data.selected_animal);

        % Apply region filter
        if ~strcmp(ui_data.selected_region, 'All') && ~isempty(ui_data.brain_regions)
            if iscell(ui_data.brain_regions)
                region_mask = strcmp(ui_data.brain_regions, ui_data.selected_region);
                mask = mask & region_mask(:);  % Ensure column vector
                fprintf('Region filter applied: %d neurons match "%s"\n', sum(region_mask), ui_data.selected_region);
            else
                region_mask = (ui_data.brain_regions == ui_data.selected_region);
                mask = mask & region_mask(:);  % Ensure column vector
                fprintf('Region filter applied: %d neurons match "%s"\n', sum(region_mask), ui_data.selected_region);
            end
        end

        % Apply animal filter
        if ~strcmp(ui_data.selected_animal, 'All') && ~isempty(ui_data.animals)
            if iscell(ui_data.animals)
                animal_mask = strcmp(ui_data.animals, ui_data.selected_animal);
                fprintf('Animal filter (cell): %d neurons match "%s"\n', sum(animal_mask), ui_data.selected_animal);
                mask = mask & animal_mask(:);  % Ensure column vector
            else
                % If animals is numeric but selected_animal is string (from dropdown)
                selected_animal_num = str2double(ui_data.selected_animal);
                if ~isnan(selected_animal_num)
                    animal_mask = (ui_data.animals == selected_animal_num);
                    fprintf('Animal filter (numeric): %d neurons match %d\n', sum(animal_mask), selected_animal_num);
                    mask = mask & animal_mask(:);  % Ensure column vector
                else
                    % Try string comparison anyway
                    animal_mask = strcmp(cellstr(num2str(ui_data.animals)), ui_data.selected_animal);
                    fprintf('Animal filter (string from numeric): %d neurons match "%s"\n', sum(animal_mask), ui_data.selected_animal);
                    mask = mask & animal_mask(:);  % Ensure column vector
                end
            end
        end

        ui_data.filtered_indices = find(mask);

        % Debug output
        fprintf('Filter changed - Region: %s, Animal: %s\n', ui_data.selected_region, ui_data.selected_animal);
        fprintf('Filtered neurons: %d out of %d\n', length(ui_data.filtered_indices), ui_data.num_neurons);

        % Jump to first neuron in filtered list
        if ~isempty(ui_data.filtered_indices)
            ui_data.current_neuron = ui_data.filtered_indices(1);
            set(ui_data.goto_text, 'String', num2str(ui_data.current_neuron));
        end

        guidata(src, ui_data);
        updatePlot(gcbf);
    end

    function toggleEditMode(src, ~)
        ui_data = guidata(src);
        ui_data.edit_mode = get(src, 'Value');
        if ui_data.edit_mode
            set(src, 'String', 'Edit Mode ON');
        else
            set(src, 'String', 'Edit Mode OFF');
        end
        guidata(src, ui_data);
    end

    function saveCorrections(~, ~)
        ui_data = guidata(gcbf);
        manual_corrections = ui_data.manual_corrections;
        trough_to_peak = ui_data.trough_to_peak;
        half_width = ui_data.half_width;
        num_neurons = ui_data.num_neurons;
        waveforms_super = ui_data.waveforms_super;
        timeaxis_super = ui_data.timeaxis_super;
        baseline_samples = ui_data.baseline_samples;
        waveforms_raw = ui_data.waveforms_raw;
        timeaxis_raw = ui_data.timeaxis_raw;
        firing_rates = ui_data.firing_rates;
        timeaxis = ui_data.timeaxis;
        brain_regions = ui_data.brain_regions;
        animals = ui_data.animals;
        batch_ids = ui_data.batch_ids;
        cell_ids = ui_data.cell_ids;
        clu_ids = ui_data.clu_ids;

        % Calculate sample duration for conversion to milliseconds
        if ~isempty(ui_data.timeaxis)
            sample_duration = ui_data.timeaxis_super(2) - ui_data.timeaxis_super(1);
        else
            sample_duration = 1; % Keep as samples if no time axis
        end

        % Create corrected values incorporating manual corrections
        trough_to_peak_corrected_ms = zeros(num_neurons, 1);
        half_width_corrected_ms = zeros(num_neurons, 1);

        for i = 1:num_neurons
            neuron_key = sprintf('n%d', i);
            if isfield(manual_corrections, neuron_key)
                % Use manual corrections
                corr = manual_corrections.(neuron_key);
                trough_to_peak_corrected_ms(i) = (corr.peak_idx - corr.trough_idx) * sample_duration;
                half_width_corrected_ms(i) = (corr.idx_ascending - corr.idx_descending) * sample_duration;
            else
                % Use automatic detection
                trough_to_peak_corrected_ms(i) = trough_to_peak(i) * sample_duration;
                half_width_corrected_ms(i) = half_width(i) * sample_duration;
            end
        end

        % Ask user whether to save new file or append to existing
        save_choice = questdlg('Save options:', 'Save Data', ...
            'Save New File', 'Append to Existing File', 'Cancel', 'Save New File');

        if strcmp(save_choice, 'Cancel') || isempty(save_choice)
            fprintf('Save cancelled by user.\n');
            return;
        end

        if strcmp(save_choice, 'Save New File')
            % Open file save dialog
            [filename, pathname] = uiputfile('*.mat', 'Save Results As', 'sorting_results.mat');

            % Check if user cancelled
            if isequal(filename, 0)
                fprintf('Save cancelled by user.\n');
                return;
            end

            % Construct full file path
            fullpath = fullfile(pathname, filename);

            % Save the data
            save(fullpath, 'manual_corrections', 'trough_to_peak', 'half_width', ...
                'trough_to_peak_corrected_ms', 'half_width_corrected_ms', ...
                'num_neurons', 'waveforms_super', 'timeaxis_super', 'baseline_samples', ...
                'waveforms_raw', 'timeaxis_raw', 'firing_rates', 'timeaxis', 'brain_regions', ...
                'animals', 'batch_ids', 'cell_ids', 'clu_ids');
            fprintf('Saved %d neurons to: %s\n', num_neurons, fullpath);

        elseif strcmp(save_choice, 'Append to Existing File')
            % Open file selection dialog
            [filename, pathname] = uigetfile('*.mat', 'Select Existing File to Append To');

            % Check if user cancelled
            if isequal(filename, 0)
                fprintf('Append cancelled by user.\n');
                return;
            end

            % Construct full file path
            fullpath = fullfile(pathname, filename);

            % Load existing file
            fprintf('Loading existing file: %s\n', fullpath);
            existing = load(fullpath);

            % Check if file has required fields
            if ~isfield(existing, 'num_neurons') || ~isfield(existing, 'animals')
                errordlg('Selected file does not have the expected format!', 'Invalid File');
                return;
            end

            % Get starting index for new neurons
            start_idx = existing.num_neurons + 1;
            end_idx = existing.num_neurons + num_neurons;

            fprintf('Appending %d new neurons (indices %d to %d) to existing %d neurons\n', ...
                num_neurons, start_idx, end_idx, existing.num_neurons);

            % Append data arrays
            combined_trough_to_peak = [existing.trough_to_peak; trough_to_peak];
            combined_half_width = [existing.half_width; half_width];
            combined_trough_to_peak_corrected_ms = [existing.trough_to_peak_corrected_ms; trough_to_peak_corrected_ms];
            combined_half_width_corrected_ms = [existing.half_width_corrected_ms; half_width_corrected_ms];
            combined_waveforms_super = [existing.waveforms_super; waveforms_super];
            combined_waveforms_raw = [existing.waveforms_raw; waveforms_raw];
            combined_firing_rates = [existing.firing_rates; firing_rates];
            combined_brain_regions = [existing.brain_regions; brain_regions];
            combined_animals = [existing.animals; animals];
            combined_batch_ids = [existing.batch_ids; batch_ids];
            combined_cell_ids = [existing.cell_ids; cell_ids];
            combined_clu_ids = [existing.clu_ids; clu_ids];

            % Merge manual corrections with updated neuron indices
            combined_manual_corrections = existing.manual_corrections;
            for i = 1:num_neurons
                old_key = sprintf('n%d', i);
                new_key = sprintf('n%d', start_idx + i - 1);
                if isfield(manual_corrections, old_key)
                    combined_manual_corrections.(new_key) = manual_corrections.(old_key);
                end
            end

            % Update total neuron count
            combined_num_neurons = existing.num_neurons + num_neurons;

            % Save combined data
            save(fullpath, ...
                'combined_manual_corrections', 'combined_trough_to_peak', 'combined_half_width', ...
                'combined_trough_to_peak_corrected_ms', 'combined_half_width_corrected_ms', ...
                'combined_num_neurons', 'combined_waveforms_super', 'combined_waveforms_raw', ...
                'combined_firing_rates', 'combined_brain_regions', 'combined_animals', ...
                'combined_batch_ids', 'combined_cell_ids', 'combined_clu_ids', ...
                'timeaxis_super', 'timeaxis_raw', 'baseline_samples', 'timeaxis');

            % Also save with standard names (overwrite old variable names)
            manual_corrections = combined_manual_corrections;
            trough_to_peak = combined_trough_to_peak;
            half_width = combined_half_width;
            trough_to_peak_corrected_ms = combined_trough_to_peak_corrected_ms;
            half_width_corrected_ms = combined_half_width_corrected_ms;
            num_neurons = combined_num_neurons;
            waveforms_super = combined_waveforms_super;
            waveforms_raw = combined_waveforms_raw;
            firing_rates = combined_firing_rates;
            brain_regions = combined_brain_regions;
            animals = combined_animals;
            batch_ids = combined_batch_ids;
            cell_ids = combined_cell_ids;
            clu_ids = combined_clu_ids;

            save(fullpath, 'manual_corrections', 'trough_to_peak', 'half_width', ...
                'trough_to_peak_corrected_ms', 'half_width_corrected_ms', ...
                'num_neurons', 'waveforms_super', 'timeaxis_super', 'baseline_samples', ...
                'waveforms_raw', 'timeaxis_raw', 'firing_rates', 'timeaxis', 'brain_regions', ...
                'animals', 'batch_ids', 'cell_ids', 'clu_ids');

            fprintf('Successfully appended to: %s\n', fullpath);
            fprintf('Total neurons in file: %d (was %d, added %d)\n', num_neurons, existing.num_neurons, num_neurons - existing.num_neurons);
        end

        fprintf('  - trough_to_peak: original values (samples)\n');
        fprintf('  - half_width: original values (samples)\n');
        fprintf('  - trough_to_peak_corrected_ms: corrected values (milliseconds)\n');
        fprintf('  - half_width_corrected_ms: corrected values (milliseconds)\n');
        fprintf('  - Metadata: animals, batch_ids, cell_ids, clu_ids, brain_regions, firing_rates\n');
    end

    function resetCurrent(~, ~)
        ui_data = guidata(gcbf);
        neuron_key = sprintf('n%d', ui_data.current_neuron);
        if isfield(ui_data.manual_corrections, neuron_key)
            ui_data.manual_corrections = rmfield(ui_data.manual_corrections, neuron_key);
            guidata(gcbf, ui_data);
            updatePlot(gcbf);
            fprintf('Reset neuron %d to automatic detection\n', ui_data.current_neuron);
        end
    end

    function invertPolarity(~, ~)
        ui_data = guidata(gcbf);

        % Invert the supersampled waveform only (used for feature detection)
        ui_data.waveforms_super(ui_data.current_neuron, :) = -ui_data.waveforms_super(ui_data.current_neuron, :);

        % DO NOT invert waveforms_raw - it shows original data from g.cell_metrics

        % Clear any manual corrections for this neuron (will recalculate)
        neuron_key = sprintf('n%d', ui_data.current_neuron);
        if isfield(ui_data.manual_corrections, neuron_key)
            ui_data.manual_corrections = rmfield(ui_data.manual_corrections, neuron_key);
        end

        guidata(gcbf, ui_data);
        updatePlot(gcbf);
        fprintf('Inverted polarity for neuron %d (supersampled only)\n', ui_data.current_neuron);
    end

    function gotoNeuron(src, ~)
        ui_data = guidata(src);

        % Get the neuron number from the text field
        neuron_num = str2double(get(ui_data.goto_text, 'String'));

        % Check if valid
        if isnan(neuron_num) || neuron_num < 1 || neuron_num > ui_data.num_neurons
            warndlg(sprintf('Please enter a valid neuron number (1-%d)', ui_data.num_neurons), 'Invalid Input');
            return;
        end

        % Jump to that neuron
        ui_data.current_neuron = round(neuron_num);
        guidata(src, ui_data);
        updatePlot(gcbf);
    end

    function addHalfWidth(~, ~)
        ui_data = guidata(gcbf);
        neuron_key = sprintf('n%d', ui_data.current_neuron);
        waveform = ui_data.waveforms_super(ui_data.current_neuron, :);

        % Get or initialize correction structure
        if isfield(ui_data.manual_corrections, neuron_key)
            corr = ui_data.manual_corrections.(neuron_key);
            trough_idx = corr.trough_idx;
        else
            % Find trough
            [trough_val, trough_idx] = min(waveform);
            corr.trough_idx = trough_idx;

            % Find peak
            [~, peak_idx] = max(waveform(trough_idx:end));
            corr.peak_idx = peak_idx + trough_idx - 1;
        end

        % Calculate a reasonable half amplitude
        baseline_mean = mean(waveform(1:ui_data.baseline_samples));
        trough_val = waveform(trough_idx);
        half_amplitude = (trough_val + baseline_mean) / 2;

        % Find half-width points
        idx_descending = find(waveform(1:trough_idx) <= half_amplitude, 1, 'first');
        idx_ascending = find(waveform(trough_idx:end) >= half_amplitude, 1, 'first');

        if isempty(idx_descending)
            % If not found, place at 1/4 way to trough
            idx_descending = round(trough_idx * 0.75);
            half_amplitude = waveform(idx_descending);
        end

        if isempty(idx_ascending)
            % If not found, place at 1/4 way after trough
            idx_ascending = round(trough_idx + (length(waveform) - trough_idx) * 0.25);
            half_amplitude = waveform(idx_ascending);
        else
            idx_ascending = idx_ascending + trough_idx - 1;
        end

        % Store the values
        corr.idx_descending = idx_descending;
        corr.idx_ascending = idx_ascending;
        corr.half_amplitude = half_amplitude;

        ui_data.manual_corrections.(neuron_key) = corr;
        guidata(gcbf, ui_data);
        updatePlot(gcbf);

        fprintf('Added half-width markers for neuron %d (desc: %d, asc: %d, amplitude: %.3f)\n', ...
            ui_data.current_neuron, idx_descending, idx_ascending, half_amplitude);
    end

    function mouseDown(src, ~)
        ui_data = guidata(gcbf);

        % Check if clicked on a waveform editing point
        if ui_data.edit_mode
            point_type = get(src, 'UserData');
            if ~isempty(point_type)
                ui_data.dragging_point = point_type;
                guidata(gcbf, ui_data);
                return;
            end
        end
    end

    function scatterMouseDown(~, ~)
        ui_data = guidata(gcbf);

        % Get selection type (normal click vs right-click)
        selection_type = get(gcbf, 'SelectionType');

        % Normal left-click for panning
        if strcmp(selection_type, 'normal')
            pt = get(ui_data.ax_scatter, 'CurrentPoint');
            ui_data.panning_scatter = true;
            ui_data.pan_start_point = pt(1,1:2);
            ui_data.pan_start_xlim = xlim(ui_data.ax_scatter);
            ui_data.pan_start_ylim = ylim(ui_data.ax_scatter);
            guidata(gcbf, ui_data);
        end
    end

    function mouseMove(~, ~)
        ui_data = guidata(gcbf);

        % Check if panning scatter plot
        if ui_data.panning_scatter
            try
                pt = get(ui_data.ax_scatter, 'CurrentPoint');
                current_point = pt(1,1:2);

                % Calculate offset
                dx = current_point(1) - ui_data.pan_start_point(1);
                dy = current_point(2) - ui_data.pan_start_point(2);

                % Apply offset to original limits (directly without redraw)
                new_xlim = ui_data.pan_start_xlim - dx;
                new_ylim = ui_data.pan_start_ylim - dy;

                set(ui_data.ax_scatter, 'XLim', new_xlim, 'YLim', new_ylim);
                drawnow limitrate; % Limit redraw rate
            catch
                % Ignore errors during panning
            end
            return;
        end

        if isempty(ui_data.dragging_point)
            return;
        end

        pt = get(ui_data.ax, 'CurrentPoint');
        x_new_data = pt(1,1);

        % Convert from time to index if using timeaxis
        if ~isempty(ui_data.timeaxis)
            [~, x_new] = min(abs(ui_data.timeaxis_super - x_new_data));
        else
            x_new = round(x_new_data);
        end

        % Get current neuron data
        neuron_key = sprintf('n%d', ui_data.current_neuron);
        if isfield(ui_data.manual_corrections, neuron_key)
            corr = ui_data.manual_corrections.(neuron_key);
        else
            % Initialize with automatic values
            [trough_val, trough_idx] = min(ui_data.waveforms_super(ui_data.current_neuron, :));
            [~, peak_idx] = max(ui_data.waveforms_super(ui_data.current_neuron, trough_idx:end));
            peak_idx = peak_idx + trough_idx - 1;

            baseline_mean = mean(ui_data.waveforms_super(ui_data.current_neuron, 1:ui_data.baseline_samples));
            half_amplitude = (trough_val + baseline_mean) / 2;
            idx_descending = find(ui_data.waveforms_super(ui_data.current_neuron, 1:trough_idx) <= half_amplitude, 1, 'first');
            idx_ascending = find(ui_data.waveforms_super(ui_data.current_neuron, trough_idx:end) >= half_amplitude, 1, 'first') + trough_idx - 1;

            corr.trough_idx = trough_idx;
            corr.peak_idx = peak_idx;
            corr.idx_descending = idx_descending;
            corr.idx_ascending = idx_ascending;
            corr.half_amplitude = half_amplitude;
        end

        waveform = ui_data.waveforms_super(ui_data.current_neuron, :);
        x_new = max(1, min(x_new, length(waveform)));

        % Update the dragged point
        switch ui_data.dragging_point
            case 'trough'
                corr.trough_idx = x_new;
            case 'peak'
                corr.peak_idx = x_new;
            case 'desc'
                % Update descending index and get amplitude at that point
                corr.idx_descending = x_new;
                corr.half_amplitude = waveform(x_new);

                % Find corresponding ascending point at same amplitude
                trough_idx = corr.trough_idx;
                idx_ascending = find(waveform(trough_idx:end) >= corr.half_amplitude, 1, 'first');
                if ~isempty(idx_ascending)
                    corr.idx_ascending = idx_ascending + trough_idx - 1;
                end

            case 'asc'
                % Update ascending index and get amplitude at that point
                corr.idx_ascending = x_new;
                corr.half_amplitude = waveform(x_new);

                % Find corresponding descending point at same amplitude
                trough_idx = corr.trough_idx;
                idx_descending = find(waveform(1:trough_idx) <= corr.half_amplitude, 1, 'first');
                if ~isempty(idx_descending)
                    corr.idx_descending = idx_descending;
                end
        end

        ui_data.manual_corrections.(neuron_key) = corr;
        guidata(gcbf, ui_data);
        updatePlot(gcbf);
    end

    function mouseUp(~, ~)
        ui_data = guidata(gcbf);
        ui_data.dragging_point = '';
        ui_data.panning_scatter = false;
        guidata(gcbf, ui_data);
    end

    function mouseScroll(~, event)
        ui_data = guidata(gcbf);

        % Get current mouse position
        figPos = get(gcbf, 'CurrentPoint');

        % Get scatter plot position in figure coordinates
        scatterPos = get(ui_data.ax_scatter, 'Position');
        figSize = get(gcbf, 'Position');

        % Convert scatter plot position to pixels
        scatterLeft = scatterPos(1) * figSize(3);
        scatterBottom = scatterPos(2) * figSize(4);
        scatterWidth = scatterPos(3) * figSize(3);
        scatterHeight = scatterPos(4) * figSize(4);

        % Check if mouse is over scatter plot
        if figPos(1) >= scatterLeft && figPos(1) <= (scatterLeft + scatterWidth) && ...
           figPos(2) >= scatterBottom && figPos(2) <= (scatterBottom + scatterHeight)

            % Get current axis limits
            xl = xlim(ui_data.ax_scatter);
            yl = ylim(ui_data.ax_scatter);

            % Get current point in data coordinates
            pt = get(ui_data.ax_scatter, 'CurrentPoint');
            x_center = pt(1,1);
            y_center = pt(1,2);

            % Zoom factor
            zoom_factor = 0.9; % Zoom in/out by 10%

            if event.VerticalScrollCount < 0
                % Scroll up - zoom in
                scale = zoom_factor;
            else
                % Scroll down - zoom out
                scale = 1 / zoom_factor;
            end

            % Calculate new limits centered on mouse position
            x_range = (xl(2) - xl(1)) * scale;
            y_range = (yl(2) - yl(1)) * scale;

            % Keep center at mouse position
            x_ratio = (x_center - xl(1)) / (xl(2) - xl(1));
            y_ratio = (y_center - yl(1)) / (yl(2) - yl(1));

            new_xl = [x_center - x_range * x_ratio, x_center + x_range * (1 - x_ratio)];
            new_yl = [y_center - y_range * y_ratio, y_center + y_range * (1 - y_ratio)];

            % Apply new limits
            xlim(ui_data.ax_scatter, new_xl);
            ylim(ui_data.ax_scatter, new_yl);
        end
    end

    function updatePlot(fig)
        ui_data = guidata(fig);
        cla(ui_data.ax);
        hold(ui_data.ax, 'on');

        % Plot waveform with time axis
        if ~isempty(ui_data.timeaxis)
            plot(ui_data.ax, ui_data.timeaxis_super, ui_data.waveforms_super(ui_data.current_neuron, :), 'k-', 'LineWidth', 1.5);
        else
            plot(ui_data.ax, ui_data.waveforms_super(ui_data.current_neuron, :), 'k-', 'LineWidth', 1.5);
        end

        % Check for manual corrections
        neuron_key = sprintf('n%d', ui_data.current_neuron);
        if isfield(ui_data.manual_corrections, neuron_key)
            corr = ui_data.manual_corrections.(neuron_key);
            trough_idx = corr.trough_idx;
            peak_idx = corr.peak_idx;
            idx_descending = corr.idx_descending;
            idx_ascending = corr.idx_ascending;
            if isfield(corr, 'half_amplitude')
                half_amplitude = corr.half_amplitude;
            else
                % Calculate half_amplitude if not stored
                waveform = ui_data.waveforms_super(ui_data.current_neuron, :);
                baseline_mean = mean(waveform(1:ui_data.baseline_samples));
                half_amplitude = (waveform(trough_idx) + baseline_mean) / 2;
            end
        else
            % Use automatic detection
            [trough_val, trough_idx] = min(ui_data.waveforms_super(ui_data.current_neuron, :));
            [~, peak_idx] = max(ui_data.waveforms_super(ui_data.current_neuron, trough_idx:end));
            peak_idx = peak_idx + trough_idx - 1;

            baseline_mean = mean(ui_data.waveforms_super(ui_data.current_neuron, 1:ui_data.baseline_samples));
            half_amplitude = (trough_val + baseline_mean) / 2;
            idx_descending = find(ui_data.waveforms_super(ui_data.current_neuron, 1:trough_idx) <= half_amplitude, 1, 'first');
            idx_ascending = find(ui_data.waveforms_super(ui_data.current_neuron, trough_idx:end) >= half_amplitude, 1, 'first') + trough_idx - 1;
        end

        waveform = ui_data.waveforms_super(ui_data.current_neuron, :);
        trough_val = waveform(trough_idx);
        peak_val = waveform(peak_idx);

        % Convert indices to time if timeaxis available
        if ~isempty(ui_data.timeaxis)
            trough_time = ui_data.timeaxis_super(trough_idx);
            peak_time = ui_data.timeaxis_super(peak_idx);

            % Mark trough (draggable)
            h_trough = plot(ui_data.ax, trough_time, trough_val, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);
            set(h_trough, 'ButtonDownFcn', @mouseDown, 'UserData', 'trough');

            % Mark peak (draggable)
            h_peak = plot(ui_data.ax, peak_time, peak_val, 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b', 'LineWidth', 2);
            set(h_peak, 'ButtonDownFcn', @mouseDown, 'UserData', 'peak');

            % Draw trough-to-peak line
            plot(ui_data.ax, [trough_time peak_time], [trough_val peak_val], 'b--', 'LineWidth', 1);

            % Draw half-width (horizontal line at same amplitude)
            if ~isempty(idx_descending) && ~isempty(idx_ascending)
                desc_time = ui_data.timeaxis_super(idx_descending);
                asc_time = ui_data.timeaxis_super(idx_ascending);

                % Draw horizontal line at half amplitude
                plot(ui_data.ax, [desc_time asc_time], [half_amplitude half_amplitude], 'g-', 'LineWidth', 2);

                % Draggable half-width points at the half amplitude level
                h_desc = plot(ui_data.ax, desc_time, half_amplitude, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 2);
                set(h_desc, 'ButtonDownFcn', @mouseDown, 'UserData', 'desc');

                h_asc = plot(ui_data.ax, asc_time, half_amplitude, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 2);
                set(h_asc, 'ButtonDownFcn', @mouseDown, 'UserData', 'asc');
            end

            xlabel(ui_data.ax, 'Time (ms)');
        else
            % Mark trough (draggable)
            h_trough = plot(ui_data.ax, trough_idx, trough_val, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);
            set(h_trough, 'ButtonDownFcn', @mouseDown, 'UserData', 'trough');

            % Mark peak (draggable)
            h_peak = plot(ui_data.ax, peak_idx, peak_val, 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b', 'LineWidth', 2);
            set(h_peak, 'ButtonDownFcn', @mouseDown, 'UserData', 'peak');

            % Draw trough-to-peak line
            plot(ui_data.ax, [trough_idx peak_idx], [trough_val peak_val], 'b--', 'LineWidth', 1);

            % Draw half-width (horizontal line at same amplitude)
            if ~isempty(idx_descending) && ~isempty(idx_ascending)
                % Draw horizontal line at half amplitude
                plot(ui_data.ax, [idx_descending idx_ascending], [half_amplitude half_amplitude], 'g-', 'LineWidth', 2);

                % Draggable half-width points at the half amplitude level
                h_desc = plot(ui_data.ax, idx_descending, half_amplitude, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 2);
                set(h_desc, 'ButtonDownFcn', @mouseDown, 'UserData', 'desc');

                h_asc = plot(ui_data.ax, idx_ascending, half_amplitude, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 2);
                set(h_asc, 'ButtonDownFcn', @mouseDown, 'UserData', 'asc');
            end

            xlabel(ui_data.ax, 'Sample Index');
        end

        ylabel(ui_data.ax, 'Normalized Amplitude');

        trough_to_peak_val = peak_idx - trough_idx;

        % Check if half-width values are valid
        if ~isempty(idx_descending) && ~isempty(idx_ascending)
            half_width_val = idx_ascending - idx_descending;
        else
            half_width_val = NaN;
        end

        manual_marker = '';
        if isfield(ui_data.manual_corrections, neuron_key)
            manual_marker = ' [MANUAL]';
        end

        % Get metadata for current neuron
        metadata_str = '';

        % Firing rate
        if ~isempty(ui_data.firing_rates) && ui_data.current_neuron <= length(ui_data.firing_rates)
            metadata_str = sprintf('%s | FR: %.2f Hz', metadata_str, ui_data.firing_rates(ui_data.current_neuron));
        end

        % Brain region
        if ~isempty(ui_data.brain_regions) && ui_data.current_neuron <= length(ui_data.brain_regions)
            if iscell(ui_data.brain_regions)
                metadata_str = sprintf('%s | %s', metadata_str, ui_data.brain_regions{ui_data.current_neuron});
            else
                metadata_str = sprintf('%s | %s', metadata_str, ui_data.brain_regions(ui_data.current_neuron));
            end
        end

        % Animal
        if ~isempty(ui_data.animals) && ui_data.current_neuron <= length(ui_data.animals)
            if iscell(ui_data.animals)
                metadata_str = sprintf('%s | %s', metadata_str, ui_data.animals{ui_data.current_neuron});
            else
                metadata_str = sprintf('%s | %s', metadata_str, ui_data.animals(ui_data.current_neuron));
            end
        end

        % Batch ID
        if ~isempty(ui_data.batch_ids) && ui_data.current_neuron <= length(ui_data.batch_ids)
            metadata_str = sprintf('%s | Batch:%d', metadata_str, ui_data.batch_ids(ui_data.current_neuron));
        end

        % Cell ID
        if ~isempty(ui_data.cell_ids) && ui_data.current_neuron <= length(ui_data.cell_ids)
            if iscell(ui_data.cell_ids)
                metadata_str = sprintf('%s | Cell:%s', metadata_str, ui_data.cell_ids{ui_data.current_neuron});
            else
                metadata_str = sprintf('%s | Cell:%d', metadata_str, ui_data.cell_ids(ui_data.current_neuron));
            end
        end

        % Cluster ID
        if ~isempty(ui_data.clu_ids) && ui_data.current_neuron <= length(ui_data.clu_ids)
            metadata_str = sprintf('%s | Clu:%d', metadata_str, ui_data.clu_ids(ui_data.current_neuron));
        end

        % Format half-width for title
        if isnan(half_width_val)
            hw_str = 'N/A';
        else
            hw_str = sprintf('%d', half_width_val);
        end

        title(ui_data.ax, sprintf('Neuron %d/%d | T2P: %d samples | HW: %s samples%s%s', ...
            ui_data.current_neuron, ui_data.num_neurons, trough_to_peak_val, hw_str, metadata_str, manual_marker));
        grid(ui_data.ax, 'on');
        hold(ui_data.ax, 'off');

        set(ui_data.neuron_text, 'String', sprintf('Neuron: %d / %d', ui_data.current_neuron, ui_data.num_neurons));

        % Plot raw waveform if available
        if ~isempty(ui_data.waveforms_raw)
            cla(ui_data.ax_raw);
            hold(ui_data.ax_raw, 'on');
            if ~isempty(ui_data.timeaxis)
                plot(ui_data.ax_raw, ui_data.timeaxis_raw, ui_data.waveforms_raw(ui_data.current_neuron, :), 'k-', 'LineWidth', 1.5);
                xlabel(ui_data.ax_raw, 'Time (ms)');
            else
                plot(ui_data.ax_raw, ui_data.waveforms_raw(ui_data.current_neuron, :), 'k-', 'LineWidth', 1.5);
                xlabel(ui_data.ax_raw, 'Sample Index');
            end
            ylabel(ui_data.ax_raw, 'Raw Amplitude');
            title(ui_data.ax_raw, 'Raw Waveform');
            grid(ui_data.ax_raw, 'on');
            hold(ui_data.ax_raw, 'off');
        end

        % Update scatter plot
        updateScatterPlot(fig, ui_data, trough_to_peak_val, half_width_val);
    end

    function updateScatterPlot(fig, ui_data, current_trough_to_peak, current_half_width)
        % Debug output
        fprintf('updateScatterPlot called - filtered neurons: %d\n', length(ui_data.filtered_indices));

        % Convert to time if timeaxis available
        if ~isempty(ui_data.timeaxis)
            sample_duration = ui_data.timeaxis_super(2) - ui_data.timeaxis_super(1);
        else
            sample_duration = 1; % Keep as samples
        end

        % Calculate values only for filtered neurons
        num_filtered = length(ui_data.filtered_indices);
        all_trough_to_peak = zeros(num_filtered, 1);
        all_half_width = zeros(num_filtered, 1);

        for idx = 1:num_filtered
            i = ui_data.filtered_indices(idx);
            neuron_key = sprintf('n%d', i);
            if isfield(ui_data.manual_corrections, neuron_key)
                corr = ui_data.manual_corrections.(neuron_key);
                all_trough_to_peak(idx) = (corr.peak_idx - corr.trough_idx) * sample_duration;
                all_half_width(idx) = (corr.idx_ascending - corr.idx_descending) * sample_duration;
            else
                all_trough_to_peak(idx) = ui_data.trough_to_peak(i) * sample_duration;
                all_half_width(idx) = ui_data.half_width(i) * sample_duration;
            end
        end

        % Convert current values to time
        current_trough_to_peak_time = current_trough_to_peak * sample_duration;
        current_half_width_time = current_half_width * sample_duration;

        % Clear and redraw scatter plot
        cla(ui_data.ax_scatter);
        hold(ui_data.ax_scatter, 'on');

        % Plot all neurons with click callback
        h_scatter = scatter(ui_data.ax_scatter, all_trough_to_peak, all_half_width, 30, 'filled', 'MarkerFaceAlpha', 0.6);
        set(h_scatter, 'ButtonDownFcn', @scatterClick);

        % Set mouse down callback on scatter axes
        set(ui_data.ax_scatter, 'ButtonDownFcn', @scatterMouseDown);

        % Plot current neuron with X marker (only if half_width is valid)
        if ~isempty(ui_data.current_marker) && ishandle(ui_data.current_marker)
            delete(ui_data.current_marker);
        end
        if ~isnan(current_half_width)
            ui_data.current_marker = plot(ui_data.ax_scatter, current_trough_to_peak_time, current_half_width_time, 'rx', 'MarkerSize', 20, 'LineWidth', 3);
        end

        if ~isempty(ui_data.timeaxis)
            xlabel(ui_data.ax_scatter, 'Trough-to-Peak (ms)');
            ylabel(ui_data.ax_scatter, 'Half-Width (ms)');
        else
            xlabel(ui_data.ax_scatter, 'Trough-to-Peak (samples)');
            ylabel(ui_data.ax_scatter, 'Half-Width (samples)');
        end
        title(ui_data.ax_scatter, 'Feature Space');
        grid(ui_data.ax_scatter, 'on');
        hold(ui_data.ax_scatter, 'off');

        % Store feature values for click detection
        ui_data.all_trough_to_peak = all_trough_to_peak;
        ui_data.all_half_width = all_half_width;

        % Update stored data
        guidata(fig, ui_data);
    end

    function scatterClick(~, ~)
        ui_data = guidata(gcbf);

        % Don't jump if currently panning
        if ui_data.panning_scatter
            return;
        end

        % Get click position
        pt = get(ui_data.ax_scatter, 'CurrentPoint');
        x_click = pt(1,1);
        y_click = pt(1,2);

        % Find closest neuron
        distances = sqrt((ui_data.all_trough_to_peak - x_click).^2 + (ui_data.all_half_width - y_click).^2);
        [~, closest_idx] = min(distances);

        % Jump to that neuron (convert from filtered index to absolute neuron number)
        ui_data.current_neuron = ui_data.filtered_indices(closest_idx);
        set(ui_data.goto_text, 'String', num2str(ui_data.current_neuron));
        guidata(gcbf, ui_data);
        updatePlot(gcbf);
    end

end
