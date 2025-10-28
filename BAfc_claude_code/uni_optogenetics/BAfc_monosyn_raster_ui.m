function BAfc_monosyn_raster_ui(g, monosyn_results, ttl)
% BAfc_monosyn_raster_ui - Interactive UI to view raster plots for each neuron
%
% Displays 4 raster plots (one for each TTL type) for each PN neuron
% Allows stepping through neurons one by one
%
% Inputs:
%   g - parameters structure from BAfc_monosyn_optogenetics
%   monosyn_results - results structure from BAfc_monosyn_optogenetics
%   ttl - cell array of TTL types

% Get list of all PN neurons
ttl_fieldname = matlab.lang.makeValidName(ttl{1});
PN_indices = monosyn_results.(ttl_fieldname).all_neuron_idx;
num_neurons = length(PN_indices);

% Create figure with wider width for list
fig = figure('Name', 'Monosynaptic Response Raster Viewer', 'Position', [50 50 1350 900]);

% Create UI data structure
ui_data.current_neuron_idx = 1;
ui_data.num_neurons = num_neurons;
ui_data.PN_indices = PN_indices;
ui_data.g = g;
ui_data.monosyn_results = monosyn_results;
ui_data.ttl = ttl;

% Create UI controls
ui_data.prev_btn = uicontrol('Style', 'pushbutton', 'String', 'Previous', ...
    'Position', [20 20 80 30], 'Callback', @prevNeuron);
ui_data.next_btn = uicontrol('Style', 'pushbutton', 'String', 'Next', ...
    'Position', [110 20 80 30], 'Callback', @nextNeuron);
ui_data.neuron_text = uicontrol('Style', 'text', ...
    'String', sprintf('Neuron: %d / %d', ui_data.current_neuron_idx, ui_data.num_neurons), ...
    'Position', [200 20 600 30], 'FontSize', 9);
ui_data.goto_label = uicontrol('Style', 'text', 'String', 'Go to Neuron_idx:', ...
    'Position', [810 20 100 30], 'FontSize', 9, 'HorizontalAlignment', 'right');
ui_data.goto_text = uicontrol('Style', 'edit', 'String', num2str(PN_indices(1)), ...
    'Position', [920 20 60 30], 'FontSize', 10);
ui_data.goto_btn = uicontrol('Style', 'pushbutton', 'String', 'Go To', ...
    'Position', [990 20 60 30], 'Callback', @gotoNeuronIdx);

% Create 4 subplot axes (2x2 grid) - shifted left to make room for list
ui_data.ax1 = axes('Position', [0.05 0.55 0.35 0.38]);
ui_data.ax2 = axes('Position', [0.42 0.55 0.35 0.38]);
ui_data.ax3 = axes('Position', [0.05 0.12 0.35 0.38]);
ui_data.ax4 = axes('Position', [0.42 0.12 0.35 0.38]);

% Create listbox for responsive neurons on the right side
ui_data.list_label = uicontrol('Style', 'text', 'String', 'Responsive Neurons:', ...
    'Position', [1060 850 260 20], 'FontSize', 10, 'FontWeight', 'bold');

% Build list of responsive neurons
responsive_list = {};
responsive_indices = [];
for i = 1:num_neurons
    neuron_idx = PN_indices(i);
    is_responsive = false;

    % Check if responsive to any TTL type
    for tt = 1:length(ttl)
        ttl_fn = matlab.lang.makeValidName(ttl{tt});
        if ismember(neuron_idx, monosyn_results.(ttl_fn).neuron_idx)
            is_responsive = true;
            break;
        end
    end

    if is_responsive
        brain_region = g.cell_metrics.brainRegion{neuron_idx};
        animal = g.cell_metrics.animal{neuron_idx};
        putative_type = g.cell_metrics.putativeCellType{neuron_idx};
        responsive_list{end+1} = sprintf('Idx:%d | %s | %s | %s', neuron_idx, brain_region, putative_type, animal);
        responsive_indices(end+1) = i;
    end
end

ui_data.responsive_indices = responsive_indices;
ui_data.listbox = uicontrol('Style', 'listbox', ...
    'String', responsive_list, ...
    'Position', [1060 100 260 750], ...
    'FontSize', 9, ...
    'Callback', @listboxCallback);

% Store data in figure
guidata(fig, ui_data);

% Initial plot
updatePlot(fig);

    function prevNeuron(src, ~)
        ui_data = guidata(src);
        if ui_data.current_neuron_idx > 1
            ui_data.current_neuron_idx = ui_data.current_neuron_idx - 1;
        end
        neuron_idx = ui_data.PN_indices(ui_data.current_neuron_idx);
        set(ui_data.goto_text, 'String', num2str(neuron_idx));
        guidata(src, ui_data);
        updatePlot(gcbf);
    end

    function nextNeuron(src, ~)
        ui_data = guidata(src);
        if ui_data.current_neuron_idx < ui_data.num_neurons
            ui_data.current_neuron_idx = ui_data.current_neuron_idx + 1;
        end
        neuron_idx = ui_data.PN_indices(ui_data.current_neuron_idx);
        set(ui_data.goto_text, 'String', num2str(neuron_idx));
        guidata(src, ui_data);
        updatePlot(gcbf);
    end

    function gotoNeuronIdx(src, ~)
        ui_data = guidata(src);
        neuron_idx_target = str2double(get(ui_data.goto_text, 'String'));

        % Find this neuron_idx in PN_indices
        idx_pos = find(ui_data.PN_indices == neuron_idx_target, 1);

        if isempty(idx_pos)
            warndlg(sprintf('Neuron_idx %d not found in PN neuron list', neuron_idx_target), 'Invalid Input');
            return;
        end

        ui_data.current_neuron_idx = idx_pos;
        guidata(src, ui_data);
        updatePlot(gcbf);
    end

    function listboxCallback(src, ~)
        ui_data = guidata(src);
        selected_idx = get(src, 'Value');

        % Get the actual neuron index from responsive_indices
        ui_data.current_neuron_idx = ui_data.responsive_indices(selected_idx);

        % Update goto text
        neuron_idx = ui_data.PN_indices(ui_data.current_neuron_idx);
        set(ui_data.goto_text, 'String', num2str(neuron_idx));

        guidata(src, ui_data);
        updatePlot(gcbf);
    end

    function updatePlot(fig)
        ui_data = guidata(fig);

        % Get current neuron index in cell_metrics
        neuron_idx = ui_data.PN_indices(ui_data.current_neuron_idx);

        % Update neuron text
        set(ui_data.neuron_text, 'String', ...
            sprintf('Neuron: %d / %d | Neuron_idx: %d | Cell ID: %d | Animal: %s', ...
            ui_data.current_neuron_idx, ui_data.num_neurons, neuron_idx, ...
            ui_data.g.cell_metrics.cellID(neuron_idx), ...
            ui_data.g.cell_metrics.animal{neuron_idx}));

        % Plot rasters for each TTL type
        axes_list = {ui_data.ax1, ui_data.ax2, ui_data.ax3, ui_data.ax4};

        for tt = 1:length(ui_data.ttl)
            ax = axes_list{tt};
            cla(ax);
            hold(ax, 'on');

            % Get spike times for this neuron and TTL type
            spike_times = ui_data.g.cell_metrics.spikes.times{neuron_idx};

            ttl_events_raw = ui_data.g.cell_metrics.general.(ui_data.ttl{tt});

            % Extract TTL events - they are per neuron in a cell array
            if iscell(ttl_events_raw) && length(ttl_events_raw) == length(ui_data.g.cell_metrics.cellID)
                % TTL events are stored per neuron
                ttl_events = ttl_events_raw{neuron_idx};
            elseif iscell(ttl_events_raw)
                % TTL events are stored per session - find which session
                neuron_session = ui_data.g.cell_metrics.sessionName{neuron_idx};
                session_idx = find(strcmp(ui_data.g.cell_metrics.general.basenames, neuron_session));
                if isempty(session_idx)
                    session_idx = 1;
                end
                ttl_events = ttl_events_raw{session_idx};
            else
                ttl_events = ttl_events_raw;
            end

            if isempty(ttl_events)
                title(ax, sprintf('%s - No Events', ui_data.ttl{tt}));
                continue;
            end

            % Use shorter time window for display (50ms baseline, 50ms post)
            plot_pre_time = g.params.pre_time_short;  % 50ms baseline
            plot_post_time = g.params.post_time_short;  % 50ms post

            % Plot raster
            for trial = 1:length(ttl_events)
                event_time = ttl_events(trial);
                trial_spikes = spike_times - event_time;
                trial_spikes = trial_spikes(trial_spikes >= -plot_pre_time & ...
                                           trial_spikes <= plot_post_time);

                % Plot spikes as vertical lines
                for spike = 1:length(trial_spikes)
                    plot(ax, trial_spikes(spike) * 1000, trial, 'k.', 'MarkerSize', 8);
                end
            end

            % Draw monosynaptic window
            fill(ax, [0 ui_data.g.params.monosyn_window*1000 ui_data.g.params.monosyn_window*1000 0], ...
                [0 0 length(ttl_events)+1 length(ttl_events)+1], ...
                'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

            % Draw stimulus onset line
            plot(ax, [0 0], [0 length(ttl_events)+1], 'r--', 'LineWidth', 2);

            % Draw latency line if responsive
            ttl_fn = matlab.lang.makeValidName(ui_data.ttl{tt});
            idx_in_results_pre = find(ui_data.monosyn_results.(ttl_fn).all_neuron_idx == neuron_idx);
            if ~isempty(idx_in_results_pre)
                zscore_pre = ui_data.monosyn_results.(ttl_fn).all_peak_zscore(idx_in_results_pre);
                prob_pre = ui_data.monosyn_results.(ttl_fn).all_response_probability(idx_in_results_pre);
                latency_pre = ui_data.monosyn_results.(ttl_fn).all_mean_latency(idx_in_results_pre);

                % Only draw if responsive (two-rule system)
                rule1_pre = zscore_pre >= ui_data.g.params.zscore_threshold && ...
                           prob_pre >= ui_data.g.params.prob_threshold;
                rule2_pre = zscore_pre >= ui_data.g.params.zscore_threshold_strict && ...
                           prob_pre >= ui_data.g.params.prob_threshold_lenient;
                is_responsive_pre = rule1_pre || rule2_pre;
                if is_responsive_pre && ~isnan(latency_pre)
                    plot(ax, [latency_pre latency_pre], [0 length(ttl_events)+1], 'b-', 'LineWidth', 1);
                end
            end

            xlim(ax, [-plot_pre_time*1000 plot_post_time*1000]);
            ylim(ax, [0 length(ttl_events)+1]);

            % Only show xlabel for bottom row
            if tt > 2
                xlabel(ax, 'Time from stimulus (ms)');
            end
            ylabel(ax, 'Trial');

            % Get response metrics for this neuron
            idx_in_results = find(ui_data.monosyn_results.(ttl_fn).all_neuron_idx == neuron_idx);
            if ~isempty(idx_in_results)
                zscore = ui_data.monosyn_results.(ttl_fn).all_peak_zscore(idx_in_results);
                prob = ui_data.monosyn_results.(ttl_fn).all_response_probability(idx_in_results);
                latency = ui_data.monosyn_results.(ttl_fn).all_mean_latency(idx_in_results);

                % Count spikes in monosynaptic window across all trials
                total_spikes_in_window = 0;
                for trial = 1:length(ttl_events)
                    event_time = ttl_events(trial);
                    trial_spikes = spike_times - event_time;
                    spikes_in_window = trial_spikes(trial_spikes > 0 & trial_spikes <= ui_data.g.params.monosyn_window);
                    total_spikes_in_window = total_spikes_in_window + length(spikes_in_window);
                end

                % Check if responsive (two-rule system)
                rule1 = zscore >= ui_data.g.params.zscore_threshold && ...
                       prob >= ui_data.g.params.prob_threshold;
                rule2 = zscore >= ui_data.g.params.zscore_threshold_strict && ...
                       prob >= ui_data.g.params.prob_threshold_lenient;
                is_responsive = rule1 || rule2;

                if is_responsive
                    responsive_str = 'RESPONSIVE';
                    title_color = 'r';
                else
                    responsive_str = 'Non-responsive';
                    title_color = 'k';
                end

                title(ax, sprintf('%s\nZ=%.1f, Spks=%d, P=%.2f, Lat=%.1fms - %s', ...
                    ui_data.ttl{tt}, zscore, total_spikes_in_window, prob, latency, responsive_str), ...
                    'Color', title_color, 'Interpreter', 'none');
            else
                title(ax, ui_data.ttl{tt}, 'Interpreter', 'none');
            end

            grid(ax, 'on');
            hold(ax, 'off');
        end
    end

end
