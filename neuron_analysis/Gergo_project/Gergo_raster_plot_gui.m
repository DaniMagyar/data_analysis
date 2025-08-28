function Gergo_raster_plot_gui(g, LA_neurons_sorted_idx)
    % GUI for plotting raster plots of neurons around stimulus times
    % Input:
    %   g - structure containing cell_metrics with spikes.times and general.shockTTL
    %   LA_neurons_sorted_idx - vector of neuron indices to examine
    
    % Validate inputs
    if nargin < 2
        error('Please provide both g structure and LA_neurons_sorted_idx');
    end
    
    % Initialize variables
    current_neuron_idx = 1;
    time_window = 0.2; % +/- 1 second
    
    % Create main figure
    fig = figure('Name', 'Neuron Raster Plot Viewer', ...
                 'NumberTitle', 'off', ...
                 'Position', [100, 100, 1000, 600], ...
                 'MenuBar', 'none', ...
                 'ToolBar', 'figure');
    
    % Create UI controls
    control_panel = uipanel('Parent', fig, ...
                           'Title', 'Controls', ...
                           'Position', [0.02, 0.85, 0.96, 0.12]);
    
    % Current neuron display
    uicontrol('Parent', control_panel, ...
              'Style', 'text', ...
              'String', 'Current Neuron:', ...
              'Position', [20, 40, 100, 20], ...
              'HorizontalAlignment', 'left');
    
    neuron_text = uicontrol('Parent', control_panel, ...
                           'Style', 'text', ...
                           'String', '', ...
                           'Position', [120, 40, 100, 20], ...
                           'HorizontalAlignment', 'left', ...
                           'FontWeight', 'bold');
    
    % Navigation buttons
    prev_btn = uicontrol('Parent', control_panel, ...
                        'Style', 'pushbutton', ...
                        'String', '< Previous', ...
                        'Position', [250, 35, 80, 30], ...
                        'Callback', @prev_neuron);
    
    next_btn = uicontrol('Parent', control_panel, ...
                        'Style', 'pushbutton', ...
                        'String', 'Next >', ...
                        'Position', [340, 35, 80, 30], ...
                        'Callback', @next_neuron);
    
    % Neuron selection
    uicontrol('Parent', control_panel, ...
              'Style', 'text', ...
              'String', 'Go to neuron #:', ...
              'Position', [450, 40, 100, 20], ...
              'HorizontalAlignment', 'left');
    
    neuron_edit = uicontrol('Parent', control_panel, ...
                           'Style', 'edit', ...
                           'String', '1', ...
                           'Position', [550, 40, 50, 25], ...
                           'Callback', @goto_neuron);
    
    total_text = uicontrol('Parent', control_panel, ...
                          'Style', 'text', ...
                          'String', sprintf('/ %d', length(LA_neurons_sorted_idx)), ...
                          'Position', [605, 40, 50, 20], ...
                          'HorizontalAlignment', 'left');
    
    % Time window control
    uicontrol('Parent', control_panel, ...
              'Style', 'text', ...
              'String', 'Time window (±sec):', ...
              'Position', [700, 40, 120, 20], ...
              'HorizontalAlignment', 'left');
    
    time_edit = uicontrol('Parent', control_panel, ...
                         'Style', 'edit', ...
                         'String', '1', ...
                         'Position', [820, 40, 40, 25], ...
                         'Callback', @update_time_window);
    
    % Main plot axis
    main_ax = axes('Parent', fig, ...
                   'Position', [0.08, 0.15, 0.85, 0.65]);
    
    % Info panel
    info_panel = uipanel('Parent', fig, ...
                        'Title', 'Information', ...
                        'Position', [0.02, 0.02, 0.96, 0.10]);
    
    info_text = uicontrol('Parent', info_panel, ...
                         'Style', 'text', ...
                         'String', '', ...
                         'Position', [10, 10, 950, 40], ...
                         'HorizontalAlignment', 'left', ...
                         'Max', 2);
    
    % Plot initial neuron
    plot_current_neuron();
    
    % Callback functions
    function prev_neuron(~, ~)
        if current_neuron_idx > 1
            current_neuron_idx = current_neuron_idx - 1;
            set(neuron_edit, 'String', num2str(current_neuron_idx));
            plot_current_neuron();
        end
    end
    
    function next_neuron(~, ~)
        if current_neuron_idx < length(LA_neurons_sorted_idx)
            current_neuron_idx = current_neuron_idx + 1;
            set(neuron_edit, 'String', num2str(current_neuron_idx));
            plot_current_neuron();
        end
    end
    
    function goto_neuron(~, ~)
        new_idx = str2double(get(neuron_edit, 'String'));
        if isnan(new_idx) || new_idx < 1 || new_idx > length(LA_neurons_sorted_idx)
            set(neuron_edit, 'String', num2str(current_neuron_idx));
            return;
        end
        current_neuron_idx = round(new_idx);
        plot_current_neuron();
    end
    
    function update_time_window(~, ~)
        new_window = str2double(get(time_edit, 'String'));
        if isnan(new_window) || new_window <= 0
            set(time_edit, 'String', num2str(time_window));
            return;
        end
        time_window = new_window;
        plot_current_neuron();
    end
    
    function plot_current_neuron()
        % Get current neuron index
        neuron_idx = LA_neurons_sorted_idx(current_neuron_idx);
        
        % Update neuron display
        set(neuron_text, 'String', sprintf('Index: %d (%d/%d)', ...
            neuron_idx, current_neuron_idx, length(LA_neurons_sorted_idx)));
        
        % Get spike times and stimulus times for this neuron
        spike_times = g.cell_metrics.spikes.times{neuron_idx};
        stim_times = g.cell_metrics.general.shockTTL{neuron_idx};
        
        % Clear the axis
        cla(main_ax);
        hold(main_ax, 'on');
        
        if isempty(stim_times) || isempty(spike_times)
            text(main_ax, 0.5, 0.5, 'No data available for this neuron', ...
                'HorizontalAlignment', 'center', 'Units', 'normalized');
            set(info_text, 'String', 'No spike or stimulus data available');
            return;
        end
        
        % Plot raster
        trial_count = 0;
        spike_count = 0;
        
        for trial = 1:length(stim_times)
            stim_time = stim_times(trial);
            
            % Find spikes within time window
            spike_window = spike_times(spike_times >= (stim_time - time_window) & ...
                                     spike_times <= (stim_time + time_window));
            
            % Convert to relative times
            relative_spikes = spike_window - stim_time;
            
            % Plot spikes for this trial
            if ~isempty(relative_spikes)
                plot(main_ax, relative_spikes, trial * ones(size(relative_spikes)), ...
                     'k|', 'MarkerSize', 8, 'LineWidth', 1);
                spike_count = spike_count + length(relative_spikes);
            end
            trial_count = trial_count + 1;
        end
        
        % Plot stimulus line
        plot(main_ax, [0 0], [0.5 trial_count + 0.5], 'r-', 'LineWidth', 2);
        
        % Set axis properties
        set(main_ax, 'XLim', [-time_window, time_window]);
        set(main_ax, 'YLim', [0.5, trial_count + 0.5]);
        xlabel(main_ax, 'Time relative to stimulus (s)');
        ylabel(main_ax, 'Trial number');
        title(main_ax, sprintf('Neuron %d - Raster Plot', neuron_idx));
        grid(main_ax, 'on');
        
        % Update info
        firing_rate = spike_count / (trial_count * 2 * time_window);
        info_str = sprintf('Trials: %d | Total spikes: %d | Avg firing rate: %.2f Hz | Time window: ±%.1f s', ...
                          trial_count, spike_count, firing_rate, time_window);
        set(info_text, 'String', info_str);
        
        hold(main_ax, 'off');
    end
    
    % Add keyboard shortcuts
    set(fig, 'KeyPressFcn', @key_press);
    
    function key_press(~, event)
        switch event.Key
            case 'leftarrow'
                prev_neuron();
            case 'rightarrow'
                next_neuron();
            case 'space'
                next_neuron();
        end
    end
end