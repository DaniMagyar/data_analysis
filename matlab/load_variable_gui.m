function selected_var_data = load_variable_gui
    % Create a figure window
    fig = figure('Name', 'Load Variable from .mat File', ...
                 'Position', [300, 300, 400, 200], ...
                 'MenuBar', 'none', ...
                 'NumberTitle', 'off', ...
                 'CloseRequestFcn', @close_gui);
    
    % Persistent variable for file path and UI controls
    mat_file_path = '';
    var_dropdown = []; % dropdown for variable names
    load_button = [];  % button to load selected variable
    
    % UI for file selection
    uicontrol('Style', 'text', ...
              'Position', [20, 150, 100, 20], ...
              'String', 'Select .mat file:');
    uicontrol('Style', 'pushbutton', ...
              'Position', [130, 150, 80, 25], ...
              'String', 'Browse', ...
              'Callback', @load_mat_file);
    
    % Dropdown menu for variables
    uicontrol('Style', 'text', ...
              'Position', [20, 100, 100, 20], ...
              'String', 'Select Variable:');
    var_dropdown = uicontrol('Style', 'popupmenu', ...
                             'Position', [130, 100, 200, 25], ...
                             'String', {' '}, ...
                             'Enable', 'off');
    
    % Load button
    load_button = uicontrol('Style', 'pushbutton', ...
                            'Position', [130, 50, 80, 25], ...
                            'String', 'Load', ...
                            'Enable', 'off', ...
                            'Callback', @load_variable);
    
    % Initialize output data variable
    selected_var_data = [];

    % Callback function for loading the .mat file and populating the dropdown
    function load_mat_file(~, ~)
        [file, path] = uigetfile('*.mat', 'Select a .mat file');
        if isequal(file, 0)
            return;
        end
        mat_file_path = fullfile(path, file);
        
        % Get variable information from the file
        mat_info = whos('-file', mat_file_path);
        var_names = {mat_info.name};
        
        % Update dropdown with variable names
        set(var_dropdown, 'String', var_names, 'Enable', 'on');
        
        % Enable the load button
        set(load_button, 'Enable', 'on');
    end

    % Callback function for loading the selected variable
    function load_variable(~, ~)
        if isempty(mat_file_path)
            errordlg('Please select a .mat file first.', 'File Error');
            return;
        end
        
        % Get selected variable name
        var_names = get(var_dropdown, 'String');
        selected_idx = get(var_dropdown, 'Value');
        selected_var_name = var_names{selected_idx};
        
        % Load the selected variable
        mat_data = load(mat_file_path, selected_var_name);
        selected_var_data = mat_data.(selected_var_name);
        
        % Display loaded data in Command Window
        disp('Loaded variable data:');
        disp(selected_var_data);
        
        % Close the GUI
        close(fig);
    end

    % Callback to handle figure close request and output data
    function close_gui(~, ~)
        % Resume execution if the GUI is closed without loading a variable
        uiresume(fig);
        delete(fig);
    end

    % Pause execution until the figure is closed
    uiwait(fig);
end