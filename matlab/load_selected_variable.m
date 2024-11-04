function selected_var_data = load_selected_variable(mat_file_path)
    % List variables in the .mat file
    mat_info = whos('-file', mat_file_path);
    
    % Display available variables to the user
    fprintf('Available variables in the .mat file:\n');
    for i = 1:length(mat_info)
        fprintf('%d. %s\n', i, mat_info(i).name);
    end
    
    % Prompt user to select a variable
    selected_idx = input('Enter the number of the variable to load: ');
    
    % Check if input is valid
    if selected_idx < 1 || selected_idx > length(mat_info)
        error('Invalid selection. Please run the function again and select a valid option.');
    end
    
    % Load only the selected variable
    selected_var_name = mat_info(selected_idx).name;
    data = load(mat_file_path, selected_var_name);
    selected_var_data = data.(selected_var_name);
    
    % Display loaded data
    disp('Loaded variable data:');
    disp(selected_var_data);
end