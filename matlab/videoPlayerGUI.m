function videoPlayerGUI
    % Create the main figure
    hFig = figure('Position', [100, 100, 700, 400], ...
                   'MenuBar', 'none', ...
                   'Name', 'Video Player', ...
                   'NumberTitle', 'off', ...
                   'CloseRequestFcn', @closeGUI);

    % Create axes for displaying the video
    hAxes = axes('Parent', hFig, 'Units', 'normalized', ...
                  'Position', [0.1, 0.3, 0.8, 0.6]);

    % Create a slider for frame navigation
    hSlider = uicontrol('Style', 'slider', ...
                        'Min', 1, 'Max', 1, 'Value', 1, ...
                        'Position', [100, 50, 400, 20], ...
                        'Callback', @sliderCallback);
    
    % Create a text display for current frame number
    hFrameNum = uicontrol('Style', 'text', ...
                          'Position', [260, 80, 80, 20], ...
                          'String', 'Frame: 1');

    % Create a button to load the video
    uicontrol('Style', 'pushbutton', ...
              'String', 'Open Video', ...
              'Position', [10, 50, 80, 30], ...
              'Callback', @openVideo);

    % Create an editable text box for frame number input
    hFrameInput = uicontrol('Style', 'edit', ...
                            'Position', [10, 100, 80, 30], ...
                            'String', '1');

    % Create a button to jump to the specified frame
    uicontrol('Style', 'pushbutton', ...
              'String', 'Go to Frame', ...
              'Position', [100, 100, 80, 30], ...
              'Callback', @jumpToFrame);

    % Create a button to load frame indices from a .mat file
    uicontrol('Style', 'pushbutton', ...
              'String', 'Load Frame Indices', ...
              'Position', [10, 150, 100, 30], ...
              'Callback', @loadIndices);

    % Create a listbox to display frame indices
    hListBox = uicontrol('Style', 'listbox', ...
                         'Position', [10, 200, 100, 120], ...
                         'Callback', @listBoxCallback);

    % Create an editable text box for slider step adjustment
    hStepInput = uicontrol('Style', 'edit', ...
                           'Position', [520, 50, 50, 20], ...
                           'String', '3');

    % Create a button to set slider step
    uicontrol('Style', 'pushbutton', ...
              'String', 'Set Step', ...
              'Position', [580, 50, 80, 30], ...
              'Callback', @setSliderStep);

    % Initialize video object and frame indices
    videoObj = [];
    frameIndices = [];

    function openVideo(~, ~)
        % Open a file dialog to select a video file
        [filename, pathname] = uigetfile({'*.mp4;*.avi;*.mov', 'Video Files (*.mp4, *.avi, *.mov)'});
        if isequal(filename, 0)
            return; % User canceled
        end

        % Create video object
        videoObj = VideoReader(fullfile(pathname, filename));

        % Set the slider range
        set(hSlider, 'Min', 1, 'Max', videoObj.NumFrames, 'Value', 1);
        setSliderStep(); % Apply initial slider step from input
        updateFrame(1);
    end

    function sliderCallback(~, ~)
        % Get the current frame number from the slider
        frameNum = round(get(hSlider, 'Value'));
        updateFrame(frameNum);
    end

    function jumpToFrame(~, ~)
        % Get the frame number from the editable text box
        frameNum = str2double(get(hFrameInput, 'String'));
        if ~isempty(videoObj) && frameNum >= 1 && frameNum <= videoObj.NumFrames
            set(hSlider, 'Value', frameNum); % Update slider position
            updateFrame(frameNum);
        else
            % Display an error message if the frame number is invalid
            msgbox('Invalid frame number', 'Error', 'error');
        end
    end

    function loadIndices(~, ~)
        % Load frame indices from a .mat file
        [matfile, matpath] = uigetfile('*.mat', 'Select .mat file with frame indices');
        if isequal(matfile, 0)
            return; % User canceled
        end
        matdata = load(fullfile(matpath, matfile));
        
        % Get variable names and let the user select one
        varnames = fieldnames(matdata);
        [indx, tf] = listdlg('PromptString', 'Select a variable:', 'SelectionMode', 'single', 'ListString', varnames);
        if ~tf
            return; % User canceled
        end
        frameIndices = matdata.(varnames{indx});
        
        % Create indexed display for list box (e.g., '1: 100', '2: 200', ...)
        indexedList = arrayfun(@(i, v) sprintf('%d: %d', i, v), (1:length(frameIndices))', frameIndices, 'UniformOutput', false);
        
        % Populate the list box with indexed frame indices
        set(hListBox, 'String', indexedList);
    end

    function listBoxCallback(~, ~)
        % Jump to the selected frame index from the list box
        selectedIdx = get(hListBox, 'Value');
        if ~isempty(frameIndices)
            frameNum = frameIndices(selectedIdx);
            if frameNum >= 1 && frameNum <= videoObj.NumFrames
                set(hSlider, 'Value', frameNum); % Update slider position
                updateFrame(frameNum);
            else
                msgbox('Frame index out of range', 'Error', 'error');
            end
        end
    end

    function setSliderStep(~, ~)
        % Get the slider step from the input box
        stepSize = str2double(get(hStepInput, 'String'));
        if ~isempty(stepSize) && stepSize > 0
            % Set the slider step for both small and large steps
            set(hSlider, 'SliderStep', [stepSize/(videoObj.NumFrames-1), stepSize/(videoObj.NumFrames-1)]);
        else
            msgbox('Invalid step size', 'Error', 'error');
        end
    end

    function updateFrame(frameNum)
        if ~isempty(videoObj) && frameNum <= videoObj.NumFrames
            % Read and display the current frame
            img = read(videoObj, frameNum);
            imshow(img, 'Parent', hAxes);
            
            % Update the frame number text
            set(hFrameNum, 'String', ['Frame: ', num2str(frameNum)]);
            set(hFrameInput, 'String', num2str(frameNum)); % Update the editable text box
        end
    end

    function closeGUI(~, ~)
        delete(hFig);
    end
end