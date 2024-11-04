% Load the CSV data
% Replace 'pupil_points.csv' with the actual filename of your CSV file
filename = 'MD278_pupilDLC_resnet50_BAfc_pupil_croppedOct25shuffle1_100000.csv';
data = readtable(filename);
data = data(:, 2:end); % Remove the first column containing the row index
frameRate = 60;
% Number of points (change if necessary)
numPoints = 8; 

% Generate new variable names and apply them
newNames = cell(1, numPoints * 3); % For x, y, likelihood columns
for j = 1:numPoints
    newNames{(j - 1) * 3 + 1} = sprintf('x%d', j);         % x1, x2, ..., x8
    newNames{(j - 1) * 3 + 2} = sprintf('y%d', j);         % y1, y2, ..., y8
    newNames{(j - 1) * 3 + 3} = sprintf('likelihood%d', j); % likelihood1, likelihood2, ..., likelihood8
end
data.Properties.VariableNames = newNames;


% Define parameters
numFrames = height(data); % Number of frames
diameters = zeros(numFrames, 1); % To store the diameter for each frame

% Likelihood threshold to filter unreliable points
likelihoodThreshold = 0.9;

% Preallocate the diameters array for speed
diameters = NaN(numFrames, 1); % Initialize with NaN

% Start parallel pool (optional, can also be set in MATLAB preferences)
% parpool;

% Start timer
tic;

% Use parfor for parallel processing
parfor i = 1:numFrames
    x = zeros(numPoints, 1);
    y = zeros(numPoints, 1);
    validPoints = 0;
    
    % Collect valid points with likelihood above threshold
    for j = 1:numPoints
        xVar = sprintf('x%d', j);
        yVar = sprintf('y%d', j);
        likelihoodVar = sprintf('likelihood%d', j);
        
        % Check if these variables exist in the table
        if ismember(xVar, data.Properties.VariableNames) && ...
           ismember(yVar, data.Properties.VariableNames) && ...
           ismember(likelihoodVar, data.Properties.VariableNames)
       
            xCoord = data{i, xVar};
            yCoord = data{i, yVar};
            likelihood = data{i, likelihoodVar};
            
            if likelihood >= likelihoodThreshold
                validPoints = validPoints + 1;
                x(validPoints) = xCoord;
                y(validPoints) = yCoord;
            end
        end
    end
    
    % Fit an ellipse only if we have enough valid points
    if validPoints >= 5 % Minimum points to reliably fit an ellipse
        try
            % Attempt to fit an ellipse to the valid points
            ellipse = fit_ellipse(x(1:validPoints), y(1:validPoints), 'n');
            
            % Calculate the pupil diameter as the major axis of the ellipse
            diameters(i) = max(ellipse.a, ellipse.b) * 2; % Diameter is twice the semi-major axis length
        catch
            % If fit_ellipse fails, set the diameter to NaN for this frame
            diameters(i) = NaN;
            warning('fit_ellipse failed for frame %d, setting diameter to NaN', i);
        end
    else
        % Set diameter to NaN if not enough valid points
        diameters(i) = NaN;
    end
    
    % Update progress
    if mod(i, 10) == 0 % Update every 10 iterations
        elapsedTime = toc; % Get elapsed time
        estimatedTime = (elapsedTime / i) * (numFrames - i); % Estimate remaining time
        fprintf('Processed frame %d of %d. Estimated time left: %.2f seconds.\n', i, numFrames, estimatedTime);
    end
end

% Plot the pupil diameter over time
time = (1:numFrames) / frameRate; % Define `frameRate` based on your video's frame rate
figure;
plot(time, diameters);
xlabel('Time (s)');
ylabel('Pupil Diameter (pixels)');
title('Pupil Diameter Over Time');

% Save diameters to a file
save([filename(1:5) '_pupil_diameters.mat'], 'diameters');