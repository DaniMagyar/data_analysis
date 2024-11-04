function [group_averages] = BAfc_plot_pupil_changes_group(animalID)

    % Load data
    %animalID = 'MD278';
    secPre = 10;
    secPost = 20;
    cd C:\Users\dmagyar\Desktop\BAfc_pupil_cropped-DM-2024-10-25\videos
    filename1 = [animalID '_crossing_indices.csv'];
    filename2 = [animalID '_pupil_diameters.mat'];
    filename3 = [animalID '_TTLsKS.mat'];
    filename4 = [animalID '.csv'];
    filename5 = ['C:\Users\dmagyar\Desktop\BA_fear_cond\' animalID '_kilosort\kilosort25preprocess\rez.mat'];
    LED_indices = table2array(readtable(filename1)) + 1; % +1 because it comes from python
    pupil_diameter_data = smoothdata(cell2mat(struct2cell(load(filename2))), 'movmean', 20); % 100fps video -> 100-500 smoothing (chatGTP)
    %plot(pupil_diameter_data)
    title('Smoothed Data')
    TTLsKS = load(filename3);
    TTL_LED = TTLsKS.LED_synch;
    trainFirsts = TTLsKS.sound_all1(1:5:end);
    TTL_3kHz = trainFirsts(1:2:end);
    TTL_12kHz = trainFirsts(2:2:end);
    
    if ismember(TTLsKS.tone_habit_first, TTL_3kHz)
        TTL_tone = TTL_3kHz;
        TTL_noise = TTL_12kHz;
    elseif ismember(TTLsKS.tone_habit_first, TTL_12kHz)
        TTL_tone = TTL_12kHz;
        TTL_noise = TTL_3kHz;
    end
    
    if numel(LED_indices) ~= 40
        error('incorrect number of TTLs detected')
    end
    
    dT_TTL_tone = TTL_LED - TTL_tone;
    dT_TTL_noise = TTL_LED - TTL_noise;
    dT_TTL_shock = TTL_LED(11:30) - TTLsKS.shocks;
    
    
    VideoData = readtable(filename4);
    VideoTS = table2cell(VideoData(:,17));
    
    % Convert the timestamp strings to datetime objects
    timestamps = datetime(VideoTS, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSXXX', 'TimeZone', 'UTC');
    
    % Calculate the time differences in seconds between each consecutive timestamp
    timeDiffs = [0; seconds(diff(timestamps))]; % 0 for length correction
    
    disp('Loading rez.mat')
    rezData = load(filename5);
    
    recordingDiff = diff([sum(timeDiffs) rezData.rez.ops.tend/30000]);
    sampling_rate = numel(timeDiffs)/(rezData.rez.ops.tend/30000);
    disp(['Recording difference:' num2str(recordingDiff)]);
    disp(['est. sampling rate = ' num2str(sampling_rate)])
    
    
    % find shock indices
    for ii = 1:numel(LED_indices(11:30))
        jj = 1;
        rTime = 0;
        while rTime < dT_TTL_shock(ii)
            jj = jj+1;
            rTime = sum(timeDiffs(LED_indices(10+ii)-jj:LED_indices(10+ii)));
        end
        shock_indices(ii,1) = LED_indices(10+ii) - jj;
    end
    
    % tone indices
    for ii = 1:numel(LED_indices)
        jj = 1;
        rTime = 0;
        while rTime < dT_TTL_tone(ii)
            jj = jj+1;
            rTime = sum(timeDiffs(LED_indices(ii)-jj:LED_indices(ii)));
        end
        tone_indices(ii,1) = LED_indices(ii) - jj;
    end
    
    TTL_tone_diff = diff(TTL_tone);
    for kk = 1:numel(TTL_tone_diff)-1
        VID_tone_diff(kk,1) = sum(timeDiffs(tone_indices(kk):tone_indices(kk+1)));
    end
    
    % noise indices        
    for ii = 1:numel(LED_indices)
        jj = 1;
        rTime = 0;
        while rTime < dT_TTL_noise(ii)
            jj = jj+1;
            rTime = sum(timeDiffs(LED_indices(ii)-jj:LED_indices(ii)));        
        end
        noise_indices(ii,1) = LED_indices(ii) - jj;
    end
    
    TTL_noise_diff = diff(TTL_noise);
    for kk = 1:numel(TTL_noise_diff)-1
        VID_noise_diff(kk,1) = sum(timeDiffs(noise_indices(kk):noise_indices(kk+1)));
    end
    
            
    mismatch_tone = TTL_tone_diff(1:38) - VID_tone_diff;        
    mismatch_noise = TTL_noise_diff(1:38) - VID_noise_diff;
    % disp('tone mismatch:   noise mismatch:')
    % disp([mismatch_tone mismatch_noise])
    if any(abs(mismatch_tone)>0.1) || any(abs(mismatch_noise)>0.1)
        error('large mismatch')
    else
        disp('mismatch is small, nice !')
    end
    
    
    stim_indices = tone_indices;
    % Parameters
    pre_stimulus = round(secPre * sampling_rate);   % 5 seconds before stimulus
    post_stimulus = round(secPost * sampling_rate); % 10 seconds after stimulus
    
    % Initialize a matrix to store the segments
    num_signals = length(stim_indices);
    segments = zeros(num_signals, pre_stimulus + post_stimulus);
    
    % Extract pupil diameter segments around each auditory signal
    for i = 1:num_signals
        index = stim_indices(i);
        
        % Extract the pupil diameter segment around each stimulus
        start_idx = max(index - pre_stimulus, 1) + 1;
        end_idx = min(index + post_stimulus, length(pupil_diameter_data));
        
        segments(i, :) = pupil_diameter_data(start_idx:end_idx);
    end
    
    % Grouping
    num_groups = 8;        % Number of groups
    rows_per_group = 5;    % Rows per group
    
    % Preallocate average array
    group_averages = zeros(num_groups, pre_stimulus + post_stimulus);
    
    % Calculate average for each group using 'segments'
    for group_idx = 1:num_groups
        start_row = (group_idx - 1) * rows_per_group + 1;
        end_row = min(group_idx * rows_per_group, num_signals); % Ensure we don't exceed num_signals
        group_averages(group_idx, :) = mean(segments(start_row:end_row, :), 1, 'omitnan'); % Omit NaNs if any
    end
    
    % % Create a UI figure with a scrollable panel
    % f = uifigure('Name', 'Pupil Diameter Segments', 'Position', [100, 100, 800, 600]); % Set window size
    % p = uipanel(f, 'Scrollable', 'on', 'Position', [20, 20, 750, 550]);
    % 
    % % Define layout for multiple subplots with scroll
    % layout = uigridlayout(p, [num_groups, 1]);
    % layout.RowHeight = repmat({'fit'}, 1, num_groups); % Dynamic row height
    % layout.Scrollable = 'on';
    % 
    % % Time axis for plotting (in seconds)
    % time_axis = linspace(-secPre, secPost, pre_stimulus + post_stimulus); 
    % 
    % % Create subplots with larger size in the scrollable panel
    % for group_idx = 1:num_groups
    %     ax = uiaxes(layout); % Create an axis in the layout
    %     plot(ax, time_axis, group_averages(group_idx, :));
    %     title(ax, ['Average Segment - Group ' num2str(group_idx)]);
    %     xlabel(ax, 'Time (s)');
    %     ylabel(ax, 'Pupil Diameter');
    %     ylim(ax, [20 100]);
    %     ax.FontSize = 12; % Increase font size for better readability
    % end
    % 
    % % Set an overall title for the figure
    % sgtitle('Pupil Diameter Segments - Grouped Averages');
    % 
    
    save([animalID '_data.mat'], 'tone_indices', 'noise_indices', 'LED_indices', 'shock_indices', 'sampling_rate')
    disp([animalID '_data.mat saved'])
end