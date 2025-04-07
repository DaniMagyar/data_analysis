clear all
filename = 'notch';
% Load LFP data
mData = load([filename '.mat']);
numChannels = mData.numChannels;
lfpFile = fopen([filename '.dat'], 'r');
lfpData = fread(lfpFile, [numChannels, Inf], 'int16');
fclose(lfpFile);
samplingRate = mData.fs; % Adjust according to your data
% Load stimulus timestamps
stimData = load('TTLsKS.mat'); 
stimTimestamps = stimData.triptest_sound_only; % Adjust variable name if needed

% Initialize variable to store all spectrograms
allSpectrograms = [];
zscoreallSpectrograms = [];

% Spectrogram parameters for 5 ms time resolution and 1 Hz frequency resolution
windowSize = 10;      % number of samples in a window (for 10ms use 300 at 30kHz sr)
overlap = round(windowSize *0.9);  % 75% overlap (for better time resolution)
nfft = 256;          % frequency resolution = samplingRate/nfft. Should be the closest power of 2 ...
% (e.g. if 10Hz resolution needed at 30kHz 3000 nfft would be the correct, but 2048 or 4096 is optimal), because Fast-Fourier Transform structure

maxFrequency = 5000;    % Maximum frequency to display (Hz)
minFrequency = 0;     % Minimum frequency to display (Hz)
preStim = 0.1;         % 0.5 sec before stimulus
postStim = 0.1;        % 0.5 sec after stimulus
loadBar = waitbar(0,'Calculating...');
for ch = 1:numChannels
    chSpectrograms = [];
    waitbar(ch/numChannels, loadBar);
    for stim = 1:length(stimTimestamps)
        stimIdx = round(stimTimestamps(stim) * samplingRate);
        
        % Define the window to include 0.5 sec before and after the stimulus
        winIdx = stimIdx - round(preStim * samplingRate) : stimIdx + round(postStim * samplingRate);
        
        % Check if the window is within bounds
        if max(winIdx) <= size(lfpData, 2) && min(winIdx) > 0
            lfpSegment = double(lfpData(ch, winIdx));
            
            if ~any(isnan(lfpSegment))
                [s, f, t] = spectrogram(lfpSegment, windowSize, overlap, nfft, samplingRate);
                
                % Select only frequencies between 20 and 120 Hz
                freqIdx = f >= minFrequency & f <= maxFrequency;
                chSpectrograms = cat(3, chSpectrograms, abs(s(freqIdx, :))); %129x242x50
            end
        end
    end
    avgSpectrogram = mean(chSpectrograms, 3); %129x242
    
    % Calculate the row-wise Z-score for the spectrogram
    zscoreSpectrogram = (avgSpectrogram - mean(avgSpectrogram, 2)) ./ std(avgSpectrogram, 0, 2);   
    zscoreallSpectrograms = cat(4, zscoreallSpectrograms, zscoreSpectrogram);

    allSpectrograms = cat(4, allSpectrograms, avgSpectrogram);


end
close(loadBar)
figure
tiledlayout(1,2)

% Average across channels for Z-score transformed spectrogram
avgSpectrogramAllChannels = mean(allSpectrograms, 4);
% Plot average Z-score transformed spectrogram for frequencies between 20 and 120 Hz
nexttile
imagesc(t, f(freqIdx), avgSpectrogramAllChannels);
axis xy;
xline(0.1, '--r', 'LineWidth', 1); 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Average Z-Score Transformed Spectrogram (0.5s Before and 0.5s After Stimuli)');
colorbar;
title('Spectrogram no CAR')

% Average across channels for Z-score transformed spectrogram
avgSpectrogramAllChannels = mean(zscoreallSpectrograms, 4);
% Plot average Z-score transformed spectrogram for frequencies between 20 and 120 Hz
nexttile
imagesc(t, f(freqIdx), avgSpectrogramAllChannels);
axis xy;
xline(0.1, '--r', 'LineWidth', 1); 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Average Z-Score Transformed Spectrogram (0.5s Before and 0.5s After Stimuli)');
colorbar;
title('Z-scored Spectrogram no CAR')


