function [trough_to_peak, half_width, waveforms_super, timeaxis_super] = calculate_waveform_features(waveforms, timeaxis, use_supersampling, supersample_factor, baseline_samples)
% CALCULATE_WAVEFORM_FEATURES Calculate trough-to-peak and half-width for waveforms
%
% Inputs:
%   waveforms - Matrix of waveforms (neurons x samples)
%   timeaxis - Time axis for waveforms
%   use_supersampling - Boolean to enable supersampling
%   supersample_factor - Factor for supersampling
%   baseline_samples - Number of samples before trough for baseline
%
% Outputs:
%   trough_to_peak - Trough-to-peak values (samples)
%   half_width - Half-width values (samples)
%   waveforms_super - Supersampled waveforms
%   timeaxis_super - Supersampled time axis

num_neurons = size(waveforms, 1);
num_samples = size(waveforms, 2);

if use_supersampling
    timeaxis_super = linspace(timeaxis(1), timeaxis(end), num_samples * supersample_factor);
    waveforms_super = zeros(num_neurons, length(timeaxis_super));
    for i = 1:num_neurons
        waveforms_super(i, :) = interp1(timeaxis, waveforms(i, :), timeaxis_super, 'spline');
    end
else
    timeaxis_super = timeaxis;
    waveforms_super = waveforms;
end

% Calculate trough-to-peak and half-width for each neuron
trough_to_peak = zeros(num_neurons, 1);
half_width = zeros(num_neurons, 1);

for i = 1:num_neurons
    [trough_val, trough_idx] = min(waveforms_super(i, :));
    [peak_val, peak_idx] = max(waveforms_super(i, trough_idx:end));
    peak_idx = peak_idx + trough_idx - 1;

    trough_to_peak(i) = peak_idx - trough_idx;

    % Calculate baseline as mean of first baseline_samples
    baseline_mean = mean(waveforms_super(i, 1:baseline_samples));

    % Half-width: distance between middle of descending and ascending parts
    half_amplitude = (trough_val + baseline_mean) / 2;

    % Find middle of descending part (before trough)
    idx_descending = find(waveforms_super(i, 1:trough_idx) <= half_amplitude, 1, 'first');

    % Find middle of ascending part (after trough)
    idx_ascending = find(waveforms_super(i, trough_idx:end) >= half_amplitude, 1, 'first') + trough_idx - 1;

    if ~isempty(idx_descending) && ~isempty(idx_ascending)
        half_width(i) = idx_ascending - idx_descending;
    else
        half_width(i) = NaN;
    end
end

end
