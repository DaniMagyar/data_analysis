function g = Gergo_firstspike_latency(g)

% Calculate first spike latencies for each neuron
% Input: g.cell_metrics.spikes.times - cell array of spike timestamps for each neuron
%        g.cell_metrics.general.shockTTL - cell array of stimulation timestamps
% Output: g.cell_metrics.firstspike_mean - mean first spike latency for each neuron
%         g.cell_metrics.firstspike_std - std of first spike latencies for each neuron

% Get number of neurons
num_neurons = length(g.cell_metrics.spikes.times);

% Initialize output arrays
g.cell_metrics.firstspike_mean = nan(num_neurons, 1);
g.cell_metrics.firstspike_std = nan(num_neurons, 1);

% Process each neuron
for neuron_idx = 1:num_neurons
    % Get spike times and stimulation times for current neuron
    spike_times = g.cell_metrics.spikes.times{neuron_idx};
    stim_times = g.cell_metrics.general.shockTTL{neuron_idx};
    
    % Skip if either array is empty
    if isempty(spike_times) || isempty(stim_times)
        fprintf('Warning: Neuron %d has empty spike times or stimulation times\n', neuron_idx);
        continue;
    end
    
    % Initialize array to store first spike latencies for this neuron
    first_spike_latencies = [];
    
    % For each stimulation, find the first spike after it
    for stim_idx = 1:length(stim_times)
        stim_time = stim_times(stim_idx);
        
        % Find spikes after this stimulation
        spikes_after_stim = spike_times(spike_times > stim_time);
        
        % If there are spikes after stimulation, record the latency to first spike
        if ~isempty(spikes_after_stim)
            first_spike_latency = spikes_after_stim(1) - stim_time;
            first_spike_latencies = [first_spike_latencies, first_spike_latency];
        end
    end
    
    % Calculate mean and std of first spike latencies for this neuron
    if ~isempty(first_spike_latencies)
        g.cell_metrics.firstspike_mean(neuron_idx) = mean(first_spike_latencies);
        g.cell_metrics.firstspike_std(neuron_idx) = std(first_spike_latencies);
    else
        fprintf('Warning: Neuron %d has no spikes after any stimulation\n', neuron_idx);
    end
end