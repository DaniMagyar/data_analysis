clear all

cell_metrics = BAfc_load_neurons;
ii = 1;



subsetPlots = [];
TTL = cell_metrics.general.TTL_tone_cond_first{ii};
AP = cell_metrics.spikes.times{ii};
pre_time = 10; % in sec     
post_time = 20;    
bin_time = 1; %in seconds
for jj = 1:numel(TTL) %Each TTL is a a column. 
    preAP{:,jj} = AP(AP>=(TTL(jj)-pre_time) & AP<TTL(jj)); % spikes before each TTL separately
    postAP{:,jj} = AP(AP>TTL(jj) & AP<(TTL(jj)+post_time)); % spikes after each TTL separately
end
for ll = 1:numel(TTL) %Each TTL is a a column. 
    preAP_norm{ll} = preAP{ll}-TTL(ll); % spikes relative to their own TTL
    postAP_norm{ll} = postAP{ll}-TTL(ll);
end
for tt = 1:numel(TTL) %Each TTL is a a column. 
    for nn = 1:(pre_time/bin_time) % number of timestamps in each bin.
         preAP_bin(nn,tt) = sum(preAP_norm{tt}>=(-pre_time+(nn-1)*bin_time) & preAP_norm{tt}<(-pre_time+nn*bin_time));
    end
    for oo = 1:(post_time/bin_time)
         postAP_bin(oo,tt) = sum(postAP_norm{tt}>=((oo-1)*bin_time) & postAP_norm{tt}<(oo*bin_time));
    end
end  
psth_spx_dani = vertcat(sum(preAP_bin,2), sum(postAP_bin,2));

psth_spx_zscore = (psth_spx_dani - mean(psth_spx_dani(1:pre_time/bin_time))) / std(psth_spx_dani(1:pre_time/bin_time));


% Define the x values from -10 to 20 with a matching length to your data
x_values = linspace(-10, 20, length(psth_spx_zscore));

% Plot the bar chart with specified x values
bar(x_values, psth_spx_zscore);
xline(0, 'LineWidth', 1.5, 'Color', 'b')
 xline(4.5, 'LineWidth', 1.5, 'Color', 'r')
% Add labels and title for clarity
xlabel('Time (s)');
ylabel('Z-scored response');
ylim([-2 10])
set(gca, 'FontSize', 14); % Set the font size for the axis labels and ticks



% cell01 ylim([-2 25])
% cell235 ylim([-3 10])
% cell159 ylim([-9 5])