clear all
mainFolder = 'Z:\HajosLab\Dani\Magyar Dániel\Analysis2\Matlab_files';
cd(mainFolder)
filename = 'md190510_bal_9230.mat';
load(filename)
rawdata = bal_9230_Ch2.values;
sampling_rate = 20000;
example_start = 27.170; % in seconds
example_end = 27.270; % in seconds
example_period = rawdata(example_start*sampling_rate+1:example_end*sampling_rate);
timewindow1 = linspace(0,example_end-example_start,(example_end-example_start)*sampling_rate+1);

filename2 = 'md190510_bal_9230_toltes.mat';
load(filename2)
filldata = bal_9230_t_lt_s_Ch2.values;
ttldata = bal_9230_t_lt_s_Ch3.values*1000;
fill_start = 1169; % in seconds
fill_end = 1171; % in seconds
fill_period = filldata(fill_start*sampling_rate+1:fill_end*sampling_rate);
ttl_period = ttldata(fill_start*sampling_rate+1:fill_end*sampling_rate)-5;
timewindow2 = linspace(0,fill_end-fill_start,(fill_end-fill_start)*sampling_rate);

fig = figure('Position', [400, 100, 800, 600]);
t = tiledlayout(fig,2,1,'TileSpacing', 'compact', 'Padding', 'none');

ax = nexttile;
plot(timewindow1,example_period, 'color',[0, 0.6, 0], 'LineWidth',2)
ylabel('Amplitude (mV)')
ylim([-0.4 0.6])
xlabel('Time (s)')
xticks([0 0.025 0.05 0.075 0.1])
box off
ax.FontSize = 12;
ax = nexttile;
hold on
yyaxis left
plot(timewindow2,fill_period, 'color',[0, 0.6, 0], 'LineWidth',2)
ylabel('Amplitude (mV)')
%ylim([-0.4 0.6])
xlabel('Time (s)')
xticks([0 0.5 1 1.5 2])
box off
ylim([-0.8 1.2])
xlim([timewindow2(1) timewindow2(end)])
yyaxis right
%a = area(timewindow2,ttl_period, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
h_fill = fill(timewindow2,ttl_period, [0.2 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
box off
ylim([-0.6 1])
ax.FontSize = 12;
ax.YAxis(1).Color = 'k';  % 'k' is shorthand for black
ax.YAxis(2).Visible = 'off';  % Hide the right y-axis
legend(ax, h_fill, {'Juxtacellular labeling'}, 'Location', 'northwest');  % or 'eastoutside', etc.




exportgraphics(gcf, 'juxta_fig.tiff', 'Resolution', 300);