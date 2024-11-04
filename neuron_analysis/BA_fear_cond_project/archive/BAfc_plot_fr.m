function [avg_fr] = BAfc_plot_fr(varargin)

cell_metrics = BAfc_load_neurons;

prs =  inputParser;
addParameter(prs,'selVariable','none',@ischar) % Variable used for neuron selection from cell_metrics. (e.g. 'labels')
addParameter(prs,'selValue','none',@ischar) % Value used for neuron selection from 'selectVariable'. (e.g. '2')
addParameter(prs,'neuronType','none',@ischar) % Options: 'Narrow Interneuron', 'Wide Interneuron', 'Pyramidal Cell'. Default 'none' loads all types. 
parse(prs,varargin{:})
g = prs.Results;

% Select neurons from cell_metrics based on structure variable
if strcmp(g.selVariable, 'none')
    cellIdx_selVariable = 1:numel(cell_metrics.cellID); % all neurons
else
    cellIdx_selVariable = find(strcmp(cell_metrics.(g.selVariable), g.selValue)); 
end
% Select neuron type
if strcmp(g.neuronType, 'none')
    cellIdx_type = 1:numel(cell_metrics.cellID); % all neurons
else
    cellIdx_type = find(strcmp(cell_metrics.putativeCellType, g.neuronType)); 
end
% Remove 'junk' 
cellIdx_keep =  1:numel(cell_metrics.cellID); % all neurons
cellIdx_keep(find(strcmp(cell_metrics.putativeCellType, 'junk'))) = []; % neurons labeled as 'junk'
% Find neurons that match selected variable and celltype and not junk
cellIdx_selected = intersect(intersect(cellIdx_selVariable, cellIdx_type), cellIdx_keep);

for ii = cellIdx_selected
    AP = cell_metrics.spikes.times{ii};
    time_active = AP(end) - AP(1);
    avg_fr(find(cellIdx_selected==ii),1) = numel(AP) / time_active;
end



% y = histogram(log(avg_fr),10,'Normalization','probability');


% h.YLim = [0 0.35];
 histfit(log(avg_fr),50); % ez akkor felel meg a histogram()-nak ha a nbins egyenlo
% data = log(avg_fr);
% Fit to a normal (or a different distribution if you choose)
% pd = fitdist(data,'normal');
% Find the pdf that spans the disribution
% x_pdf = linspace(min(data),max(data));
% y_pdf = pdf(pd,x_pdf);
% Plot
% histogram(data,10,'Normalization','pdf')
% hold on
% line(x_pdf,y_pdf,'LineWidth',2, 'Color', 'r')
% hold off

h = gca;
h.XLim = [-6 6];
h.XTick = [log(0.1) log(1) log(10) log(100)];
h.XTickLabels = exp(h.XTick);
% figure()
% myhist = histogram(avg_fr,2000);
% set(gca,'xscale','log')
% ylim([0 25])
% figure()

% sorted_fr = sort(avg_fr);
% pd  = fitdist(sorted_fr, 'normal');
% pd = makedist('normal', 'mu', pd.mu, 'sigma', pd.sigma);
% y = pdf(pd, sorted_fr);
% plot(sorted_fr, y)
%set(gca,'xscale','log')
