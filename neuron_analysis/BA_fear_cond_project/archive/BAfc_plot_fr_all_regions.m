function BAfc_plot_fr_all_regions

tiledlayout(1,6)
regions = {'LA', 'BA_med', 'BA_lat', 'CeA', 'Astria'};
for ii = 1:length(regions)
    nexttile
    avg_fr = [];
    avg_fr = BAfc_plot_fr('selVariable', 'brainRegion', 'selValue', regions{ii});
    title(regions{ii})
    subtitle( "mean: " + string(mean(avg_fr)) + ",  std: " + string(std(avg_fr)))
    all_regions_fr{ii} = avg_fr;
end
nexttile
color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]};
for jj = 1:length(regions)
    h = histfit(log(all_regions_fr{jj}));
    set(h(2),'color',color{jj})
    delete(h(1))
    hold on
end
h = gca;
h.XLim = [-6 6];
%h.YLim = [0 35];
h.XTick = [log(0.1) log(1) log(10) log(100)];
h.XTickLabels = exp(h.XTick);
legend(regions)