function BAfc_plot_Zscore_all_regions(varargin)
prs =  inputParser;
addParameter(prs,'Sorted',0,@isnumeric) % If 1, Load 'SortIDX' from 'BAparams.mat'. Default 0.
parse(prs,varargin{:})
g = prs.Results;


regions = {'LA', 'BA_med', 'BA_lat', 'CeA', 'Astria'};

allStim = {'TTL_tone_habit_first', 'TTL_noise_habit_first',...
    'TTL_tone_cond_first', 'TTL_noise_cond_first',...
    'TTL_tone_recall_first', 'TTL_noise_recall_first'};

% allStim = {'TTL_tone_habit_all', 'TTL_noise_habit_all',...
%     'TTL_tone_cond_all', 'TTL_noise_cond_all',...
%     'TTL_tone_recall_all', 'TTL_noise_recall_all'};

newcolors = [0.8500 0.3250 0.0980; 0 0.4470 0.7410; 0.9 0.9 0.9];

for jj = 1:6 % Cycle through stims
    currStim = allStim(jj);
    figure(jj)
    tiledlayout(5,2)
    for ii = 1:length(regions)
        nexttile
        if jj == 1
            response_dir = BAfc_plot_Zscore('psth_bin',600, 'int',[-5 5], 'Wcx_win',[-4.5 4.5], 'Stims', currStim,...
                'SelVariable', 'brainRegion', 'selValue', regions{ii});
        elseif jj ~= 1
            response_dir = BAfc_plot_Zscore('psth_bin',600, 'int',[-5 5], 'Wcx_win',[-4.5 4.5], 'Stims', currStim,...
                'SelVariable', 'brainRegion', 'selValue', regions{ii}, 'Sorted', 1);
        end
        title(regions{ii}, 'FontSize',10)
        cb = colorbar('FontSize',10);
        cb.Ticks = linspace(-20,20,41);
        cb.Label.String = 'Z-score change';
        cb.Label.FontSize = 10;
        ax(ii) = nexttile;
        piedata = [numel(response_dir(response_dir==1)),numel(response_dir(response_dir==-1)),...
            numel(response_dir(response_dir==0))];
        labels = {['Excited (' num2str(round(piedata(1)/sum(piedata)*100)) '%)'],...
            ['Inhibited (' num2str(round(piedata(2)/sum(piedata)*100)) '%)'],...
            ['Non-responsive (' num2str(round(piedata(3)/sum(piedata)*100)) '%)']};
        pie(ax(ii), piedata, labels);
    end
    for ii = 1:length(regions)
        ax(ii).Colormap = newcolors;
    end
end


figure(7)
tiledlayout(5,2)
regions = {'LA', 'BA_med', 'BA_lat', 'CeA', 'Astria'};
for ii = 1:length(regions)
    nexttile
    response_dir = BAfc_plot_Zscore('psth_bin',600, 'int',[-0.5 0.5], 'Wcx_win',[-0.2 0.2], 'Stims', ...
        {'TTL_shocks'}, 'SelVariable', 'brainRegion', 'selValue', regions{ii});
    title(regions{ii})
    cb = colorbar('FontSize',10);
    cb.Ticks = linspace(-20,20,41);
    cb.Label.String = 'Z-score change';
    cb.Label.FontSize = 10;
    ax2(ii) = nexttile;
    piedata = [numel(response_dir(response_dir==1)),numel(response_dir(response_dir==-1)),...
        numel(response_dir(response_dir==0))];
    labels = {['Excited (' num2str(round(piedata(1)/sum(piedata)*100)) '%)'],...
        ['Inhibited (' num2str(round(piedata(2)/sum(piedata)*100)) '%)'],...
        ['Non-responsive (' num2str(round(piedata(3)/sum(piedata)*100)) '%)']};
    pie(ax2(ii), piedata, labels)
end
% newcolors = [0.8500 0.3250 0.0980; 0 0.4470 0.7410; 0.9 0.9 0.9];
% ax(1).Colormap = newcolors;
% ax(2).Colormap = newcolors;
% ax(3).Colormap = newcolors;
% ax(4).Colormap = newcolors;
% ax(5).Colormap = newcolors;
% ax2(1).Colormap = newcolors;
% ax2(2).Colormap = newcolors;
% ax2(3).Colormap = newcolors;
% ax2(4).Colormap = newcolors;
% ax2(5).Colormap = newcolors;