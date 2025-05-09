cell_metrics = BAfc_load_neurons('NP_BAfc_triptest');
responseTypes = {'exc', 'inh'};
smoothwindow = 20; % excitaciora kb ugyanolyan eredmenyt ad a 20 es az 50-es zscore is, tehat jobb a 20-at hasznalni
for rT = 1:numel(responseTypes)
    brainRegions = {'LA','BA'};
    for bR = 1:numel(brainRegions)    
        %% Calculate response magnitude (z-score)
        % Find the index or responseive units
        [~, idx_shocks, ~] =   BAfc_find_response('NP_BAfc_triptest', 'shock','TTL_triptest_shocks_only',       responseTypes{rT}, 0.5, 0.5, [0.012 0.05], 0.001, 'artefactLength', 0.012, 'brainRegion', brainRegions{bR});
        [~, idx_sounds, ~] =   BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_sound_only',    responseTypes{rT}, 0.5, 0.5, [0.012 0.05], 0.001, 'artefactLength',     0, 'brainRegion', brainRegions{bR});
        [~, idx_both, ~] =     BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_both',          responseTypes{rT}, 0.5, 0.5, [0.012 0.05], 0.001, 'artefactLength', 0.012, 'brainRegion', brainRegions{bR});
        
        idx_all = sort(intersect(idx_shocks,idx_sounds));
        % calculate the zscore change in 0.1 s after the stim. Baseline time 1 sec, bin time 0.1 sec
        [resp_shocks, ~, ~] =  BAfc_find_response('NP_BAfc_triptest', 'shock', 'TTL_triptest_shocks_only',     responseTypes{rT}, 1, 0.1, [0.012 0.1], 0.1); % smoothing does not affect psth_spx
        [resp_sounds, ~, ~] =  BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_sound_only',   responseTypes{rT}, 1, 0.1, [0.012 0.1], 0.1);
        [resp_both, ~, ~] =    BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_both',         responseTypes{rT}, 1, 0.1, [0.012 0.1], 0.1);
        z_resp_shocks = zscore(resp_shocks,0,2);
        z_resp_sounds = zscore(resp_sounds,0,2);
        z_resp_both = zscore(resp_both,0,2);
        fig = figure;
        mainTile = tiledlayout(fig,1,7); % 1 row, 3 column
        ax = nexttile(mainTile); %Line plot with points
        hold on;
        for i = idx_all
            % Plot individual data points
            plot([1, 2, 3], [z_resp_sounds(i,end), z_resp_shocks(i,end), z_resp_both(i,end)], 'k-'); % Line connecting paired points
            plot(1, z_resp_sounds(i,end), 'ro'); % Data point from vector1
            plot(2, z_resp_shocks(i,end), 'bo'); % Data point from vector2
            plot(3, z_resp_both(i,end), 'go'); % Data point from vector2
        end
        hold off;
        title(ax,[brainRegions{bR} ' Lineplots of response magnitude']);
        % Combine the data into one matrix
        data = [z_resp_sounds(idx_all,end), z_resp_shocks(idx_all,end), z_resp_both(idx_all,end)];
        % Create boxplot
        ax = nexttile(mainTile); %Line plot with points
        boxplot(data, 'Labels', {'Vector 1', 'Vector 2', 'Vector 3'});
        title([brainRegions{bR} ' Boxplots of response magnitude']);
        ylabel('Values');

        %% Calculate excitatory response latencies (ms)
        [lat_shocks] =  BAfc_calc_latency('NP_BAfc_triptest', 'TTL_triptest_shocks_only',    1, 0.05, 0.001, 'abs', 1);
        [lat_sounds] =  BAfc_calc_latency('NP_BAfc_triptest', 'TTL_triptest_sound_only',     1, 0.05, 0.001, 'abs', 1);
        [lat_both] =    BAfc_calc_latency('NP_BAfc_triptest', 'TTL_triptest_both',           1, 0.05, 0.001, 'abs', 1);
        ax = nexttile(mainTile); %Line plot with points %Line plot with points
        hold on;
        for i = idx_all
            % Plot individual data points
            plot([1, 2, 3], [lat_sounds(i), lat_shocks(i), lat_both(i)], 'k-'); % Line connecting paired points
            plot(1, lat_sounds(i), 'ro'); % Data point from vector1
            plot(2, lat_shocks(i), 'bo'); % Data point from vector2
            plot(3, lat_both(i), 'go'); % Data point from vector2
        end
        hold off;
        title([brainRegions{bR} ' Lineplots of latencies']);
        % Combine the data
        data2 = {lat_sounds(idx_sounds), lat_shocks(idx_shocks), lat_both(idx_both)};
        % Determine the maximum length among the data arrays
        maxLength = max(cellfun(@length, data2));
        % Pad shorter arrays with NaNs
        paddedData = nan(maxLength, length(data2)); % Initialize with NaNs
        for i = 1:length(data2)
            len = length(data2{i});
            paddedData(1:len, i) = data2{i}; % Fill with actual data
        end
        % Create boxplot
        ax = nexttile(mainTile); %Line plot with points;
        paddedData(find(paddedData<10)) = NaN;
        h = boxplot(paddedData, {'Sounds', 'Shocks', 'Both'}, 'Symbol', '', 'Widths', 0.3);
        set(findobj(h, 'Tag', 'Box'), 'Color', 'k'); % Set box color to black
        set(h,{'linew'},{2})
        ylim([0 50])
        ylabel('Latency (ms)');
        title([brainRegions{bR} ' Latency Distributions']);
        % Histograms
        stimnames = {'Sound','Shock','Both'};
        for ii = 1:numel(stimnames)
            ax = nexttile(mainTile);
            histogram(data2{ii}(find(data2{ii}>10)),'BinWidth',2) % remove <3ms responses
            xlim([0 50])
            ylim([0 30])
            xlabel('Time(ms)')
            ylabel('Number of neurons')
            set(gca,'XTick',0:20:100);
            set(gca,'TickLength',[0.02, 0.02])
            set(gca, 'TickDir', 'out')
            set(gca,'box','off')
            title([stimnames{ii} ' latency(ms)'])
        end
        set(findall(gcf,'-property','FontSize'),'FontSize',15)

        responses.(brainRegions{bR}).shocks.(responseTypes{rT}) = idx_shocks;
        responses.(brainRegions{bR}).sounds.(responseTypes{rT}) = idx_sounds;
    end
end

cell_metrics.ResponseShocks(1:numel(cell_metrics.cellID)) = 0;
cell_metrics.ResponseSounds(1:numel(cell_metrics.cellID)) = 0;
cell_metrics.ResponseShocks([responses.LA.shocks.exc;responses.BA.shocks.exc]) = 1;
cell_metrics.ResponseShocks([responses.LA.shocks.inh;responses.BA.shocks.inh]) = -1;
cell_metrics.ResponseSounds([responses.LA.sounds.exc;responses.BA.sounds.exc]) = 1;
cell_metrics.ResponseSounds([responses.LA.sounds.inh;responses.BA.sounds.inh]) = -1;

CellExplorer('metrics',cell_metrics);


disp(['sound excited: ' num2str(sum(cell_metrics.ResponseSounds==1))])
disp(['sound inhibited: ' num2str(sum(cell_metrics.ResponseSounds==-1))])
disp(['shock excited: ' num2str(sum(cell_metrics.ResponseShocks==1))])
disp(['shock inhibited: ' num2str(sum(cell_metrics.ResponseShocks==-1))])