function plot_PsthZscore(Stim, PSTHall, Wilcoxon_results, SortIDX, testMean, preferences)

figure; 
subplot(1,4,1:2)
imagesc(preferences.time/preferences.fs, 1:size(PSTHall,1), PSTHall(SortIDX,:)); % plotting sorted psth matrix 
clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
colormap(preferences.mycolormap); 
hold on;
plot([-0.0125 -0.0125],[1 size(PSTHall,1)],'r', 'LineWidth', 2);
% plot(int,[min(find(testMean(SortIDX)<0))-0.5 min(find(testMean(SortIDX)<0))-0.5], 'r') % finding the first negative Z-socre mean
hold off;
ylabel('# Cell');
xlabel('Time (s)');
set(gca,'FontSize',35);

Wilcox_sorted = str2num(cell2mat(Wilcoxon_results(:,2)));
Wilcox_sorted = Wilcox_sorted(SortIDX); 
significant_idx = find(Wilcox_sorted == 1);
non_significant_idx = find(Wilcox_sorted == 0);
subplot(1,4,3:4)
plot(zeros(1,length(non_significant_idx)), non_significant_idx, 'b.', 'MarkerSize', 30)
hold on
plot(zeros(1,length(significant_idx)), significant_idx, 'r.', 'MarkerSize', 30)
hold off
set(gca,'ydir','reverse')
set(gca,'visible','off')
set(subplot(1,4,3:4), 'Position', [0.5, 0.115, 0.01, 0.805])
xlim([-0.5 0.5])
ylim([1 length(Wilcox_sorted)])

testMeanSorted = testMean(SortIDX);
for qq = 1:length(testMean)
    if Wilcox_sorted(qq) == 1
        if testMeanSorted(qq) < 0
            response_dir(qq) = -1;
        elseif testMeanSorted(qq) >= 0
            response_dir(qq) = 1;
        end
    else
        response_dir(qq) = 0;
    end
end    

annotation('textbox', [0.55, 0.9, 1, 0], 'string',...
    ['Significantly modulated units:' num2str(numel(significant_idx))], 'LineStyle','none', 'FontSize', 35)
annotation('textbox', [0.55, 0.7, 1, 0], 'string',...
    ['Excited: ' num2str(numel(find(response_dir == 1))) '  ('...
    num2str(numel(find(response_dir == 1))/length(testMean)*100) '%)'], 'LineStyle','none', 'FontSize', 35)
annotation('textbox', [0.55, 0.5, 1, 0], 'string',...
    ['Inhibited: ' num2str(numel(find(response_dir == -1))) '  ('...
    num2str(numel(find(response_dir == -1))/length(testMean)*100) '%)'], 'LineStyle','none', 'FontSize', 35)

set(gcf,'Position',[500 50 1600 1300])
% saveas(gcf, [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) '.png'])


%     % Plot the NeuronIDs
%     MyNewOrderStrings=cellfun(@num2str,MyNewOrder,'un',0);
%     for pp = 1:length(MyNewOrderStrings)
%         CurrentID=strjoin(MyNewOrderStrings(pp,:));
%         text(1,pp,[CurrentID(1:5) CurrentID(17:end)], 'FontSize', 8);
%     end
    switch Stim
        case {'BA_25_5Hz', 'BA_25_10Hz', 'BA_250_5Hz', 'shock_only'}
            save ([preferences.mainFolder '\BAparams.mat'], 'SortIDX', 'response_dir')
        case {'TO_25_5Hz', 'TO_25_10Hz', 'TO_250_5Hz', 'shock_inh'}
            load ([preferences.mainFolder '\BAparams.mat'], 'response_dir') % for average plots
    end

%% Plot the averages
if preferences.average == 1

    PSTHallSorted = PSTHall(SortIDX,:);
    Clustmeans(1,:) = mean(PSTHallSorted((find(response_dir==1)),:),1); % excited avg
    Clustmeans(2,:) = mean(PSTHallSorted((find(response_dir==-1)),:),1); % inhibited avg
    Clustmeans(1,:) = Clustmeans(1,:) - mean(Clustmeans(1,1:abs(preferences.int(1)*0.8/(preferences.psth_bin/preferences.fs)))); %offset, 0.8: PV start excluded
    Clustmeans(2,:) = Clustmeans(2,:) - mean(Clustmeans(2,1:abs(preferences.int(1)*0.8/(preferences.psth_bin/preferences.fs)))); %offset
    for tt = 1:2
        figure(tt+1)
        hold on
        plot(preferences.time/preferences.fs+preferences.psth_bin/preferences.fs, Clustmeans(tt,:), 'LineWidth', 2.5) % time shifted (different timing is necessary from heatmap's)
        hold off
        ylabel('Z-scored firing rate');
        xlabel('Time (s)');
        set(gca,'FontSize',35);
        set(gcf,'Position',[500 50 1600 1300])
        if tt == 1
            title('Average response of excited neurons')
        elseif tt == 2
            title('Average response of inhibited neurons')
        end
%         switch Stim
%             case {'BA_25_5Hz', 'shock_only'}
%                 if tt == 1
%                     saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) Stim '_AVG_exc.png'])
%                 elseif tt == 2
%                     saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) Stim '_AVG_inh.png'])
%                 end
%             case {'TO_25_5Hz', 'shock_inh'}
%                 if tt == 1
%                     saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) Stim '_AVG_exc.png'])
%                 elseif tt == 2
%                     saveas(figure(tt+1), [mainFolder '\figures\' matlab.lang.makeValidName(cell2mat(nucleus)) Stim(1:2) num2str(int(2)) Stim '_AVG_inh.png'])
%                 end
%         end
    end
end