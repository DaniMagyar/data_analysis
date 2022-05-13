function plot_PsthPCA(PSTHall, leafOrder, Dend, Clusters, preferences)

% PCA Figure --------------------------------------------------------------
figure; 
subplot(1,4,1:2)
imagesc(preferences.time/preferences.fs, 1:size(PSTHall,1), PSTHall(leafOrder,:)); %plotting sorted psth matrix 
clim([-max(max(PSTHall,[],1))/2.5 max(max(PSTHall,[],1))]);
colormap(preferences.mycolormap); 
hold on;
plot([-0.00125 -0.00125],[1 size(PSTHall,1)],'r');
hold off;
ylabel('# Cell');
xlabel('Time (s)')

subplot(1,4,3:4)
cutoff = Dend(end-preferences.clustnum+2,3); % cutting the tree to get 'clustnum' clusters
h = dendrogram(Dend,0,'Orientation','right','reorder',leafOrder,...
    'ColorThreshold',cutoff); %plotting dendrogram
set(h,'LineWidth',1)
set(gca,'Ydir','reverse');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'visible','off')


% The first cluster is the 'deepest'
if preferences.average == 1
    for ll = 1:preferences.clustnum
        Clustmeans(ll,:) = mean(PSTHall((find(Clusters==ll)),:),1);
    end
% 
%     averageFig = tiledlayout (clustnum,1);
%     for mm = 1:clustnum
%         nexttile
%         plot(time/fs, Clustmeans(manualclusterorder(mm),:))
%     end
    for mm = 1:preferences.clustnum
        figure(mm+1)
        hold on
        plot(preferences.time/preferences.fs, Clustmeans(mm,:), 'LineWidth', 1.5)
        title(['n neurons = ' num2str(numel(find(Clusters==mm)))])
        hold off
    end
end