function [StimResponses, FirstSpxEst, SmoothedZscore, psth_spx_dani] =  BAfc_find_sound_v03(inputTTL)

cell_metrics = BAfc_load_neurons;

for ii = 1:numel(cell_metrics.cellID)
    disp(ii)
    TTL = cell_metrics.general.(inputTTL){ii};
    AP = cell_metrics.spikes.times{ii};
    pre_time = 10; % in sec     
    post_time = 4.5;    
    bin_time = 0.001; %in seconds
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
    psth_spx_dani(ii,:) = vertcat(sum(preAP_bin,2), sum(postAP_bin,2));

    sigma = 20; % standard deviation of the Gaussian kernel (in ms)
    binWidth = 1; % bin width for the spike train (in ms)
    edges = -3*sigma:binWidth:3*sigma; % create kernel edges
    gaussKernel = exp(-edges.^2 / (2*sigma^2)); % Gaussian formula
    gaussKernel = gaussKernel / sum(gaussKernel); % normalize kernel
    smoothed_spike_train = conv(psth_spx_dani(ii,:), gaussKernel, 'same'); % convolve
    smoothed_spike_train_zscore(ii,:) = (smoothed_spike_train - mean(smoothed_spike_train(1:10000)))/std(smoothed_spike_train(1:10000));
    clearvars -except smoothed_spike_train_zscore cell_metrics ii psth_spx_dani inputTTL
end

% Find responses in 5 sec time window
[idx_exc,~] = find(smoothed_spike_train_zscore(:,10001:14500)>=5);
idx_exc = sort(idx_exc);
repeats_exc = diff(idx_exc) == 0; 
idx_exc(repeats_exc) = [];
% Find responses in 200ms time window
window200 = horzcat(...
    smoothed_spike_train_zscore(:,10001:10200), ...
    smoothed_spike_train_zscore(:,11001:11200), ...
    smoothed_spike_train_zscore(:,12001:12200),...
    smoothed_spike_train_zscore(:,13001:13200),...
    smoothed_spike_train_zscore(:,14001:14200));

[idx_exc200,~] = find(window200>=5);
idx_exc200 = sort(idx_exc200);
repeats_exc = diff(idx_exc200) == 0; 
idx_exc200(repeats_exc) = [];
% Find lowFR neurons in 200ms time window
window200_spx = horzcat(...
    psth_spx_dani(:,10001:10200), ...
    psth_spx_dani(:,11001:11200), ...
    psth_spx_dani(:,12001:12200),...
    psth_spx_dani(:,13001:13200),...
    psth_spx_dani(:,14001:14200));

window200_spx = sum(window200_spx,2);
[idx_lowFR,~] = find(window200_spx<10);
idx_lowFR = sort(idx_lowFR);
repeats_lowFR = diff(idx_lowFR) == 0; 
idx_lowFR(repeats_lowFR) = [];

% Find inhibitory responses 200ms time window
[idx_inh200,~] = find(window200<=-3);
idx_inh200 = sort(idx_inh200);
repeats_inh = diff(idx_inh200) == 0; 
idx_inh200(repeats_inh) = [];

% smaller inhibitory responses
window1mean = mean(smoothed_spike_train_zscore(:,10001:10200),2);
window2mean = mean(smoothed_spike_train_zscore(:,11001:11200),2);
window3mean = mean(smoothed_spike_train_zscore(:,12001:12200),2);
window4mean = mean(smoothed_spike_train_zscore(:,13001:13200),2);
window5mean = mean(smoothed_spike_train_zscore(:,14001:14200),2);

meanmatrix_inh(:,1) = window1mean<=-1;
meanmatrix_inh(:,2) = window2mean<=-1;
meanmatrix_inh(:,3) = window3mean<=-1;
meanmatrix_inh(:,4) = window4mean<=-1;
meanmatrix_inh(:,5) = window5mean<=-1;

meanmatrix_inhSum = sum(meanmatrix_inh,2);
idx_inhZ1 = find(meanmatrix_inhSum>=3);

%smaller excitatory responses
meanmatrix_exc(:,1) = window1mean>=3;
meanmatrix_exc(:,2) = window2mean>=3;
meanmatrix_exc(:,3) = window3mean>=3;
meanmatrix_exc(:,4) = window4mean>=3;
meanmatrix_exc(:,5) = window5mean>=3;

meanmatrix_excSum = sum(meanmatrix_exc,2);
idx_excZ3 = find(meanmatrix_excSum>=3);

idx_biphasic200 = intersect(idx_exc200, idx_inh200);
idx_exc200filt = setdiff(idx_exc200, idx_lowFR);

% Firstspike calculation
% Smoothdata
psth_spx_sound = ...
    psth_spx_dani(:,10001:10200) + ...
    psth_spx_dani(:,11001:11200) + ...
    psth_spx_dani(:,12001:12200) + ...
    psth_spx_dani(:,13001:13200) + ...
    psth_spx_dani(:,14001:14200);

for ii = 1:size(psth_spx_sound,1)
    smoothed_spx(ii,:) = smoothdata(psth_spx_sound(ii,:),'movmean', 10); 
    PeakBin = find(smoothed_spx(ii,:) == max(smoothed_spx(ii,:)), 1, 'first');
    if isempty(PeakBin)
        PeakSpike(ii,1) = 0;
        Estimate(ii,1) = 0;
    else
        PeakSpike(ii,1) = PeakBin;
        Estimate(ii,1) = find(smoothed_spx(ii,:) >= mean([max(smoothed_spx(ii,:)) min(smoothed_spx(ii,:))]), 1, 'first');
    end
end

psth_spx_sound_plotra = ...
    psth_spx_dani(:, 9751:10250) + ...
    psth_spx_dani(:,10751:11250) + ...
    psth_spx_dani(:,11751:12250) + ...
    psth_spx_dani(:,12751:13250) + ...
    psth_spx_dani(:,13751:14250);
cell_metrics.general.psth_spx_dani = psth_spx_sound_plotra;


% % Detailed classification
% cell_metrics.FirstSpikeEst = Estimate';
% cell_metrics.labels(1:end) = {'neutral'};
% cell_metrics.labels(idx_inhZ1) = {'inhZ1'};
% cell_metrics.labels(idx_excZ3) = {'excZ3'};
% cell_metrics.labels(idx_exc200) = {'exc200'};
% cell_metrics.labels(idx_inh200) = {'inh200'};
% cell_metrics.labels(intersect(idx_exc200,idx_lowFR)) = {'exc200_lowFR'};
% if numel(idx_biphasic200) >= 1
%     error('idx_biphasic200 present, but not solved')
% end


% Simple classification
cell_metrics.FirstSpikeEst = Estimate';
cell_metrics.labels(1:end) = {'neutral'};
cell_metrics.labels(idx_inhZ1) = {'neutral'};
cell_metrics.labels(idx_excZ3) = {'neutral'};
cell_metrics.labels(idx_exc200) = {'exc'};
cell_metrics.labels(idx_inh200) = {'inh'};
cell_metrics.labels(intersect(idx_exc200,idx_lowFR)) = {'neutral'};
if numel(idx_biphasic200) >= 1
    error('idx_biphasic200 present, but not solved')
end

FirstSpxEst = cell_metrics.FirstSpikeEst;
StimResponses = cell_metrics.labels;
SmoothedZscore = smoothed_spike_train_zscore;




% Ez masolhato majd a main scriptbe
bR_order_flip = {'LA','BA_med','BA_lat','CeA','Astria'};
bR_order_flipnames = {'LA','BA Med','BA Lat','CeA','Astria'};
figure(4)
t = tiledlayout(1,5);
for ii = 1:numel(bR_order_flip)
    nexttile
    histogram(Estimate(intersect(find(Estimate>3),intersect(idx_exc200filt, find(contains(cell_metrics.brainRegion,bR_order_flip{ii}))))),'BinWidth',2) % remove <3ms responses
    xlim([0 100])
    xlabel('Time(ms)')
    ylabel('Number of neurons')
    set(gca,'XTick',0:20:100);
    set(gca,'TickLength',[0.02, 0.02])
    set(gca, 'TickDir', 'out')
    set(gca,'box','off')
    title(bR_order_flipnames{ii})
end
set(findall(gcf,'-property','FontSize'),'FontSize',15)
title(t,'Spike latencies after CS delivery', 'Fontsize', 20)




