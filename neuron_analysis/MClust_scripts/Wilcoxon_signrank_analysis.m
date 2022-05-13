function Wilcoxon_signrank_analysis(File, nucleus, Stim)
% Default params ----------------------------------------------------------
int      = [-2 2]; % psth interval, even numbers are better for plot
mainFolder = 'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december';
% Select neurons assigned to the selected NUCLEUS and Set the folders -----
Tab = File;
for ii = size(File,1):-1:1
    if startsWith(File.Location(ii), nucleus) == 0
        Tab(ii,:) = [];
    end
end


folder     = Tab.Recording;
nNeurons   = size(folder,1);

for kk = 1:nNeurons
    % Go to folder
    cd([mainFolder filesep folder{kk,:}]);
    
    % Load TTLs & Choose variable to use as stimulus
    load('TTLs.mat'); %#ok<LOAD>
    switch Stim
        case 'BA_500'
            ttl = BA_500;
        case 'BA_500_5Hz'
            ttl = BA_500_5Hz;
        case 'BA_50'
            ttl = BA_50;
        case 'TTL_500'
            ttl = TTL;
        case 'TTL_50'
            ttl = TTL(1:10:end);
        case 'BA_25'
            ttl = BA_25;
        case 'BA_25_5Hz'
            ttl = BA_25_5Hz;
%            ttl = BA_25_5Hz(11:25);
        case 'BA_25_10Hz'
            ttl = BA_25_10Hz;
        case 'TO_25'
            ttl = TO_25;
        case 'TO_25_5Hz'
            ttl = TO_25_5Hz;
%            ttl = TO_25_5Hz(1:10);
        case 'TO_25_10Hz'
            ttl = TO_25_10Hz;
        case 'BA_250'
            ttl = BA_250;
        case 'BA_250_5Hz'
            ttl = BA_250_5Hz;
%           ttl = BA_250_5Hz(1:125);
%             ttl = [ttl; BA_250_5Hz(7:10:end)];
%             ttl = [ttl; BA_250_5Hz(8:10:end)];
%             ttl = [ttl; BA_250_5Hz(9:10:end)];
%             ttl = [ttl; BA_250_5Hz(10:10:end)];
%             ttl = sort(ttl);
        case 'BA_250_10Hz'
            ttl = BA_250_10Hz;
        case 'TO_250'
            ttl = TO_250;
        case 'TO_250_5Hz'
            ttl = TO_250_5Hz;
%            ttl = TO_250_5Hz(1:125);
        case 'TO_250_10Hz'
            ttl = TO_250_10Hz;
        case 'BA_500_5Hz_10Hz'
            ttl = BA_250_5Hz;
            ttl = [ttl; BA_250_10Hz];
            ttl = sort(ttl);
        case 'SK'
            ttl = TTL; % The SK TTL is already called TTL.
    end
    % Select the neuron
    NeuronID = dir(['GR',Tab.Group{kk},'_',num2str(Tab.Neuron(kk)),...
        '.mat']);
    NeuronID = NeuronID.name;
    
    load(NeuronID,'TS')
    TT = TS(:,1)/10000;

    %% Generating preAP and postAP matrices
    bin_time = int(2);     
    pre_time = abs(int(1));      
    post_time = int(2);      
    AP = TT;
    TTL=ttl;
    for jj = 1:numel(TTL)
         preAP{:,jj} = AP(AP>=(TTL(jj)-pre_time) & AP<TTL(jj));
         postAP{:,jj} = AP(AP>TTL(jj) & AP<(TTL(jj)+post_time));
     end
     for ll = 1:numel(TTL)
         preAP_norm{ll} = preAP{ll}-TTL(ll);
         postAP_norm{ll} = postAP{ll}-TTL(ll);
     end
     for mm = 1:numel(TTL)
         for nn = 1:(pre_time/bin_time)
             preAP_bin(nn,mm) = sum(preAP_norm{mm}>=(-pre_time+(nn-1)*bin_time) & preAP_norm{mm}<(-pre_time+nn*bin_time));
         end
         for oo = 1:(post_time/bin_time)
             postAP_bin(oo,mm) = sum(postAP_norm{mm}>=((oo-1)*bin_time) & postAP_norm{mm}<(oo*bin_time));
         end
     end
     preAP_bin_freq = preAP_bin/bin_time;
     postAP_bin_freq = postAP_bin/bin_time;
     %% paired Wilcoxon signed rank test: equal length of pre and post data. 
     % Each post point cointain all spikes during the whole illumination period. (Like one 2 secs long bin)
    
     preAP_bin_freq_vector = reshape(preAP_bin_freq,[],1);
     postAP_bin_freq_vector = reshape(postAP_bin_freq,[],1);
     bin_freq_Origin = horzcat(preAP_bin_freq_vector, postAP_bin_freq_vector);
     [~,h] = signrank(preAP_bin_freq_vector, postAP_bin_freq_vector);
     Wilcoxon_results{kk,1} = NeuronID;
     Wilcoxon_results{kk,2} = num2str(h);
end

Wilcox_values = str2num(cell2mat(Wilcoxon_results(:,2)));
significant_idx = find(Wilcox_values == 1);
non_significant_idx = find(Wilcox_values == 0);
figure 2
plot(ones(1,length(significant_idx)), significant_idx, 'r.', 'MarkerSize', 20)
hold on
plot(ones(1,length(non_significant_idx)), non_significant_idx, 'b.', 'MarkerSize', 20)















