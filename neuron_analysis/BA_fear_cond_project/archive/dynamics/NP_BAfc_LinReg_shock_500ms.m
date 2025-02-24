function [p] = NP_BAfc_LinReg_shock_500ms


cell_metrics = NP_BAfc_load_neurons;

for zz = 1:2
    if zz == 1
        inputTTL = 'TTL_shocks_predicted';
        session = 'predicted';
    elseif zz == 2
        inputTTL = 'TTL_shocks_nonpredicted';
        session = 'NONpredict';
    end 

    for ii = 1:numel(cell_metrics.cellID)
        disp(ii)
        TTL = cell_metrics.general.(inputTTL){ii};
        AP = cell_metrics.spikes.times{ii};
        pre_time = 0.25; % in sec     
        post_time = 15;    
        bin_time = 0.001; %in seconds
        timewin = 500; %in ms
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
    
        preAP_bin = transpose(preAP_bin);
        postAP_bin = transpose(postAP_bin);
    
        psth_spx = [preAP_bin postAP_bin];
        psth_spx_zscore = zscore(psth_spx, 0, 'all');
        psth_spx_zscore_mean = smoothdata(mean(psth_spx_zscore,1));
        resp_mean_zscore(ii) = mean(psth_spx_zscore_mean(pre_time/bin_time+1:pre_time/bin_time+timewin));
    
        for kk = 1:size(psth_spx_zscore,1)
            aucShock_2sec(kk) = trapz(1:timewin, psth_spx_zscore(kk,pre_time/bin_time+1:pre_time/bin_time+timewin));
        end
        x = 1:20;
        y = aucShock_2sec;
        % b = bar(x, y);
        % b.FaceColor = "#0072BD";
        % hold on
        % Perform linear regression (1st degree polynomial fitting)
        p(ii,1:2) = polyfit(x, y, 1);  % p contains the slope and intercept of the line
        
        % Generate y-values for the regression line
        y_fit = polyval(p(ii,1:2), x);
    
        clear preAP_bin postAP_bin preAP postAP preAP_norm postAP_norm % reset matrices, because different TTL num    
    end
    cell_metrics.(['Resp_mean_zscore' session]) = resp_mean_zscore;
    cell_metrics.(['SlopeShock' session]) = p(:,1)';
    cell_metrics.(['InterceptShock' session]) = p(:,2)';

end

[~, idx] =  BAfc_find_shocks_v05('NP_BAfc', 500);
cell_metrics.labels(1:end) = {'0'};
cell_metrics.labels(idx) = {'idx'};

respmean_diff = NP_BAfc_response_mean_diff;
cell_metrics.Respmean_diff = respmean_diff';

cell_metrics = CellExplorer('metrics',cell_metrics);