function [psth_spikes] = calculate_psth(AP, TTL, preferences)
% Input:
% -AP: in seconds, from recording start point
% -TTL: in seconds, from recording start point
% -int: in seconds, pre and post TTL e.g.:[-2 2]
% -bin: in seconds

bin_time = preferences.psth_bin/preferences.fs;   
pre_time = abs(preferences.int(1));      
post_time = preferences.int(2);      
for jj = 1:numel(TTL) %Each TTL is a a column. 
     preAP{:,jj} = AP(AP>=(TTL(jj)-pre_time) & AP<TTL(jj)); % spikes before each TTL separately
     postAP{:,jj} = AP(AP>TTL(jj) & AP<(TTL(jj)+post_time)); % spikes after each TTL separately
     %postAP{:,jj} = AP(AP>(TTL(jj)+0.005) & AP<(TTL(jj)+0.005+post_time)); %BA_250
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

 psth_spikes = vertcat(sum(preAP_bin,2), sum(postAP_bin,2));
 %psth_spikes(abs(preferences.int(1)*preferences.fs/preferences.psth_bin)+1)=psth_spikes(abs(preferences.int(1)*preferences.fs/preferences.psth_bin));%mean(sum(preAP_bin,2)); %change zero bin to pre_mean


%  % Collecting normalized timestamps into vector
%  all_postAP_norm =[];
% for mm = 1:numel(TTL)
%     all_postAP_norm = vertcat(all_postAP_norm, postAP_norm{mm}); 
% end
% disp('end')