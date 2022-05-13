function [h] = calculate_Wilcoxon(AP, TTL, preferences)

 %% paired Wilcoxon signed rank test: equal length of pre and post data. 
 % Each post point cointain all spikes during the whole illumination period. (Like one 2 secs long bin)
 for wcx = 1:numel(TTL) %Each TTL is a a column. 
     preAP_wcx{:,wcx} = AP(AP>=(TTL(wcx)-abs(preferences.Wcx_win(1))) & AP<TTL(wcx)); % spikes before each TTL separately
     postAP_wcx{:,wcx} = AP(AP>TTL(wcx) & AP<(TTL(wcx)+preferences.Wcx_win(2))); % spikes after each TTL separately
     %postAP_wcx{:,wcx} = AP(AP>(TTL(wcx)+0.005) & AP<(TTL(wcx)+0.005+preferences.Wcx_win(2))); %BA_250
 end
 preAP_wcx_num = cellfun(@numel, preAP_wcx);
 postAP_wcx_num = cellfun(@numel, postAP_wcx);
 preAP_wcx_freq = preAP_wcx_num/abs(preferences.Wcx_win(1));
 postAP_wcx_freq = postAP_wcx_num/preferences.Wcx_win(2);
 [~,h] = signrank(preAP_wcx_freq, postAP_wcx_freq);
 h = num2str(h);