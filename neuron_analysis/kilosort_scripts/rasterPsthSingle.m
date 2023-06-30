function rasterPsthSingle(varargin)

prs =  inputParser;
addRequired(prs,'AP',@isnumeric) 
addRequired(prs,'TTL',@isnumeric)
addOptional(prs,'window',[-1 1],@(s)isnumeric(s)&isequal(length(s),2)) 
addOptional(prs, 'bin_time', 0.02,@isnumeric)
addOptional(prs,'dt',0.001) 
addOptional(prs,'eventline',1) % Alternative 1 (yes)
parse(prs,varargin{:})
g = prs.Results;

%% Raster Plot
figure()
hold on 
time = g.window(1):g.dt:g.window(2); 
stimes       = g.AP;
event_spikes = cell(1,length(g.TTL));
ntt = 1:length(g.TTL);
for nTT = 1:length(g.TTL)
    event_spikes{nTT} = stimes(stimes>=g.TTL(nTT)+g.window(1) & ...
        stimes<g.TTL(nTT)+g.window(2)) - g.TTL(nTT);
    ind_spikes        = round((event_spikes{nTT}-time(1))/g.dt) + 1;
    if size(ind_spikes,1) == 1
        line([ind_spikes; ind_spikes],[(ntt(nTT)-0.4)*ones(1,length(ind_spikes));...
            (ntt(nTT)+0.4)*ones(1,length(ind_spikes))],'Color','k','Linewidth',0.05);
    elseif isempty(ind_spikes)
        ind_spikes = 1;
        line([ind_spikes; ind_spikes],[(ntt(nTT)-0.4)*ones(1,length(ind_spikes));...
            (ntt(nTT)+0.4)*ones(1,length(ind_spikes))],'Color','w','Linewidth',0.05);
    else
        line([ind_spikes ind_spikes]',[(ntt(nTT)-0.4)*ones(1,length(ind_spikes));...
            (ntt(nTT)+0.4)*ones(1,length(ind_spikes))],'Color','k','Linewidth',0.05);
    end
end
xlim([0 length(time)-1]);
xt = xticks;
xt2 = round(linspace(time(1),time(end),length(xt)),2);
xticklabels(xt2);
yt = yticks;
yticks(yt(end));
yticklabels(ntt(end));
ylim([0 yt(end)]);
if g.eventline == 1
    hold on; 
    X = find(time==0);
    plot([X-1 X-1],[0 yt(end)],'r','LineWidth',1);  
    hold off;
end
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1),pos(2),pos(3),150]);
xlabel('Time (s)')
ylabel('#Trials')
hold off

%% Psth Plot
bin_time = g.bin_time;   
pre_time = abs(g.window(1));      
post_time = g.window(2);      
for jj = 1:numel(g.TTL) %Each TTL is a a column. 
     preAP{:,jj} = g.AP(g.AP>=(g.TTL(jj)-pre_time) & g.AP<g.TTL(jj)); % spikes before each TTL separately
     postAP{:,jj} = g.AP(g.AP>g.TTL(jj) & g.AP<(g.TTL(jj)+post_time)); % spikes after each TTL separately
     %postAP{:,jj} = AP(AP>(TTL(jj)+0.005) & AP<(TTL(jj)+0.005+post_time)); %BA_250
 end
 for ll = 1:numel(g.TTL) %Each TTL is a a column. 
     preAP_norm{ll} = preAP{ll}-g.TTL(ll); % spikes relative to their own TTL
     postAP_norm{ll} = postAP{ll}-g.TTL(ll);
 end
 for tt = 1:numel(g.TTL) %Each TTL is a a column. 
     for nn = 1:(pre_time/bin_time) % number of timestamps in each bin.
         preAP_bin(nn,tt) = sum(preAP_norm{tt}>=(-pre_time+(nn-1)*bin_time) & preAP_norm{tt}<(-pre_time+nn*bin_time));
     end
     for oo = 1:(post_time/bin_time)
         postAP_bin(oo,tt) = sum(postAP_norm{tt}>=((oo-1)*bin_time) & postAP_norm{tt}<(oo*bin_time));
     end
 end

 psth_spikes = vertcat(sum(preAP_bin,2), sum(postAP_bin,2));
 figure(2)
 hold on
 bar((g.window(1)+bin_time:bin_time:g.window(2)),psth_spikes, 'FaceAlpha', 0.75)
 hold off

 %psth_spikes = (psth_spikes/sum(psth_spikes))*100;
 psth_spikesZsc = zscore(psth_spikes);
 psth_spikesZsc = psth_spikesZsc-(mean(psth_spikesZsc(1:pre_time/bin_time)));

 figure(3)
 hold on
 xlabel('Time (s)')
 ylabel('Z-score')
 plot((g.window(1)+bin_time:bin_time:g.window(2)),psth_spikesZsc,'LineWidth',1)
 hold off