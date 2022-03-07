function raster_plot(neuronID,event,varargin)

% Generates a raster plot for the specified unit for a specific time
% window.
%
% INPUTS 
% - neuronID: vector with timestamps of spike events (units as recorded)
% - event: times of the event (TTL) to use as 0 values.  
% Varargin 
% - window: time [before after] event to include in raster. Default [-1 1]s. 
% - 'dt': time resolution of the raster, in seconds. Default 0.01s.
% - 'folder': specify folder where the unit is located. Default pwd.
% - 'eventline': plot a line in the times of the event. Default 0 = no.
%
% OUTPUT
% - plot
% 
% Examples 
% raster_plot('GR5_3',TTL,[-3 1])
% raster_plot('GR5_3',event,'eventline',1,'window',[-1 3])
%
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

% Default arguments
prs = inputParser;
addRequired(prs,'neuronID',@(s)ischar(s)) 
addRequired(prs,'event',@(s)isnumeric(s))

addParameter(prs,'window',[-1 1],@(s)isnumeric(s)&isequal(length(s),2)) 
addParameter(prs,'dt',0.01) 
addParameter(prs,'folder',pwd)
addParameter(prs,'eventline',0) % Alternative 1 (yes)
parse(prs,neuronID,event,varargin{:})
g = prs.Results;

% Load spikes
load([neuronID '.mat'],'TS')

time = g.window(1):g.dt:g.window(2); 
stimes       = ((TS(:,1)/10000));
event_spikes = cell(1,length(g.event));
ntt = 1:length(g.event);
for nTT = 1:length(g.event)
    event_spikes{nTT} = stimes(stimes>=g.event(nTT)+g.window(1) & ...
        stimes<g.event(nTT)+g.window(2)) - g.event(nTT);
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

