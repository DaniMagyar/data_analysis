function [psth, ts, psth_1st, ts_1st] = ttl_psth_dani(spktimes,trigtimes,bins,varargin)

% Generates a peri-stimulus time histogram (psth) and returns the psth and
% the spike timestamps relative to trigger times. IMPORTANT: all timestamp
% inputs (spktimes, trigtimes) must be in the same units, ideally secinds.
%
% INPUTS 
% - spktimes: vector with timestamps of spike events (units as recorded)
% - trigtimes: vector with timestamps of trigger events (units as recorded)
% - bins: resolution of the histogram. 
% Varargin 
% - 'pre': time before trigger to include in psth. Default 1s. 
% - 'post': time after trigger to include in psth. Default 1s. 
% - 'fs': sampling frequency. Deafult 30000.
% - 'binsz': bin size of psth. Default 1 ms. 
% - 'chart': if '0' (default), no plot will be generated
%            if '1', a PSTH will be generated 
%            if '2', a PSTH together with a raster plot will be generated
%
% OUTPUT
% - psth: peri-stimulus time histogram
% - psth_1st: peri-stimulus time histogram up to the first after the trigger
% - ts: spike timestamps relative to trigger times
% - ts_1st: spike timestamps up to the first after the trigger 
%
% Examples 
% [psth, ts] = ttl_psth(timestamps, ttl, 1000);
% [psth, ts, psth_1st, ts_1st] =  ttl_psth(timestamps, ttl, 1000,'pre',0.5,...
%       'post',0.5,'chart',2);
%
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% Based on mpsth by Maik C. Stttgen, Feb 2011
% -------------------------------------------------------------------------

% Deafult params
fs = 30000;
pre   = 1 * fs;
post  = 1 * fs;
binsz = 1;
chart = 0;

% Varargin params
if nargin
  for ii=1:2:size(varargin,2)
    switch varargin{ii}
      case 'pre'
        pre = varargin{ii+1};
        if pre >= 5000
        elseif pre < 5000
            pre = pre * fs;
        end
      case 'post'
        post = varargin{ii+1};
        if post >= 5000
        elseif post < 5000
            post = post * fs;
        end
      case 'fs'
        fs = varargin{ii+1};
      case 'binsz'
        binsz = varargin{ii+1};
      case 'chart'
        chart = varargin{ii+1};
      otherwise
        fprintf(['Unknown argument:' varargin(ii,:) '\n']);
    end
  end
else
end

% pre-allocate for speed
if binsz>1
  psth = zeros(ceil(pre/binsz +post/binsz),2); % one extra chan for timebase
  psth (:,1) = (-1*pre:binsz:post-1); % time base
  psth_1st = zeros(ceil(pre/binsz+post/binsz)+1,2); 
  psth_1st (:,1) = (-1*pre:binsz:post);
elseif binsz==1
  % In this case, pre+post+1 bins are generated (ranging from pre:1:post)
  psth = zeros(pre+post+1,2);
  psth (:,1) = (-1*pre:1:post); % time base
  psth_1st = zeros(pre+post+1,2); 
  psth_1st (:,1) = (-1*pre:1:post);
end

% construct psth & trialspx
ts = cell(numel(trigtimes),1);
for ii = 1:numel(trigtimes)
  clear spikes
  spikes = spktimes - trigtimes(ii); % all spikes relative to current trigtime
  ts{ii} = round(spikes(spikes>=-pre & spikes<=post)); % spikes close to current trigtime
  if binsz==1 % just to make sure...
    psth(ts{ii}+pre+1,2) = psth(ts{ii}+pre+1,2)+1; % markers just add up
    % previous line works fine as long as not more than one spike occurs in the same ms bin
    % in the same trial - else it's omitted
  elseif binsz>1
    try
      for j = 1:numel(ts{ii})
        psth(floor(ts{ii}(j)/binsz+pre/binsz+1),2) = ...
            psth(floor(ts{ii}(j)/binsz+pre/binsz+1),2)+1;
      end
    catch
    end
  end
end

% construct psth & trialspx up to the first spike after the trigger
ts_1st = cell(numel(trigtimes),1);
for ii = 1:numel(trigtimes)
  clear spikes
  spikes = spktimes - trigtimes(ii); % all spikes relative to current trigtime
  st = find(spikes>0,1);
  st = (spikes(st));
    if isempty(st)
      st = 0;
    end
  SPpre = spikes>=-pre;
  SPpost = spikes<=post;
  SPst = spikes<=st;
  SP = SPpost & SPst;
  ts_1st{ii} = round(spikes(SPpre & SP)); % spikes close to trigtime up to 1st after trigger
  if binsz==1 % just to make sure...
     psth_1st(ts_1st{ii}+pre+1,2) = psth_1st(ts_1st{ii}+pre+1,2)+1;  
  elseif binsz>1
    try
      for j = 1:numel(ts{ii})
        psth_1st(floor(ts{ii}(j)/binsz+pre/binsz+1),2) = ...
            psth_1st(floor(ts{ii}(j)/binsz+pre/binsz+1),2)+1;
      end
    catch
    end
  end
  
end

% plot
if chart==1
  figure('name','Peri-stimulus time histogram','units','normalized','position',[0.3 0.4 0.4 0.2])
  
  [psth_spx, psth_t] = psth_hist(psth, bins);
  psth_t = psth_t + 0.5*bins;
  bar((psth_t/fs),psth_spx);
  h= get(gca,'ylim');
  axis([min(psth_1st(:,1)/fs) max(psth_1st(:,1)/fs) 0 h(2)+1]);
  hold on; 
  h= get(gca,'ylim');
  plot([0 0],[h(1) h(2)],'b');  
  %plot([-1 -1],[h(1) h(2)],'b');  %ezt írtam be
  hold off;
  ylabel('Spikes count');
  xlabel('Peri-stimulus time (s)');

elseif chart==2
  figure('name','Peri-stimulus time histogram','units','normalized','position',[0.3 0.3 0.4 0.3])
  
  subplot(223)
  [psth_spx, psth_t] = psth_hist(psth, bins);
  psth_t = psth_t + 0.5*bins;
  bar((psth_t/fs),psth_spx);
  h= get(gca,'ylim');
  axis([min(psth(:,1)/fs) max(psth(:,1)/fs) 0 h(2)+1]);
  hold on; 
  h= get(gca,'ylim');
  plot([0 0],[h(1) h(2)],'r');  
  %plot([-1 -1],[h(1) h(2)],'b');
  hold off;
  ylabel('Spikes count');
  xlabel('Peri-stimulus time (s)');
  
  subplot(221)
  rastmat = zeros(numel(ts),pre+1+post);
  timevec = -pre:1:post;  
  % generate raster
    for ii = 1:numel(ts)
        rastmat(ii,ts{ii}+pre+1) = 1;
    end
    
  for ii = 1:numel(ts)
    plot(timevec/fs,rastmat(ii,:)*ii,'Color','k','Marker','.','MarkerSize',2,'LineStyle','none');
    hold on;
  end
  hold on; h= get(gca,'ylim');
  plot([0 0],[h(1) numel(ts)+0.5],'b');  
  hold off;
  axis([min(timevec/fs) max(timevec/fs) 0.5 numel(ts)+0.5]);
  ylabel('Trials');
  
%   subplot(222)
%   rastmat = zeros(numel(ts_1st),pre+1+post);
%   timevec = -pre:1:post;  
%   % generate raster
%     for ii = 1:numel(ts_1st)
%         rastmat(ii,ts_1st{ii}+pre+1) = 1;
%     end
%     
%   for ii = 1:numel(ts_1st)
%     plot(timevec/fs,rastmat(ii,:)*ii,'Color','k','Marker','.','MarkerSize',2,'LineStyle','none');
%     hold on;
%   end
%   hold on; h= get(gca,'ylim'); 
%   plot([0 0],[h(1) numel(ts_1st)]+0.5,'r');  
%   %plot([-1 -1],[h(1) h(2)],'b');
%   hold off;
%   axis([min(timevec/fs) max(timevec/fs) 0.5 numel(ts_1st)+0.5]);
%   ylabel('Trials');
  
%   subplot(224)
%   [psth_spx, psth_t] = psth_hist(psth_1st, bins);
%   psth_t = psth_t + 0.5*bins;
%   bar((psth_t/fs),psth_spx);
%   h= get(gca,'ylim');
%   axis([min(psth_1st(:,1)/fs) max(psth_1st(:,1)/fs) 0 h(2)+1]);
%   hold on; 
%   h= get(gca,'ylim');
%   plot([0 0],[h(1) h(2)],'b');  
%   hold off;
%   ylabel('Spikes count');
%   xlabel('Peri-stimulus time (s)');
  
end
