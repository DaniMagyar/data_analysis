function ccg_HB(cellids,window,varargin)
%CCG   Cross-correlation.
%   CCG(CELLIDS,WIN) calculates cross-correlations at 1 ms resolution for a
%   given window size (WIN). For details on the algorithm, see
%   SOMCCG_CONF_FILTER. Minimal shift for shuffled CCGs is set to 1100 ms.
%
%   Required input arguments:
%       CELLIDS: pairs of cells within CELLIDS are subjected to CCG.
%           Selection is determined by 'whichcells' and 'include'
%           properties (see below).
%       WINDOW: CCG window size in seconds.
%
%   Optional input parameter-value pairs:
%       'issave', false - controls saving behavior; plots and
%           cross-correlation matrices with confidence intervals are saved 
%           only if 'issave' is set to true
%       'resdir' - results directory. A 'CCG' directory is created and used
%           on the analysis path if 'issave' is set to 'true' but no
%           results directory is specified
%       'whichcells', 'nontetrodepairs' - method of pair selection;
%           'nontetrodepairs' selects cells from other tetrodes,
%           'tetrodepairs' selects cells from the same tetrode and
%           'allpairs' selects all cells from the session
%       'include', 'list' - by default, only pairs for which both cells are
%           included in CELLIDS are analyzed; if 'include' is set to 
%           'cellbase', all cells in CellBase that are paired with the ones
%           in CELLIDS according to 'whichcells' are analyzed
%       'minspikeno', 100 - minimal spike number to perform CCG calculation
%       'maxspikeno', 5000 - maximal spike number; the first N spikes are
%           included in the CCG calculation
%       'segfilter', 'none' - filter recording segments (see FINDSEGS3)
%           If set to 'none' (default), the entire recording is used
%       'filterinput', [] - additional input for FINDSEGS3
%       'longsegments', false - if true, only the longest segment is
%           analyzed after filtering
%       'seglim', 0.3 - if 'longsegments' is 'true', CCG is calculated if
%           the longest segment reaches 'seglim' (in seconds)
%
%   Examples for possible segment filters:
%         tseg = findSegs3(cell1,'segfilter','stimfb_excl_nb',...
%             'light_activation_duration',[-5 5],'margins',[0 0]);
%         tseg = findSegs3(cell1,'segfilter','prestim3');
%         tseg = findSegs3(cell1,'segfilter','fb_incl_nb',...
%             'feedback_duration',[-0.5 0.5],'margins',[0 0],'min_int',0);
%         tseg = findSegs3(cell1,'segfilter','cue_incl_nb',...
%             'feedback_duration',[-1.5 0],'margins',[0 0],'min_int',0);
%
%   See also ACG, SOMCCG_CONF_FILTER and XCORR_WRAND_FILTER.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   14-June-2013

% Input arguments
prs = inputParser;
addRequired(prs,'cellids',@(s)iscell(s)|iscellstr(s)|ischar(s))
addRequired(prs,'window',@isscalar)  % CCG window, in seconds
addParameter(prs,'issave',false,@islogical)   % control saving behavior
addParameter(prs,'resdir',[],@(s)isdir(s)|isempty(s))   % results directory
addParameter(prs,'longsegments',false,@islogical)   % use only the longest segment after segment filtering
addParameter(prs,'seglim',0.3,@isnumeric);   % minimal segment length (s) to perform CCG calculation if 'longsegments' is 'true'
addParameter(prs,'segfilter','none',@(s)ischar(s)|iscellstr(s))   % filter segments
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'minspikeno',100,@isnumeric)   % calculate CCG above minimal spike number
addParameter(prs,'maxspikeno',5000,@isnumeric)   % maximize included spikes
addParameter(prs,'whichcells','nontetrodepairs',...
    @(s)ismember(s,{'nontetrodepairs','tetrodepairs','allpairs'}))   % which cells to include
addParameter(prs,'include','list',...
    @(s)ismember(s,{'list','cellbase'}))   % cell pair selection behavior: only from input list or full CellBase
parse(prs,cellids,window,varargin{:})
g = prs.Results;

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
if isempty(g.resdir)
    resdir = fullfile(DATAPATH,'CCG');  % results directory
    if ~isdir(resdir)
        mkdir(resdir)
    end
end
fnmm = 'CCG_matrices.mat';   % filename for saving the result matrices

% Cell pairs
PairOfCells = cell(0,2);
numCells = length(cellids);  % number of cells
for iC = 1:numCells
    cellid = cellids{iC};
    [animalID, sessionID, tetrode1, unit1] = cellid2tags(cellid);
    switch g.whichcells
        case 'nontetrodepairs'
            [nm, ps] = nontetrodepairs(cellid);   % cells on other tetrodes
        case 'tetrodepairs'
            [~, ps] = tetrodepairs(cellid);   % cells on the same tetrode
            ps = setdiff(ps,cellid);
            nm = length(ps);
        case 'allpairs'
            ps = findcell('rat',animalID,'session',sessionID);   % all concurrently recorded cells
            ps = setdiff(ps,cellid);
            nm = length(ps);
    end
    if isequal(g.include,'list')
        ps = intersect(cellids,ps);   % include only those that are in the input list
        nm = length(ps);
    end
    for k = 1:nm
        [~, ~, tetrode2, unit2] = cellid2tags(ps(k));
        if (tetrode2*10+unit2) > (tetrode1*10+unit1)   % prevent duplicates
            PairOfCells(end+1,1:2) = {cellid ps{k}}; %#ok<AGROW>
        end
    end
end
numPairs = size(PairOfCells,1);

% Determine time window
sr = 1000;      % sampling rate
wn = g.window * sr;    % CCG window in data points

% Cell pair loop for CCG
wb = waitbar(0,'Please wait...','Name','Running CCG...');  % progress indicator
global WB
WB(end+1) = wb;
limit_spikes = [g.minspikeno g.maxspikeno];   % include max 50000 spikes; calculate only if min 100 spikes by default
[CCR, LCCR, UCCR] = deal(zeros(numPairs,2*wn+1));
[MeanH0, SDH0] = deal(nan(numPairs,1));
SegmentLength = nan(numPairs,1);
for iP = 1:numPairs   % loop through pairs of cells
    cell1 = PairOfCells{iP,1};
    cell2 = PairOfCells{iP,2};
    if isequal(g.segfilter,'none')
        ncc1 = loadcb(cell1,'SPIKES');   % use all spikes
        ncc2 = loadcb(cell2,'SPIKES');
    else
        tseg = findSegs3(cell1,'segfilter',g.segfilter,...
            g.filterinput{:});  % find time segments
        ltseg = tseg(2,:) - tseg(1,:);  % length of the segments
        if g.longsegments   % use the longest segment if it's longer than the threshold ('seglim')
            seginx = find(ltseg==max(ltseg));
            tseg = tseg(:,seginx(1));   % find the longest segment
            ltseg = ltseg(seginx(1));
            if tseg < g.seglim * sr
                continue
            end
        end
        SegmentLength(iP) = sum(ltseg);  % cumulative length of the segments
        ncc1 = extractSegSpikes(cell1,tseg);   % find spikes in the time segments
        ncc2 = extractSegSpikes(cell2,tseg);
    end
    
    % Implement upper spike number limits
    if length(ncc1) > limit_spikes(2)     % crop if too long to avoid out of memory
        ncc1 = ncc1(1:limit_spikes(2));
    end
    if length(ncc2) > limit_spikes(2)
        ncc2 = ncc2(1:limit_spikes(2));
    end
    
    if length(ncc1) > limit_spikes(1) && length(ncc2) > limit_spikes(1)     % implement minimum number of spikes
        [H1, ccr, lwr, upr, rccg] = somccg_conf_filter(ncc1,ncc2,wn,1100);    % 1->2
        xl = xlim;   % put the cell IDs in the plot
        yl = ylim;
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.9,regexprep(cell1,'_',' '))
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.85,regexprep(cell2,'_',' '))
        if g.issave   % save figure
            ncl1 = regexprep(cell1,'\.','_');
            ncl2 = regexprep(cell2,'\.','_');
            fnm = ['CCG_' ncl1 '_' ncl2 '.fig'];
            saveas(H1,fullfile(g.resdir,fnm))   % save CCG plot
            close(H1)
        end
        CCR(iP,:) = ccr;   % cross-correlogram
        LCCR(iP,:) = lwr;  % lower significance limit
        UCCR(iP,:) = upr;  % upper significance limit
        MeanH0(iP) = mean(mean(rccg,2),1);   % surrogate mean
        SDH0(iP) = mean(std(rccg,[],2),1);   % surrogate SD
        disp(['Pair #' num2str(iP) ' / ' num2str(numPairs) ' done......'])
    end
    
    % Save
    if g.issave
        if isequal(mod(iP,50),0)   % save after every 50 pairs to prevent data loss
            save(fullfile(g.resdir,fnm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength')
            disp('Autosave done.')
        end
    end
    
    waitbar(iP/numPairs)
end
close(wb)   % eliminate progress indicator

% Save
if g.issave
    save(fullfile(g.resdir,fnmm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength')
end