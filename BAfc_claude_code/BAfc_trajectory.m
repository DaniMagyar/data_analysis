%% ========================================================================
% POPULATION TRAJECTORIES BY REGION (BA, LA, Astria, CeA)
% ========================================================================
% This standalone script plots ONLY the neural population trajectories in PC
% space (PC1–PC3) for four categories:
%   1) LA–PN, 2) LA–IN, 3) BA–PN, 4) BA–IN
% using PSTHs computed from Tone, Shock, and Combined trials.
%
% It also prints a quantitative/statistical interpretation to the console:
% - Variance explained by PC1–PC3
% - Trajectory length in the response window
% - Baseline→response displacement
% - Pairwise separations (mean Euclidean) with bootstrap 95% CIs
% - Mahalanobis distances between conditions
% - Angles between displacement vectors
%
% REQUIREMENTS:
% - A data struct `g` in workspace or loaded from a .mat file with fields:
%     g.cell_metrics.cellID
%     g.cell_metrics.putativeCellType{i}  -> 'PN' / 'IN'
%     g.cell_metrics.brainRegion{i}       -> 'LA' / 'BA'
%     g.cell_metrics.spikes.times{i}      -> spike times (s)
%     g.cell_metrics.general.triptest_sound_only{i} -> tone times (s)
%     g.cell_metrics.general.triptest_shocks_only{i} -> shock times (s)
%     g.cell_metrics.general.triptest_both{i} -> combined times (s)
%
% If needed, uncomment the load line and point to your data file.
% ========================================================================

clear; clc; close all;
% recordings = {...
%     'MD292_002_kilosort',...
%     'MD293_001_kilosort',...
%     'MD294_001_kilosort',...
%     'MD295_001_kilosort',...
%     'MD296_001_kilosort',...
%     'MD297_001_kilosort',...
%     'MD298_001_kilosort',...
%     'MD299_001_kilosort',...
%     'MD300_001_kilosort',...
%     'MD304_001_kilosort',...
%     'MD307_001_kilosort',...
%     'MD309_001_kilosort',...
%     'MD311_002_kilosort',...
%     'MD312_001_kilosort',...
%     'MD313_001_kilosort',...
%     'MD314_001_kilosort',...
%     'MD315_001_kilosort',...
%     'MD316_002_kilosort',...
%     'MD317_001_kilosort',...
%     'MD318_001_kilosort',...
%     'MD318_002_kilosort',...
%     'MD319_003_kilosort'};

recordings = {...
   'MD298_001_kilosort',...
   'MD299_001_kilosort',...
   'MD300_001_kilosort',...
   'MD304_001_kilosort',...
   'MD307_001_kilosort',...
   'MD309_001_kilosort',...
   'MD315_001_kilosort',...
   'MD316_002_kilosort',...
   'MD317_001_kilosort',...
   'MD318_001_kilosort',...
   'MD318_002_kilosort',...
   'MD319_003_kilosort'};

%ttl = {'triptest_sound_only','triptest_shocks_only', 'triptest_both'};

ttl = {'triptest_sound_only','triptest_shocks_only', 'triptest_both'};


hmptitles = {'CS', 'US', 'CS + US'}; % Titles for plots

g.cell_metrics = BAfc_load_neurons('recordings', recordings, 'ttl', ttl);

assert(exist('g','var')==1, 'This script expects a struct `g` in the workspace.Please load it (uncomment the load line) before running.');

%% PARAMETERS
params.pre_stim        = 1.0;      % seconds before stimulus
params.post_stim       = 2.0;      % seconds after stimulus
params.bin_size        = 0.01;     % sec (10 ms)
params.smooth_sigma    = 3;        % Gaussian smoothing sigma (in bins)
params.min_trials      = 5;        % min trials per condition to include a cell
params.min_cells_pca   = 3;        % min cells per category to run PCA
params.baseline_window = [-1.0, -0.1]; % baseline window (s)
params.response_window = [0, 0.5];     % response window (s)

% Time vectors
params.time_vec  = -params.pre_stim:params.bin_size:params.post_stim; % edges
params.n_bins    = numel(params.time_vec) - 1;
time_centers     = params.time_vec(1:end-1) + params.bin_size/2;     % bin centers
[~, zero_idx]    = min(abs(time_centers));                             % closest index to 0 s

%% EXTRACT & PREPARE
n_cells = numel(g.cell_metrics.cellID);
cellType    = cell(n_cells,1);
brainRegion = cell(n_cells,1);

% PSTH matrices (neurons x time)
psth_tone  = zeros(n_cells, params.n_bins);
psth_shock = zeros(n_cells, params.n_bins);
psth_both  = zeros(n_cells, params.n_bins);

% Trial counts
n_tone  = zeros(n_cells,1);
n_shock = zeros(n_cells,1);
n_both  = zeros(n_cells,1);

for i = 1:n_cells
    cellType{i}    = g.cell_metrics.putativeCellType{i};
    brainRegion{i} = g.cell_metrics.brainRegion{i};
    spike_times    = g.cell_metrics.spikes.times{i};

    tone_times  = g.cell_metrics.general.triptest_sound_only{i};
    shock_times = g.cell_metrics.general.triptest_shocks_only{i};
    both_times  = g.cell_metrics.general.triptest_both{i};

    % tone_times  = g.cell_metrics.general.triptest_sound_only_light{i};
    % shock_times = g.cell_metrics.general.triptest_shocks_only_light{i};
    % both_times  = g.cell_metrics.general.triptest_both_light{i};


    n_tone(i)  = numel(tone_times);
    n_shock(i) = numel(shock_times);
    n_both(i)  = numel(both_times);

    % Compute PSTHs (FR in Hz) for each condition
    psth_tone(i,:)  = compute_psth(spike_times, tone_times,  params);
    psth_shock(i,:) = compute_psth(spike_times, shock_times, params);
    psth_both(i,:)  = compute_psth(spike_times, both_times,  params);
end

% Cell-quality filter: sufficient trials in all 3 conditions
valid_cells = (n_tone >= params.min_trials) & (n_shock >= params.min_trials) & (n_both >= params.min_trials);
if ~any(valid_cells)
    warning('No cells meet the min trial criterion (min_trials=%d). Nothing to plot.', params.min_trials);
    return;
end

%% DEFINE CATEGORIES (Regions: BA, LA, Astria, CeA)
cat(1).name = 'LA';      cat(1).mask = valid_cells & strcmpi(brainRegion,'LA');
cat(2).name = 'BA';      cat(2).mask = valid_cells & strcmpi(brainRegion,'BA');
cat(3).name = 'Astria';  cat(3).mask = valid_cells & strcmpi(brainRegion,'Astria');
cat(4).name = 'CeA';     cat(4).mask = valid_cells & strcmpi(brainRegion,'CeA');

n_cat = numel(cat);

%% PLOT 3D TRAJECTORIES (PC1–PC3) PER CATEGORY
figure('Position',[80 80 1350 850]);
nb = params.n_bins;  % bins per condition
nrows = ceil(n_cat/2); ncols = min(2,n_cat);
for k = 1:n_cat
    idx = cat(k).mask; n_k = sum(idx);
    subplot(nrows, ncols, k);

    if n_k < params.min_cells_pca
        text(0.05, 0.5, sprintf('%s: not enough cells (n=%d)', cat(k).name, n_k), 'FontSize', 12);
        axis off; title(cat(k).name); continue;
    end

    % Concatenate time across the three conditions for this category
    % Rows = neurons in cat; Cols = time (tone | shock | both)
    psth_concat = [psth_tone(idx,:), psth_shock(idx,:), psth_both(idx,:)];

    % PCA across neurons (observations=time bins; variables=neurons)
    [~, score, ~, ~, explained] = pca(psth_concat');

    % Indices for each condition block
    t1 = 1:nb; t2 = nb + (1:nb); t3 = 2*nb + (1:nb);

    % 3D trajectories
    hold on;
    h1 = plot3(score(t1,1), score(t1,2), score(t1,3), '.-','LineWidth',1.5,'MarkerSize',10);
    h2 = plot3(score(t2,1), score(t2,2), score(t2,3), '.-','LineWidth',1.5,'MarkerSize',10);
    h3 = plot3(score(t3,1), score(t3,2), score(t3,3), '.-','LineWidth',1.5,'MarkerSize',10);

    % Mark t = 0 for each trajectory
    plot3(score(t1(zero_idx),1), score(t1(zero_idx),2), score(t1(zero_idx),3), 'o', 'MarkerFaceColor', get(h1,'Color'), 'MarkerEdgeColor','none');
    plot3(score(t2(zero_idx),1), score(t2(zero_idx),2), score(t2(zero_idx),3), 'o', 'MarkerFaceColor', get(h2,'Color'), 'MarkerEdgeColor','none');
    plot3(score(t3(zero_idx),1), score(t3(zero_idx),2), score(t3(zero_idx),3), 'o', 'MarkerFaceColor', get(h3,'Color'), 'MarkerEdgeColor','none');

    grid on; view(3);
    xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
    zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
    title(sprintf('%s (n=%d)', cat(k).name, n_k));
    legend([h1,h2,h3], {'Tone','Shock','Combined'}, 'Location','best');
end

sgtitle('Population Trajectories by Region (PC1–PC3)', 'FontSize', 14);

%% SECTION: QUANTITATIVE INTERPRETATION & VISUALIZATION
% ========================================================================
% Compute numerical metrics on PCA trajectories and both print and PLOT
% them for interpretability. Metrics use PC1–PC3 scores.
% - EV(1–3): cumulative variance explained by first 3 PCs
% - Trajectory length in the response window
% - Baseline→response displacement
% - Pairwise separation (mean Euclidean) between conditions in response window
%   with bootstrap 95% CI (resampling time bins)
% - Mahalanobis distance between condition response clouds
% - Angle between displacement vectors (baseline→response centroids)

fprintf('================= QUANTITATIVE INTERPRETATION =================');
resp_loc = find(time_centers >= params.response_window(1) & time_centers <= params.response_window(2));
base_loc = find(time_centers >= params.baseline_window(1) & time_centers <= params.baseline_window(2));
nb = params.n_bins;  % bins per condition

% Colors for Tone / Shock / Combined
colTone = [0 0.4470 0.7410];
colShock = [0.8500 0.3250 0.0980];
colComb = [0 0 0];

% Storage for plotting across categories
metricsCat = struct('ok',false,'name',[], 'EV',NaN, ...
    'L',[NaN NaN NaN], 'D',[NaN NaN NaN], ...
    'Sep',[NaN NaN NaN], 'CI_low',[NaN NaN NaN], 'CI_high',[NaN NaN NaN], ...
    'MD',[NaN NaN NaN], 'ANG',[NaN NaN NaN], ...
    'cBase',[], 'cResp',[], 'score',[]);
metricsCat = repmat(metricsCat, 1, n_cat);

for k = 1:n_cat
    metricsCat(k).name = cat(k).name;
    idx = cat(k).mask; n_k = sum(idx);
    if n_k < params.min_cells_pca
        fprintf('[%s] n=%d  |  Not enough cells for PCA.', cat(k).name, n_k);
        continue;
    end

    % PCA for this category
    psth_concat = [psth_tone(idx,:), psth_shock(idx,:), psth_both(idx,:)];
    [~, score, ~, ~, explained] = pca(psth_concat');

    t1 = 1:nb; t2 = nb + (1:nb); t3 = 2*nb + (1:nb);
    S_t = score(t1,1:3); S_s = score(t2,1:3); S_b = score(t3,1:3);

    trajlen = @(M,ix) sum( sqrt(sum(diff(M(ix,:),1,1).^2,2)) );
    meandist = @(A,B,ix) mean( sqrt(sum( (A(ix,:) - B(ix,:)).^2, 2)) );
    centroid = @(M,ix) mean(M(ix,:),1);
    angdeg = @(u,v) ( acosd( max(-1,min(1, (dot(u,v)/((norm(u)*norm(v))+eps)) )) ) );

    % Lengths (response window)
    L_t = trajlen(S_t, resp_loc); L_s = trajlen(S_s, resp_loc); L_b = trajlen(S_b, resp_loc);

    % Centroids baseline/response
    cTb = centroid(S_t, base_loc); cTr = centroid(S_t, resp_loc);
    cSb = centroid(S_s, base_loc); cSr = centroid(S_s, resp_loc);
    cBb = centroid(S_b, base_loc); cBr = centroid(S_b, resp_loc);

    % Displacements
    D_t = norm(cTr - cTb); D_s = norm(cSr - cSb); D_b = norm(cBr - cBb);

    % Separations (Euclidean)
    Sep_TS = meandist(S_t, S_s, resp_loc);
    Sep_TB = meandist(S_t, S_b, resp_loc);
    Sep_SB = meandist(S_s, S_b, resp_loc);

    % Bootstrap CIs
    nboot = 500; rng('default');
    bootfun = @(A,B) meandist(A,B, resp_loc(randi(numel(resp_loc), numel(resp_loc), 1)));
    boots = @(A,B) quantile( arrayfun(@(~) bootfun(A,B), 1:nboot), [0.025 0.975]);
    CI_TS = boots(S_t,S_s); CI_TB = boots(S_t,S_b); CI_SB = boots(S_s,S_b);

    % Mahalanobis
    mdist = @(A,B) sqrt( (mean(A(resp_loc,:),1)-mean(B(resp_loc,:),1)) / (cov([A(resp_loc,:);B(resp_loc,:)])+eye(3)*1e-6) * (mean(A(resp_loc,:),1)-mean(B(resp_loc,:),1))' );
    MD_TS = mdist(S_t,S_s); MD_TB = mdist(S_t,S_b); MD_SB = mdist(S_s,S_b);

    % Angles between displacement vectors
    vT = cTr - cTb; vS = cSr - cSb; vB = cBr - cBb;
    ANG_TS = angdeg(vT, vS); ANG_TB = angdeg(vT, vB); ANG_SB = angdeg(vS, vB);

    % Printout
%     fprintf('
% [%s] n=%d  |  EV(1–3)=%.1f%%
% ', cat(k).name, n_k, sum(explained(1:3)) );
%     fprintf('  Trajectory length (response window):   Tone %.3f, Shock %.3f, Combined %.3f
% ', L_t, L_s, L_b);
%     fprintf('  Displacement |baseline→response|:      Tone %.3f, Shock %.3f, Combined %.3f
% ', D_t, D_s, D_b);
%     fprintf('  Separation mean Euclid (95%% CI) resp:  T–S %.3f [%.3f–%.3f], T–C %.3f [%.3f–%.3f], S–C %.3f [%.3f–%.3f]
% ', ...
%             Sep_TS, CI_TS(1), CI_TS(2), Sep_TB, CI_TB(1), CI_TB(2), Sep_SB, CI_SB(1), CI_SB(2));
%     fprintf('  Mahalanobis distance (resp window):    T–S %.3f, T–C %.3f, S–C %.3f
% ', MD_TS, MD_TB, MD_SB);
%     fprintf('  Angle between displacement vectors °:  T–S %.1f, T–C %.1f, S–C %.1f
% ', ANG_TS, ANG_TB, ANG_SB);

    % Store
    metricsCat(k).ok = true;
    metricsCat(k).EV = sum(explained(1:3));
    metricsCat(k).L  = [L_t L_s L_b];
    metricsCat(k).D  = [D_t D_s D_b];
    metricsCat(k).Sep = [Sep_TS Sep_TB Sep_SB];
    metricsCat(k).CI_low  = [CI_TS(1) CI_TB(1) CI_SB(1)];
    metricsCat(k).CI_high = [CI_TS(2) CI_TB(2) CI_SB(2)];
    metricsCat(k).MD  = [MD_TS MD_TB MD_SB];
    metricsCat(k).ANG = [ANG_TS ANG_TB ANG_SB];
    metricsCat(k).cBase = [cTb; cSb; cBb];
    metricsCat(k).cResp = [cTr; cSr; cBr];
    metricsCat(k).score = score; %#ok<STRNU>
end

% =====================
% Visualization of metrics
% =====================
figure('Position',[40 40 1500 800]);
for k = 1:n_cat
    % Column per category: 3 rows (vectors, Euclidean sep+CI, MD & lengths)
    if ~metricsCat(k).ok
        % Add placeholders
        subplot(n_cat,3,(k-1)*3+1); axis off; title(sprintf('%s: insufficient n', metricsCat(k).name));
        subplot(n_cat,3,(k-1)*3+2); axis off;
        subplot(n_cat,3,(k-1)*3+3); axis off;
        continue;
    end

    % Row 1: PC1–PC2 centroid displacement vectors
    subplot(n_cat,3,(k-1)*3+1);
    hold on;
    cB = metricsCat(k).cBase;  % 3x3 rows = Tone, Shock, Combined
    cR = metricsCat(k).cResp;
    % Arrows (quiver) Tone/Shock/Combined
    quiver(cB(1,1),cB(1,2), cR(1,1)-cB(1,1), cR(1,2)-cB(1,2), 0, 'Color', colTone, 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(cB(2,1),cB(2,2), cR(2,1)-cB(2,1), cR(2,2)-cB(2,2), 0, 'Color', colShock,'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(cB(3,1),cB(3,2), cR(3,1)-cB(3,1), cR(3,2)-cB(3,2), 0, 'Color', colComb, 'LineWidth', 2, 'MaxHeadSize', 0.5);
    % Baseline and response centroids
    scatter(cB(:,1), cB(:,2), 40, [colTone;colShock;colComb], 'o');
    scatter(cR(:,1), cR(:,2), 60, [colTone;colShock;colComb], 'filled');
    grid on; axis equal;
    xlabel('PC1'); ylabel('PC2');
    title(sprintf('%s: Centroid displacements (EV1–3=%.1f%%)', metricsCat(k).name, metricsCat(k).EV));
    legend({'Tone','Shock','Combined'}, 'Location','bestoutside');
    ylim([-60 60])

    % Row 2: Euclidean separations with 95% CIs
    subplot(n_cat,3,(k-1)*3+2);
    barVals = metricsCat(k).Sep; hold on;
    b = bar(1:3, barVals);
    % Error bars from CI
    ebLower = barVals - metricsCat(k).CI_low;
    ebUpper = metricsCat(k).CI_high - barVals;
    errorbar(1:3, barVals, ebLower, ebUpper, 'k', 'LineStyle','none', 'LineWidth',1.2);
    set(gca,'XTick',1:3,'XTickLabel',{'T–S','T–C','S–C'});
    ylabel('Mean Euclidean sep (resp win)'); grid on;
    title('Pairwise separations (95% CI)');
    ylim([0 200])

    % Row 3: Mahalanobis + trajectory length
    subplot(n_cat,3,(k-1)*3+3);
    yyaxis left;
    bar(1:3, metricsCat(k).MD); ylabel('Mahalanobis distance');
    set(gca,'XTick',1:6,'XTickLabel',[]); hold on; grid on;
    % overlay trajectory lengths as second axis
    yyaxis right;
    bar(4:6, metricsCat(k).L, 0.9); ylim([0 max([metricsCat(k).L, 1])]);
    set(gca,'XTick',1:6,'XTickLabel',{'MD T–S','MD T–C','MD S–C','Len T','Len S','Len C'});
    ylabel('Trajectory length (resp win)');
    title('Mahalanobis & Trajectory lengths');
    ylim([0 800])
end
sgtitle('Interpretation Metrics by Category', 'FontSize', 14);

%% ========================= HELPER FUNCTIONS ========================= %%
% function psth = compute_psth(spike_times, stim_times, params)
%     % Compute PSTH (Hz) aligned to stimulus times
%     n_trials = numel(stim_times);
%     if n_trials == 0
%         psth = zeros(1, params.n_bins);
%         return;
%     end
% 
%     edges = params.time_vec;    % time edges (s)
%     counts = zeros(n_trials, params.n_bins);
% 
%     for tr = 1:n_trials
%         rel_spikes = spike_times - stim_times(tr);
%         keep = (rel_spikes >= -params.pre_stim) & (rel_spikes <= params.post_stim);
%         tr_spk = rel_spikes(keep);
%         counts(tr,:) = histcounts(tr_spk, edges);
%     end
% 
%     fr = mean(counts,1) / params.bin_size;        % Hz
%     psth = smooth_gaussian(fr, params.smooth_sigma);
% end
% 
% function y = smooth_gaussian(x, sigma)
%     % 1D Gaussian smooth over samples; sigma in samples (bins)
%     if sigma <= 0
%         y = x; return;
%     end
%     w = ceil(sigma*6);
%     if mod(w,2)==0, w = w+1; end
%     half = (w-1)/2;
%     xx = -half:half;
%     k = exp(-(xx.^2)/(2*sigma^2));
%     k = k/sum(k);
%     y = conv(x, k, 'same');
% end
function psth = compute_psth(spike_times, stim_times, params)
    % Compute PSTH (Hz) aligned to stimulus times
    n_trials = numel(stim_times);
    if n_trials == 0
        psth = zeros(1, params.n_bins);
        return;
    end

    edges = params.time_vec;    % time edges (s)
    counts = zeros(n_trials, params.n_bins);

    for tr = 1:n_trials
        rel_spikes = spike_times - stim_times(tr);
        keep = (rel_spikes >= -params.pre_stim) & (rel_spikes <= params.post_stim);
        tr_spk = rel_spikes(keep);
        counts(tr,:) = histcounts(tr_spk, edges);
    end

    fr = mean(counts,1) / params.bin_size;        % Hz
    psth = smooth_gaussian(fr, params.smooth_sigma);
end

function y = smooth_gaussian(x, sigma)
    % 1D Gaussian smooth over samples; sigma in samples (bins)
    if sigma <= 0
        y = x; return;
    end
    w = ceil(sigma*6);
    if mod(w,2)==0, w = w+1; end
    half = (w-1)/2;
    xx = -half:half;
    k = exp(-(xx.^2)/(2*sigma^2));
    k = k/sum(k);
    y = conv(x, k, 'same');
end
