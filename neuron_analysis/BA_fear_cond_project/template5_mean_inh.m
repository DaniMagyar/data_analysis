cell_metrics = BAfc_load_neurons('NP_BAfc_triptest');
% Find location of neurons
idx_LA = find(contains(cell_metrics.brainRegion, 'LA'));
idx_BA = find(contains(cell_metrics.brainRegion, 'BA'));

% Find excited neurons
[~, idx_sounds_100, ~] =   BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_sound_only',    'inh', 2, 0.5, [0.013 0.1], 0.001, 'artefactLength',     0);
[~, idx_shocks_100, ~] =   BAfc_find_response('NP_BAfc_triptest', 'shock','TTL_triptest_shocks_only',       'inh', 2, 0.5, [0.013 0.1], 0.001, 'artefactLength', 0.012);
[~, idx_both_100, ~] =     BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_both',          'inh', 2, 0.5, [0.013 0.1], 0.001, 'artefactLength', 0.012);

[~, idx_sounds_200, ~] =   BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_sound_only',    'inh', 2, 0.5, [0.1 0.2], 0.001, 'artefactLength',     0);
[~, idx_shocks_200, ~] =   BAfc_find_response('NP_BAfc_triptest', 'shock','TTL_triptest_shocks_only',       'inh', 2, 0.5, [0.1 0.2], 0.001, 'artefactLength', 0.012);
[~, idx_both_200, ~] =     BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_both',          'inh', 2, 0.5, [0.1 0.2], 0.001, 'artefactLength', 0.012);

[~, idx_sounds_300, ~] =   BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_sound_only',    'inh', 2, 0.5, [0.2 0.3], 0.001, 'artefactLength',     0);
[~, idx_shocks_300, ~] =   BAfc_find_response('NP_BAfc_triptest', 'shock','TTL_triptest_shocks_only',       'inh', 2, 0.5, [0.2 0.3], 0.001, 'artefactLength', 0.012);
[~, idx_both_300, ~] =     BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_both',          'inh', 2, 0.5, [0.2 0.3], 0.001, 'artefactLength', 0.012);

[~, idx_sounds_400, ~] =   BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_sound_only',    'inh', 2, 0.5, [0.3 0.4], 0.001, 'artefactLength',     0);
[~, idx_shocks_400, ~] =   BAfc_find_response('NP_BAfc_triptest', 'shock','TTL_triptest_shocks_only',       'inh', 2, 0.5, [0.3 0.4], 0.001, 'artefactLength', 0.012);
[~, idx_both_400, ~] =     BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_both',          'inh', 2, 0.5, [0.3 0.4], 0.001, 'artefactLength', 0.012);

[~, idx_sounds_500, ~] =   BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_sound_only',    'inh', 2, 0.5, [0.4 0.5], 0.001, 'artefactLength',     0);
[~, idx_shocks_500, ~] =   BAfc_find_response('NP_BAfc_triptest', 'shock','TTL_triptest_shocks_only',       'inh', 2, 0.5, [0.4 0.5], 0.001, 'artefactLength', 0.012);
[~, idx_both_500, ~] =     BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_both',          'inh', 2, 0.5, [0.4 0.5], 0.001, 'artefactLength', 0.012);

idx_sounds = unique([idx_sounds_100; idx_sounds_200; idx_sounds_300; idx_sounds_400; idx_sounds_500]);
idx_shocks = unique([idx_shocks_100; idx_shocks_200; idx_shocks_300; idx_shocks_400; idx_shocks_500]);
idx_both = unique([idx_both_100; idx_both_200; idx_both_300; idx_both_400; idx_both_500]);

[psth_spx_sounds, ~, ~] =   BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_sound_only',    'inh', 0, 1, [0.4 0.5], 0.001, 'artefactLength',     0);
[psth_spx_shocks, ~, ~] =   BAfc_find_response('NP_BAfc_triptest', 'shock','TTL_triptest_shocks_only',       'inh', 0, 1, [0.4 0.5], 0.001, 'artefactLength', 0.012);
[psth_spx_both, ~, ~] =     BAfc_find_response('NP_BAfc_triptest', 'allsound', 'TTL_triptest_both',          'inh', 0, 1, [0.4 0.5], 0.001, 'artefactLength', 0.012);

psth_spx_sounds = zscore(psth_spx_sounds,[],2);
psth_spx_shocks = zscore(psth_spx_shocks,[],2);
psth_spx_both = zscore(psth_spx_both,[],2);

tiledlayout(1,2)
%% LA
nexttile
% Extract the selected rows
selected_rows_sounds = psth_spx_sounds(intersect(idx_sounds, idx_LA), :);
selected_rows_shocks = psth_spx_shocks(intersect(idx_shocks, idx_LA), :);
selected_rows_both = psth_spx_both(intersect(idx_both, idx_LA), :);

% Compute the mean across selected rows
mean_sounds = smoothdata(mean(selected_rows_sounds, 1), 'movmean', 5);
mean_shocks = smoothdata(mean(selected_rows_shocks, 1), 'movmean', 5);
mean_both = smoothdata(mean(selected_rows_both, 1), 'movmean', 5);

% Plot the mean values
hold on;
plot(mean_sounds, 'b', 'LineWidth', 2, 'DisplayName', 'Sounds');
plot(mean_shocks, 'r', 'LineWidth', 2, 'DisplayName', 'Shocks');
plot(mean_both, 'g', 'LineWidth', 2, 'DisplayName', 'Both');

% Formatting
xlabel('Time Bins');
ylabel('Mean Response');
title('LA');
legend;
grid on;
hold off;

%% BA
nexttile
% Extract the selected rows
selected_rows_sounds = psth_spx_sounds(intersect(idx_sounds, idx_BA), :);
selected_rows_shocks = psth_spx_shocks(intersect(idx_shocks, idx_BA), :);
selected_rows_both = psth_spx_both(intersect(idx_both, idx_BA), :);

% Compute the mean across selected rows
mean_sounds = smoothdata(mean(selected_rows_sounds, 1), 'movmean', 5);
mean_shocks = smoothdata(mean(selected_rows_shocks, 1), 'movmean', 5);
mean_both = smoothdata(mean(selected_rows_both, 1), 'movmean', 5);

% Plot the mean values

hold on;
plot(mean_sounds, 'b', 'LineWidth', 2, 'DisplayName', 'Sounds');
plot(mean_shocks, 'r', 'LineWidth', 2, 'DisplayName', 'Shocks');
plot(mean_both, 'g', 'LineWidth', 2, 'DisplayName', 'Both');

% Formatting
xlabel('Time Bins');
ylabel('Mean Response');
title('BA');
legend;
grid on;
hold off;
