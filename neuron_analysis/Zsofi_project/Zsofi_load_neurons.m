clear all
mainFolder = 'Z:\HajosLab\Dani\Magyar Dániel\Analysis2\Matlab_files';
cd(mainFolder)
ii = 1;
% MD180506_R - no baseline activity

% MD190211_L_1283 - no baseline activity, only 1-1 spikes after each TTL

% MD190212_R_1056 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190212_R_1056';
cell_metrics.brainRegion{ii} = 'LAp';
filename = 'md190212_jobb_1056_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb_1056_nohexa_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb_1056_nohexa_Ch31.times;
clearvars -except cell_metrics

% MD190212_R_588 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190212_R_588';
cell_metrics.brainRegion{ii} = 'LAp';
filename = 'md190212_jobb_588_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb_588_nohexa_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb_588_nohexa_Ch31.times;
clearvars -except cell_metrics

% MD190322_L_3062 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190322_L_3062';
cell_metrics.brainRegion{ii} = 'LAp';
filename = 'md190322_bal_3062_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_3062_nohexa_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_3062_nohexa_Ch31.times;
clearvars -except cell_metrics

% MD190306_R_-4010 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190306_R_-4010';
cell_metrics.brainRegion{ii} = 'LAp';
filename = 'md190306_jobb_-4010_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb__4010_nohexa_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb__4010_nohexa_Ch31.times;
clearvars -except cell_metrics

% MD190209_R_4452 - 2 units, I don't know which one was filled

% MD190218_R_9789 - OK
ii = ii+1;
disp(ii)
cell_metrics.cellID{ii} = 'MD190218_R_9789';
cell_metrics.brainRegion{ii} = 'BAa';
filename = 'md190218_jobb_9789.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb_9789_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb_9789_Ch31.times;
clearvars -except cell_metrics

% MD190315_L_14528 - 2 units, I don't know which one was filled

% MD181017_L_5510 - OK, jo pelda sejt lehet, van baseline es toltes, de nem rekonstrualt
ii = ii+1;
cell_metrics.cellID{ii} = 'MD181017_L_5510';
cell_metrics.brainRegion{ii} = 'BAa';
filename = 'md181017_bal5510_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_5510_no_hexa_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_5510_no_hexa_Ch31.times;
clearvars -except cell_metrics

% MD190613_R_12383 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190613_R_12383';
cell_metrics.brainRegion{ii} = 'BAa';
filename = 'md190613_jobb_12383.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb_12383_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb_12383_Ch31.times;
clearvars -except cell_metrics

% MD180710_L_5211 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD180710_L_5211';
cell_metrics.brainRegion{ii} = 'BAa';
filename = 'md180710_bal5211_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_5211_no_hexa_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_5211_no_hexa_Ch31.times;
clearvars -except cell_metrics

% MD180504_R_34109 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD180504_R_34109';
cell_metrics.brainRegion{ii} = 'BAa';
filename = 'md180504_jobb34109_stimulus.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb_34109_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb_34109_Ch31.times;
clearvars -except cell_metrics

% MD190612_R_8475 - OK
ii = ii+1;
disp(ii)
cell_metrics.cellID{ii} = 'MD190612_R_8475';
cell_metrics.brainRegion{ii} = 'BAa';
filename = 'md190612_jobb_8475.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb_8475_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb_8475_Ch31.times;
clearvars -except cell_metrics

% MD190511_R_9149 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190511_R_9149';
cell_metrics.brainRegion{ii} = 'BAa';
filename = 'md190511_jobb_9149.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb_9149_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb_9149_Ch31.times;
clearvars -except cell_metrics

% MD180426_R_5700 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD180426_R_5700';
cell_metrics.brainRegion{ii} = 'BAa';
filename = 'md180426_jobb5700_stimulus.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb__5700_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb__5700_Ch31.times;
clearvars -except cell_metrics

% MD181016_L_4296 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD181016_L_4296';
cell_metrics.brainRegion{ii} = 'BAa';
filename = 'md181016_bal4296_EZAJOO.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_4296_csak_shockolgat_s2_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_4296_csak_shockolgat_s2_Ch31.times;
clearvars -except cell_metrics

% MD190218_L_11020 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190218_L_11020';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md190218_bal_11020_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_11020_nohexa_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_11020_nohexa_Ch31.times;
clearvars -except cell_metrics

% MD180327_L_3250 - OK
ii = ii+1;
disp(ii)
cell_metrics.cellID{ii} = 'MD180327_L_3250';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md180327_bal_3250_stimulus.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_3250_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_3250_Ch31.times;
clearvars -except cell_metrics

% MD180731_L_6383 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD180731_L_6383';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md180731_bal6383_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_6383_no_hexa_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_6383_no_hexa_Ch31.times;
clearvars -except cell_metrics

% MD180426_L_5555 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD180426_L_5555';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md180426_bal5555_stimulus.mat';
load(filename)
cell_metrics.spikes.times{ii} = ba__5555_Ch13.times;
cell_metrics.stimulus.times{ii} = ba__5555_Ch31.times;
clearvars -except cell_metrics

% MD180504_L_34120 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD180504_L_34120';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md180504_bal34120_stimulus.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_34120_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_34120_Ch31.times;
clearvars -except cell_metrics

% MD190305_L_7885 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190305_L_7885';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md190305_bal_7885_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_7885_nohexa_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_7885_nohexa_Ch31.times;
clearvars -except cell_metrics

% MD190328_L_9300 - OK
ii = ii+1;
disp(ii)
cell_metrics.cellID{ii} = 'MD190328_L_9300';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md190328_bal_9300.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_9300_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_9300_Ch31.times;
clearvars -except cell_metrics

% MD190316_L_10056 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190316_L_10056';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md190316_bal_10056_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_10056_nohexa_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_10056_nohexa_Ch31.times;
clearvars -except cell_metrics

% MD190617_L_9716 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190617_L_9716';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md190617_bal_9716.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_9716_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_9716_Ch31.times;
clearvars -except cell_metrics

% MD190209_L_-1849 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190209_L_-1849';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md190209_bal_-1849_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal__1849_nohexa_Ch13.times;
cell_metrics.stimulus.times{ii} = bal__1849_nohexa_Ch31.times;
clearvars -except cell_metrics

% MD190212_L_-778 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190212_L_-778';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md190212_bal_-778_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal__778_nohexa_Ch13.times;
cell_metrics.stimulus.times{ii} = bal__778_nohexa_Ch31.times;
clearvars -except cell_metrics

% MD190510_L_9230 - OK, EZ IS JO, talan ez jobb pelda sejtnek
ii = ii+1;
disp(ii)
cell_metrics.cellID{ii} = 'MD190510_L_9230';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md190510_bal_9230.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_9230_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_9230_Ch31.times;
clearvars -except cell_metrics

% MD190512_R_6666 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190512_R_6666';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md190512_jobb_6666.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb_6666_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb_6666_Ch31.times;
clearvars -except cell_metrics

% MD180510_R_29921 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD180510_R_29921';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md180510_jobb29921_stimulus.mat';
load(filename)
cell_metrics.spikes.times{ii} = jobb_29921_Ch13.times;
cell_metrics.stimulus.times{ii} = jobb_29921_Ch31.times;
clearvars -except cell_metrics

% MD190208_L_1909 - OK
ii = ii+1;
cell_metrics.cellID{ii} = 'MD190208_L_1909';
cell_metrics.brainRegion{ii} = 'BAp';
filename = 'md190208_bal1909_nohexa.mat';
load(filename)
cell_metrics.spikes.times{ii} = bal_1909_nohexa_Ch13.times;
cell_metrics.stimulus.times{ii} = bal_1909_nohexa_Ch31.times;
clearvars -except cell_metrics

savefolder = 'C:\Users\dmagyar\My Drive\Papers\Zsofi';
save([savefolder '\cell_metrics.mat'], 'cell_metrics')