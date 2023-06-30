session = sessionTemplate(cd);
%% CHETA PFC
load('TTLsKS.mat', 'ChETA_50_20Hz')
laserChETA_50_20Hz.timestamps = ChETA_50_20Hz;
laserChETA_50_20Hz.timestamps(:,2) = ChETA_50_20Hz+5;
save('temp_wh.laserChETA_50_20Hz.manipulation.mat', 'laserChETA_50_20Hz')
cell_metrics = ProcessCellMetrics('session', session);
% Cheta Running 
session = sessionTemplate(cd);
load('TTLsKS.mat', 'ChETA_50_20Hz_Running')
laserChETA_50_20Hz_Running.timestamps = ChETA_50_20Hz_Running;
laserChETA_50_20Hz_Running.timestamps(:,2) = ChETA_50_20Hz_Running+5;
save('temp_wh.laserChETA_50_20Hz_Running.manipulation.mat', 'laserChETA_50_20Hz_Running')
cell_metrics = ProcessCellMetrics('session', session);

% Cheta Resting
session = sessionTemplate(cd);
load('TTLsKS.mat', 'ChETA_50_20Hz_Resting')
laserChETA_50_20Hz_Resting.timestamps = ChETA_50_20Hz_Resting;
laserChETA_50_20Hz_Resting.timestamps(:,2) = ChETA_50_20Hz_Resting+5;
save('temp_wh.laserChETA_50_20Hz_Resting.manipulation.mat', 'laserChETA_50_20Hz_Resting')
cell_metrics = ProcessCellMetrics('session', session);

%% CHRIMSON + STGTACR PFC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('TTLsKS.mat', 'BA_25_5Hz')
laserBA_25_5Hz_short.timestamps = BA_25_5Hz;
laserBA_25_5Hz_short.timestamps(:,2) = BA_25_5Hz+2;
save('temp_wh.laserBA_25_5Hz_short.manipulation.mat', 'laserBA_25_5Hz_short')
cell_metrics = ProcessCellMetrics('session', session);

load('TTLsKS.mat', 'TO_25_5Hz')
laserTO_25_5Hz_short.timestamps = TO_25_5Hz;
laserTO_25_5Hz_short.timestamps(:,2) = TO_25_5Hz+2;
save('temp_wh.laserTO_25_5Hz_short.manipulation.mat', 'laserTO_25_5Hz_short')
cell_metrics = ProcessCellMetrics('session', session);

load('TTLsKS.mat', 'BA_25_5Hz')
laserBA_25_5Hz_long.timestamps = BA_25_5Hz;
laserBA_25_5Hz_long.timestamps(:,2) = BA_25_5Hz+10;
save('temp_wh.laserBA_25_5Hz_long.manipulation.mat', 'laserBA_25_5Hz_long')
cell_metrics = ProcessCellMetrics('session', session);

load('TTLsKS.mat', 'TO_25_5Hz')
laserTO_25_5Hz_long.timestamps = TO_25_5Hz;
laserTO_25_5Hz_long.timestamps(:,2) = TO_25_5Hz+10;
save('temp_wh.laserTO_25_5Hz_long.manipulation.mat', 'laserTO_25_5Hz_long')
cell_metrics = ProcessCellMetrics('session', session);

% Running all
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_25_all_Running')
laserBA_25_all_Running.timestamps = BA_25_all_Running;
laserBA_25_all_Running.timestamps(:,2) = BA_25_all_Running+5;
save('temp_wh.laserBA_25_all_Running.manipulation.mat', 'laserBA_25_all_Running')
cell_metrics = ProcessCellMetrics('session', session);

% Running long
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_25_5Hz_Running')
laserBA_25_5Hz_Running_long.timestamps = BA_25_5Hz_Running;
laserBA_25_5Hz_Running_long.timestamps(:,2) = BA_25_5Hz_Running+10;
save('temp_wh.laserBA_25_5Hz_Running_long.manipulation.mat', 'laserBA_25_5Hz_Running_long')
cell_metrics = ProcessCellMetrics('session', session);

session = sessionTemplate(cd);
load('TTLsKS.mat', 'TO_25_5Hz_Running')
laserTO_25_5Hz_Running_long.timestamps = TO_25_5Hz_Running;
laserTO_25_5Hz_Running_long.timestamps(:,2) = TO_25_5Hz_Running+10;
save('temp_wh.laserTO_25_5Hz_Running_long.manipulation.mat', 'laserTO_25_5Hz_Running_long')
cell_metrics = ProcessCellMetrics('session', session);

% Resting all
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_25_all_Resting')
laserBA_25_all_Resting.timestamps = BA_25_all_Resting;
laserBA_25_all_Resting.timestamps(:,2) = BA_25_all_Resting+5;
save('temp_wh.laserBA_25_all_Resting.manipulation.mat', 'laserBA_25_all_Resting')
cell_metrics = ProcessCellMetrics('session', session);

% Resting long
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_25_5Hz_Resting')
laserBA_25_5Hz_Resting_long.timestamps = BA_25_5Hz_Resting;
laserBA_25_5Hz_Resting_long.timestamps(:,2) = BA_25_5Hz_Resting+10;
save('temp_wh.laserBA_25_5Hz_Resting_long.manipulation.mat', 'laserBA_25_5Hz_Resting_long')
cell_metrics = ProcessCellMetrics('session', session);

session = sessionTemplate(cd);
load('TTLsKS.mat', 'TO_25_5Hz_Resting')
laserTO_25_5Hz_Resting_long.timestamps = TO_25_5Hz_Resting;
laserTO_25_5Hz_Resting_long.timestamps(:,2) = TO_25_5Hz_Resting+10;
save('temp_wh.laserTO_25_5Hz_Resting_long.manipulation.mat', 'laserTO_25_5Hz_Resting_long')
cell_metrics = ProcessCellMetrics('session', session);

% Resting supershort 
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_25_5Hz_Resting')
laserBA_25_5Hz_Resting_supershort.timestamps = BA_25_5Hz_Resting;
laserBA_25_5Hz_Resting_supershort.timestamps(:,2) = BA_25_5Hz_Resting+0.1;
save('temp_wh.laserBA_25_5Hz_Resting_supershort.manipulation.mat', 'laserBA_25_5Hz_Resting_supershort')
cell_metrics = ProcessCellMetrics('session', session);

session = sessionTemplate(cd);
load('TTLsKS.mat', 'TO_25_5Hz_Resting')
laserTO_25_5Hz_Resting_supershort.timestamps = TO_25_5Hz_Resting;
laserTO_25_5Hz_Resting_supershort.timestamps(:,2) = TO_25_5Hz_Resting+0.1;
save('temp_wh.laserTO_25_5Hz_Resting_supershort.manipulation.mat', 'laserTO_25_5Hz_Resting_supershort')
cell_metrics = ProcessCellMetrics('session', session);

%250
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_250_5Hz')
laserBA_250_5Hz_all.timestamps = BA_250_5Hz;
laserBA_250_5Hz_all.timestamps(:,2) = BA_250_5Hz+0.1;
save('temp_wh.laserBA_250_5Hz_all.manipulation.mat', 'laserBA_250_5Hz_all')
cell_metrics = ProcessCellMetrics('session', session);

load('TTLsKS.mat', 'TO_250_5Hz')
laserTO_250_5Hz_all.timestamps = TO_250_5Hz;
laserTO_250_5Hz_all.timestamps(:,2) = TO_250_5Hz+0.1;
save('temp_wh.laserTO_250_5Hz_all.manipulation.mat', 'laserTO_250_5Hz_all')
cell_metrics = ProcessCellMetrics('session', session);





%% BA_50_MD130-133 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_50_5Hz')
laserBA_50_5Hz.timestamps = BA_50_5Hz;
laserBA_50_5Hz.timestamps(:,2) = BA_50_5Hz+5;
save('temp_wh.laserBA_50_5Hz.manipulation.mat', 'laserBA_50_5Hz')
cell_metrics = ProcessCellMetrics('session', session);

% BA_50_MD130-133 Running 
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_50_5Hz_Running')
laserBA_50_5Hz_Running.timestamps = BA_50_5Hz_Running;
laserBA_50_5Hz_Running.timestamps(:,2) = BA_50_5Hz_Running+5;
save('temp_wh.laserBA_50_5Hz_Running.manipulation.mat', 'laserBA_50_5Hz_Running')
cell_metrics = ProcessCellMetrics('session', session);

% BA_50_MD130-133 Resting
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_50_5Hz_Resting')
laserBA_50_5Hz_Resting.timestamps = BA_50_5Hz_Resting;
laserBA_50_5Hz_Resting.timestamps(:,2) = BA_50_5Hz_Resting+5;
save('temp_wh.laserBA_50_5Hz_Resting.manipulation.mat', 'laserBA_50_5Hz_Resting')
cell_metrics = ProcessCellMetrics('session', session);



%% M2 footshock pre1post3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session = sessionTemplate(cd);
load('TTLsKS.mat', 'shocks')
shocksTTL.timestamps = shocks;
shocksTTL.timestamps(:,2) = shocks+5;
save('temp_wh.shocksTTL.manipulation.mat', 'shocksTTL')
cell_metrics = ProcessCellMetrics('session', session);

%% PFC Footshock
session = sessionTemplate(cd);
load('TTLsKS.mat', 'shock_only', 'shock_inh')
shocks_onlyTTL.timestamps = shock_only;
%save('temp_wh.shocks_onlyTTL.manipulation.mat', 'shocks_onlyTTL')
%load('TTLsKS.mat', 'shocks_inh')
shocks_inhTTL.timestamps = shock_inh;
save('temp_wh.shocks_onlyTTL.manipulation.mat', 'shocks_onlyTTL')
save('temp_wh.shocks_inhTTL.manipulation.mat', 'shocks_inhTTL')
cell_metrics = ProcessCellMetrics('session', session);

% PFC Footshock Running -Resting
session = sessionTemplate(cd);
load('TTLsKS.mat', 'shock_only_Running') % shock only Running 
lasershock_only_Running.timestamps = shock_only_Running;
lasershock_only_Running.timestamps(:,2) = shock_only_Running + 5;
save('temp_wh.lasershock_only_Running.manipulation.mat', 'lasershock_only_Running')
load('TTLsKS.mat', 'shock_only_Resting') % shock only resting
lasershock_only_Resting.timestamps = shock_only_Resting;
lasershock_only_Resting.timestamps(:,2) = shock_only_Resting + 5;
save('temp_wh.lasershock_only_Resting.manipulation.mat', 'lasershock_only_Resting')

load('TTLsKS.mat', 'shock_inh_Running') % shock inh Running 
lasershock_inh_Running.timestamps = shock_inh_Running;
lasershock_inh_Running.timestamps(:,2) = shock_inh_Running + 5;
save('temp_wh.lasershock_inh_Running.manipulation.mat', 'lasershock_inh_Running')
load('TTLsKS.mat', 'shock_inh_Resting') % shock inh resting
lasershock_inh_Resting.timestamps = shock_inh_Resting;
lasershock_inh_Resting.timestamps(:,2) = shock_inh_Resting + 5;
save('temp_wh.lasershock_inh_Resting.manipulation.mat', 'lasershock_inh_Resting')
cell_metrics = ProcessCellMetrics('session', session);

%%PFC shock run rest OEPv063
session = sessionTemplate(cd);
load('TTLsKS.mat', 'shock_motor_rest') % shock only Running 
TTLshock_motor_rest.timestamps = shock_motor_rest;
TTLshock_motor_rest.timestamps(:,2) = shock_motor_rest + 2.5;
save('temp_wh.TTLshock_motor_rest.manipulation.mat', 'TTLshock_motor_rest')
load('TTLsKS.mat', 'shock_motor_run') % shock only resting
TTLshock_motor_run.timestamps = shock_motor_run;
TTLshock_motor_run.timestamps(:,2) = shock_motor_run + 2.5;
save('temp_wh.TTLshock_motor_run.manipulation.mat', 'TTLshock_motor_run')
cell_metrics = ProcessCellMetrics('session', session);

%%PFC shock run rest compare OEPv063
session = sessionTemplate(cd);
load('TTLsKS.mat', 'shock_motor_rest') % shock only Running 
TTLshock_motor_rest.timestamps = shock_motor_rest;
TTLshock_motor_rest.timestamps(:,2) = shock_motor_rest + 2.5;
save('temp_wh.TTLshock_motor_rest.manipulation.mat', 'TTLshock_motor_rest')
load('TTLsKS.mat', 'shock_motor_run_small') % shock only resting
TTLshock_motor_run_small.timestamps = shock_motor_run_small;
TTLshock_motor_run_small.timestamps(:,2) = shock_motor_run_small + 2.5;
save('temp_wh.TTLshock_motor_run_small.manipulation.mat', 'TTLshock_motor_run_small')
load('TTLsKS.mat', 'shock_motor_run_high') % shock only resting
TTLshock_motor_run_high.timestamps = shock_motor_run_high;
TTLshock_motor_run_high.timestamps(:,2) = shock_motor_run_high + 2.5;
save('temp_wh.TTLshock_motor_run_high.manipulation.mat', 'TTLshock_motor_run_high')
cell_metrics = ProcessCellMetrics('session', session);

%%PFC BAopto run rest OEPv063
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_5Hz_motor_rest_first') 
TTLBA_5Hz_motor_rest_first.timestamps = BA_5Hz_motor_rest_first;
TTLBA_5Hz_motor_rest_first.timestamps(:,2) = BA_5Hz_motor_rest_first + 2.5;
save('temp_wh.TTLBA_5Hz_motor_rest_first.manipulation.mat', 'TTLBA_5Hz_motor_rest_first')

load('TTLsKS.mat', 'BA_20Hz_motor_rest_first')
TTLBA_20Hz_motor_rest_first.timestamps = BA_20Hz_motor_rest_first;
TTLBA_20Hz_motor_rest_first.timestamps(:,2) = BA_20Hz_motor_rest_first + 2.5;
save('temp_wh.TTLBA_20Hz_motor_rest_first.manipulation.mat', 'TTLBA_20Hz_motor_rest_first')

load('TTLsKS.mat', 'BA_5Hz_motor_run_first')
TTLBA_5Hz_motor_run_first.timestamps = BA_5Hz_motor_run_first;
TTLBA_5Hz_motor_run_first.timestamps(:,2) = BA_5Hz_motor_run_first + 2.5;
save('temp_wh.TTLBA_5Hz_motor_run_first.manipulation.mat', 'TTLBA_5Hz_motor_run_first')

load('TTLsKS.mat', 'BA_20Hz_motor_run_first') 
TTLBA_20Hz_motor_run_first.timestamps = BA_20Hz_motor_run_first;
TTLBA_20Hz_motor_run_first.timestamps(:,2) = BA_20Hz_motor_run_first + 2.5;
save('temp_wh.TTLBA_20Hz_motor_run_first.manipulation.mat', 'TTLBA_20Hz_motor_run_first')

load('TTLsKS.mat', 'BA_5Hz_motor_rest_all') 
TTLBA_5Hz_motor_rest_all.timestamps = BA_5Hz_motor_rest_all;
TTLBA_5Hz_motor_rest_all.timestamps(:,2) = BA_5Hz_motor_rest_all + 0.1;
save('temp_wh.TTLBA_5Hz_motor_rest_all.manipulation.mat', 'TTLBA_5Hz_motor_rest_all')

load('TTLsKS.mat', 'BA_5Hz_motor_run_all') 
TTLBA_5Hz_motor_run_all.timestamps = BA_5Hz_motor_run_all;
TTLBA_5Hz_motor_run_all.timestamps(:,2) = BA_5Hz_motor_run_all + 0.1;
save('temp_wh.TTLBA_5Hz_motor_run_all.manipulation.mat', 'TTLBA_5Hz_motor_run_all')

load('TTLsKS.mat', 'BA_20Hz_motor_rest_all') 
TTLBA_20Hz_motor_rest_all.timestamps = BA_20Hz_motor_rest_all;
TTLBA_20Hz_motor_rest_all.timestamps(:,2) = BA_20Hz_motor_rest_all + 0.025;
save('temp_wh.TTLBA_20Hz_motor_rest_all.manipulation.mat', 'TTLBA_20Hz_motor_rest_all')

load('TTLsKS.mat', 'BA_20Hz_motor_run_all') 
TTLBA_20Hz_motor_run_all.timestamps = BA_20Hz_motor_run_all;
TTLBA_20Hz_motor_run_all.timestamps(:,2) = BA_20Hz_motor_run_all + 0.025;
save('temp_wh.TTLBA_20Hz_motor_run_all.manipulation.mat', 'TTLBA_20Hz_motor_run_all')

cell_metrics = ProcessCellMetrics('session', session);





% chrimson + stgtacr
basenames = {'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december\Kilosort_v2\MD098_kilosort\kilosort3preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december\Kilosort_v2\MD099_kilosort\kilosort3preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december\Kilosort_v2\MD101_kilosort\kilosort3preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december\Kilosort_v2\MD108_kilosort\kilosort3preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);

% chrimson + stgtacr KS 25
basenames = {'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\PFC_project_base\Experiments\ChrimsonStGtACR\MD098_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_project_base\Experiments\ChrimsonStGtACR\MD099_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_project_base\Experiments\ChrimsonStGtACR\MD101_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_project_base\Experiments\ChrimsonStGtACR\MD108_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);

%%PFC shock run rest OEPv063
basenames = {'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\PFC_shock_run_rest\MD156_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_shock_run_rest\MD157_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_shock_run_rest\MD158_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_shock_run_rest\MD159_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);

%%PFC BAopto run rest OEPv063
basenames = {'temp_wh','temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\BAopto_run_rest\MD173_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BAopto_run_rest\MD174_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BAopto_run_rest\MD175_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BAopto_run_rest\MD176_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BAopto_run_rest\MD177_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PFC FOOTSHOCK
basenames = {'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD138_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD139_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD140_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD141_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD142_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_PFC_shock_waveform\kilosort\MD143_kilosort\kilosort3preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);



% PFC CHRIMSON ONLY
basenames = {'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\2022_august_Chrimson\kilosort\MD130_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_Chrimson\kilosort\MD131_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\2022_august_Chrimson\kilosort\MD132_kilosort\kilosort3preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);

% PFC ChETA
basenames = {'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\ChETA\2022_august\kilosort\MD133_kilosort\kilosort3preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\ChETA\2022_august\kilosort\MD134_kilosort\kilosort3preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\ChETA\2022_august\kilosort\MD135_kilosort\kilosort3preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\ChETA\2022_august\kilosort\MD136_kilosort\kilosort3preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\ChETA\2022_august\kilosort\MD137_kilosort\kilosort3preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);

% M2 fottshock only
basenames = {'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\M2_shock_waveform\kilosort\MD127_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_shock_waveform\kilosort\MD128_kilosort\kilosort3preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_shock_waveform\kilosort\MD129_kilosort\kilosort3preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);


cell_metrics.animal(1:216)={'MD098'};
save('temp_wh.cell_metrics.cellinfo.mat', 'cell_metrics')

cell_metrics.animal(1:89)={'MD099'};
save('temp_wh.cell_metrics.cellinfo.mat', 'cell_metrics')

cell_metrics.animal(1:137)={'MD101'};
save('temp_wh.cell_metrics.cellinfo.mat', 'cell_metrics')

cell_metrics.animal(1:96)={'MD108'};
save('temp_wh.cell_metrics.cellinfo.mat', 'cell_metrics')
