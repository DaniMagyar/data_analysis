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
shocksTTL.timestamps(:,2) = shocks+1;
save('temp_wh.shocksTTL.manipulation.mat', 'shocksTTL')
load('TTLsKS.mat', 'shocks_drug')
shocks_drugTTL.timestamps = shocks_drug;
shocks_drugTTL.timestamps(:,2) = shocks_drug+1;
save('temp_wh.shocks_drugTTL.manipulation.mat', 'shocks_drugTTL')


cell_metrics = ProcessCellMetrics('session', session);

%% M2 footshock light(Arch) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session = sessionTemplate(cd);
load('TTLsKS.mat', 'shocks_only')
shocks_onlyTTL.timestamps = shocks_only;
shocks_onlyTTL.timestamps(:,2) = shocks_only+1;
save('temp_wh.shocks_onlyTTL.manipulation.mat', 'shocks_onlyTTL')
load('TTLsKS.mat', 'shocks_light')
shocks_lightTTL.timestamps = shocks_light;
shocks_lightTTL.timestamps(:,2) = shocks_light+1;
save('temp_wh.shocks_lightTTL.manipulation.mat', 'shocks_lightTTL')


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
load('TTLsKS.mat', 'motorStart') % shock only resting
TTLmotorStart.timestamps = motorStart;
TTLmotorStart.timestamps(:,2) = motorStart + 2.5;
save('temp_wh.TTLmotorStart.manipulation.mat', 'TTLmotorStart')
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

load('TTLsKS.mat', 'BA_5Hz_all_first') 
TTLBA_5Hz_all_first.timestamps = BA_5Hz_all_first;
TTLBA_5Hz_all_first.timestamps(:,2) = BA_5Hz_all_first + 2.5;
save('temp_wh.TTLBA_5Hz_all_first.manipulation.mat', 'TTLBA_5Hz_all_first')


cell_metrics = ProcessCellMetrics('session', session);

% PFC single pulse and shock OEPv063
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_singles_all') 
TTLBA_singles_all.timestamps = BA_singles_all;
TTLBA_singles_all.timestamps(:,2) = BA_singles_all + 0.1;
save('temp_wh.TTLBA_singles_all.manipulation.mat', 'TTLBA_singles_all')

load('TTLsKS.mat', 'TO_singles_all') 
TTLTO_singles_all.timestamps = TO_singles_all;
TTLTO_singles_all.timestamps(:,2) = TO_singles_all + 0.1;
save('temp_wh.TTLTO_singles_all.manipulation.mat', 'TTLTO_singles_all')

load('TTLsKS.mat', 'BA_only_trains_first') 
TTLBA_only_trains_first.timestamps = BA_only_trains_first;
TTLBA_only_trains_first.timestamps(:,2) = BA_only_trains_first + 2.5;
save('temp_wh.TTLBA_only_trains_first.manipulation.mat', 'TTLBA_only_trains_first')

load('TTLsKS.mat', 'shock_only') 
TTLshock_only.timestamps = shock_only;
TTLshock_only.timestamps(:,2) = shock_only + 2.5;
save('temp_wh.TTLshock_only.manipulation.mat', 'TTLshock_only')

load('TTLsKS.mat', 'BA_shock_trains_first') 
TTLBA_shock_trains_first.timestamps = BA_shock_trains_first;
TTLBA_shock_trains_first.timestamps(:,2) = BA_shock_trains_first + 2.5;
save('temp_wh.TTLBA_shock_trains_first.manipulation.mat', 'TTLBA_shock_trains_first')

load('TTLsKS.mat', 'TO_shock_trains_first') 
TTLTO_shock_trains_first.timestamps = TO_shock_trains_first;
TTLTO_shock_trains_first.timestamps(:,2) = TO_shock_trains_first + 2.5;
save('temp_wh.TTLTO_shock_trains_first.manipulation.mat', 'TTLTO_shock_trains_first')



cell_metrics = ProcessCellMetrics('session', session);

% MD202, MD262-265
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_trains_all') 
TTLBA_trains_all.timestamps = BA_trains_all;
TTLBA_trains_all.timestamps(:,2) = BA_trains_all + 0.1;
save('temp_wh.TTLBA_trains_all.manipulation.mat', 'TTLBA_trains_all')

load('TTLsKS.mat', 'TO_trains_all') 
TTLTO_trains_all.timestamps = TO_trains_all;
TTLTO_trains_all.timestamps(:,2) = TO_trains_all + 0.1;
save('temp_wh.TTLTO_trains_all.manipulation.mat', 'TTLTO_trains_all')

load('TTLsKS.mat', 'BA_trains_first') 
TTLBA_trains_first.timestamps = BA_trains_first;
TTLBA_trains_first.timestamps(:,2) = BA_trains_first + 0.1;
save('temp_wh.TTLBA_trains_first.manipulation.mat', 'TTLBA_trains_first')

load('TTLsKS.mat', 'TO_trains_first') 
TTLTO_trains_first.timestamps = TO_trains_first;
TTLTO_trains_first.timestamps(:,2) = TO_trains_first + 0.1;
save('temp_wh.TTLTO_trains_first.manipulation.mat', 'TTLTO_trains_first')

load('TTLsKS.mat', 'BLUE_only') 
TTLBLUE_only.timestamps = BLUE_only;
TTLBLUE_only.timestamps(:,2) = BLUE_only + 0.1;
save('temp_wh.TTLBLUE_only.manipulation.mat', 'TTLBLUE_only')

cell_metrics = ProcessCellMetrics('session', session);




% MD236_rec001
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_trains_all') 
TTLBA_trains_all.timestamps = BA_trains_all;
TTLBA_trains_all.timestamps(:,2) = BA_trains_all + 0.05;
save('temp_wh.TTLBA_trains_all.manipulation.mat', 'TTLBA_trains_all')

load('TTLsKS.mat', 'BA_trains_first') 
TTLBA_trains_first.timestamps = BA_trains_first;
TTLBA_trains_first.timestamps(:,2) = BA_trains_first + 2.5;
save('temp_wh.TTLBA_trains_first.manipulation.mat', 'TTLBA_trains_first')

load('TTLsKS.mat', 'BAstimON_BeforeShock') 
TTLBAstimON_BeforeShock.timestamps = BAstimON_BeforeShock;
TTLBAstimON_BeforeShock.timestamps(:,2) = BAstimON_BeforeShock + 0.05;
save('temp_wh.TTLBAstimON_BeforeShock.manipulation.mat', 'TTLBAstimON_BeforeShock')

load('TTLsKS.mat', 'BAstimON_BeforeShock_first') 
TTLBAstimON_BeforeShock_first.timestamps = BAstimON_BeforeShock_first;
TTLBAstimON_BeforeShock_first.timestamps(:,2) = BAstimON_BeforeShock_first + 2.5;
save('temp_wh.TTLBAstimON_BeforeShock_first.manipulation.mat', 'TTLBAstimON_BeforeShock_first')

load('TTLsKS.mat', 'BAstimON_AfterShock') 
TTLBAstimON_AfterShock.timestamps = BAstimON_AfterShock;
TTLBAstimON_AfterShock.timestamps(:,2) = BAstimON_AfterShock + 0.05;
save('temp_wh.TTLBAstimON_AfterShock.manipulation.mat', 'TTLBAstimON_AfterShock')

load('TTLsKS.mat', 'BAstimON_AfterShock_first') 
TTLBAstimON_AfterShock_first.timestamps = BAstimON_AfterShock_first;
TTLBAstimON_AfterShock_first.timestamps(:,2) = BAstimON_AfterShock_first + 2.5;
save('temp_wh.TTLBAstimON_AfterShock_first.manipulation.mat', 'TTLBAstimON_AfterShock_first')

cell_metrics = ProcessCellMetrics('session', session);





% MD237
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_trains_all') 
TTLBA_trains_all.timestamps = BA_trains_all;
TTLBA_trains_all.timestamps(:,2) = BA_trains_all + 0.05;
save('temp_wh.TTLBA_trains_all.manipulation.mat', 'TTLBA_trains_all')

load('TTLsKS.mat', 'BA_trains_first') 
TTLBA_trains_first.timestamps = BA_trains_first;
TTLBA_trains_first.timestamps(:,2) = BA_trains_first + 2.5;
save('temp_wh.TTLBA_trains_first.manipulation.mat', 'TTLBA_trains_first')

load('TTLsKS.mat', 'shock_only') 
TTLshock_only.timestamps = shock_only;
TTLshock_only.timestamps(:,2) = shock_only + 2.5;
save('temp_wh.TTLshock_only.manipulation.mat', 'TTLshock_only')

load('TTLsKS.mat', 'shock_BAstim') 
TTLshock_BAstim.timestamps = shock_BAstim;
TTLshock_BAstim.timestamps(:,2) = shock_BAstim + 2.5;
save('temp_wh.TTLshock_BAstim.manipulation.mat', 'TTLshock_BAstim')

load('TTLsKS.mat', 'BA_trains_all') 
TTLBA_trains_all.timestamps = BA_trains_all;
TTLBA_trains_all.timestamps(:,2) = BA_trains_all + 0.05;
save('temp_wh.TTLBA_trains_all.manipulation.mat', 'TTLBA_trains_all')

load('TTLsKS.mat', 'BAstimOFF_AfterShock') 
TTLBAstimOFF_AfterShock.timestamps = BAstimOFF_AfterShock;
TTLBAstimOFF_AfterShock.timestamps(:,2) = BAstimOFF_AfterShock + 0.05;
save('temp_wh.TTLBAstimOFF_AfterShock.manipulation.mat', 'TTLBAstimOFF_AfterShock')

load('TTLsKS.mat', 'BAstimOFF_BeforeShock') 
TTLBAstimOFF_BeforeShock.timestamps = BAstimOFF_BeforeShock;
TTLBAstimOFF_BeforeShock.timestamps(:,2) = BAstimOFF_BeforeShock + 0.05;
save('temp_wh.TTLBAstimOFF_BeforeShock.manipulation.mat', 'TTLBAstimOFF_BeforeShock')

load('TTLsKS.mat', 'BAstimON_AfterShock') 
TTLBAstimON_AfterShock.timestamps = BAstimON_AfterShock;
TTLBAstimON_AfterShock.timestamps(:,2) = BAstimON_AfterShock + 0.05;
save('temp_wh.TTLBAstimON_AfterShock.manipulation.mat', 'TTLBAstimON_AfterShock')

load('TTLsKS.mat', 'BAstimON_BeforeShock') 
TTLBAstimON_BeforeShock.timestamps = BAstimON_BeforeShock;
TTLBAstimON_BeforeShock.timestamps(:,2) = BAstimON_BeforeShock + 0.05;
save('temp_wh.TTLBAstimON_BeforeShock.manipulation.mat', 'TTLBAstimON_BeforeShock')

cell_metrics = ProcessCellMetrics('session', session);


% MD238
session = sessionTemplate(cd);

load('TTLsKS.mat', 'shocks')
shocksTTL.timestamps = shocks;
shocksTTL.timestamps(:,2) = shocks + 0.5;
save('temp_wh.shocksTTL.manipulation.mat', 'shocksTTL')

load('TTLsKS.mat', 'tone_habit_all') 
TTL_tone_habit_all.timestamps = tone_habit_all;
TTL_tone_habit_all.timestamps(:,2) = tone_habit_all + 0.5;
save('temp_wh.TTL_tone_habit_all.manipulation.mat', 'TTL_tone_habit_all')

load('TTLsKS.mat', 'noise_habit_all') 
TTL_noise_habit_all.timestamps = noise_habit_all;
TTL_noise_habit_all.timestamps(:,2) = noise_habit_all + 0.5;
save('temp_wh.TTL_noise_habit_all.manipulation.mat', 'TTL_noise_habit_all')

load('TTLsKS.mat', 'tone_cond_all') 
TTL_tone_cond_all.timestamps = tone_cond_all;
TTL_tone_cond_all.timestamps(:,2) = tone_cond_all + 0.5;
save('temp_wh.TTL_tone_cond_all.manipulation.mat', 'TTL_tone_cond_all')

load('TTLsKS.mat', 'noise_cond_all') 
TTL_noise_cond_all.timestamps = noise_cond_all;
TTL_noise_cond_all.timestamps(:,2) = noise_cond_all + 0.5;
save('temp_wh.TTL_noise_cond_all.manipulation.mat', 'TTL_noise_cond_all')

load('TTLsKS.mat', 'tone_recall_all') 
TTL_tone_recall_all.timestamps = tone_recall_all;
TTL_tone_recall_all.timestamps(:,2) = tone_recall_all + 0.5;
save('temp_wh.TTL_tone_recall_all.manipulation.mat', 'TTL_tone_recall_all')

load('TTLsKS.mat', 'noise_recall_all') 
TTL_noise_recall_all.timestamps = noise_recall_all;
TTL_noise_recall_all.timestamps(:,2) = noise_recall_all + 0.5;
save('temp_wh.TTL_noise_recall_all.manipulation.mat', 'TTL_noise_recall_all')

load('TTLsKS.mat', 'tone_habit_first') 
TTL_tone_habit_first.timestamps = tone_habit_first;
TTL_tone_habit_first.timestamps(:,2) = tone_habit_first + 0.5;
save('temp_wh.TTL_tone_habit_first.manipulation.mat', 'TTL_tone_habit_first')

load('TTLsKS.mat', 'noise_habit_first') 
TTL_noise_habit_first.timestamps = noise_habit_first;
TTL_noise_habit_first.timestamps(:,2) = noise_habit_first + 0.5;
save('temp_wh.TTL_noise_habit_first.manipulation.mat', 'TTL_noise_habit_first')

load('TTLsKS.mat', 'tone_cond_first') 
TTL_tone_cond_first.timestamps = tone_cond_first;
TTL_tone_cond_first.timestamps(:,2) = tone_cond_first + 0.5;
save('temp_wh.TTL_tone_cond_first.manipulation.mat', 'TTL_tone_cond_first')

load('TTLsKS.mat', 'noise_cond_first') 
TTL_noise_cond_first.timestamps = noise_cond_first;
TTL_noise_cond_first.timestamps(:,2) = noise_cond_first + 0.5;
save('temp_wh.TTL_noise_cond_first.manipulation.mat', 'TTL_noise_cond_first')

load('TTLsKS.mat', 'tone_recall_first') 
TTL_tone_recall_first.timestamps = tone_recall_first;
TTL_tone_recall_first.timestamps(:,2) = tone_recall_first + 0.5;
save('temp_wh.TTL_tone_recall_first.manipulation.mat', 'TTL_tone_recall_first')

load('TTLsKS.mat', 'noise_recall_first') 
TTL_noise_recall_first.timestamps = noise_recall_first;
TTL_noise_recall_first.timestamps(:,2) = noise_recall_first + 0.5;
save('temp_wh.TTL_noise_recall_first.manipulation.mat', 'TTL_noise_recall_first')

cell_metrics = ProcessCellMetrics('session', session);

%%PFC BAopto shock
session = sessionTemplate(cd);
load('TTLsKS.mat', 'BA_20Hz_only_first') 
TTLBA_20Hz_only_first.timestamps = BA_20Hz_only_first;
TTLBA_20Hz_only_first.timestamps(:,2) = BA_20Hz_only_first + 2.5;
save('temp_wh.TTLBA_20Hz_only_first.manipulation.mat', 'TTLBA_20Hz_only_first')

load('TTLsKS.mat', 'BA_20Hz_shock_first') 
TTLBA_20Hz_shock_first.timestamps = BA_20Hz_shock_first;
TTLBA_20Hz_shock_first.timestamps(:,2) = BA_20Hz_shock_first + 2.5;
save('temp_wh.TTLBA_20Hz_shock_first.manipulation.mat', 'TTLBA_20Hz_shock_first')

load('TTLsKS.mat', 'BA_20Hz_only_first') 
TTLBA_20Hz_only_first_short.timestamps = BA_20Hz_only_first;
TTLBA_20Hz_only_first_short.timestamps(:,2) = BA_20Hz_only_first + 0.2;
save('temp_wh.TTLBA_20Hz_only_first_short.manipulation.mat', 'TTLBA_20Hz_only_first_short')

load('TTLsKS.mat', 'BA_20Hz_shock_first') 
TTLBA_20Hz_shock_first_short.timestamps = BA_20Hz_shock_first;
TTLBA_20Hz_shock_first_short.timestamps(:,2) = BA_20Hz_shock_first + 0.2;
save('temp_wh.TTLBA_20Hz_shock_first_short.manipulation.mat', 'TTLBA_20Hz_shock_first_short')


load('TTLsKS.mat', 'BA_20Hz_only') 
TTLBA_20Hz_only_all.timestamps = BA_20Hz_only;
TTLBA_20Hz_only_all.timestamps(:,2) = BA_20Hz_only + 0.05;
save('temp_wh.TTLBA_20Hz_only_all.manipulation.mat', 'TTLBA_20Hz_only_all')

load('TTLsKS.mat', 'BA_20Hz_shock') 
TTLBA_20Hz_shock_all.timestamps = BA_20Hz_shock;
TTLBA_20Hz_shock_all.timestamps(:,2) = BA_20Hz_shock + 0.05;
save('temp_wh.TTLBA_20Hz_shock_all.manipulation.mat', 'TTLBA_20Hz_shock_all')

cell_metrics = ProcessCellMetrics('session', session);


% MD279

session = sessionTemplate(cd);

load('TTLsKS.mat', 'shocks')
shocksTTL.timestamps = shocks;
shocksTTL.timestamps(:,2) = shocks + 0.5;
save('temp_wh.shocksTTL.manipulation.mat', 'shocksTTL')

load('TTLsKS.mat', 'sound')
soundTTL.timestamps = sound;
soundTTL.timestamps(:,2) = sound + 0.5;
save('temp_wh.soundTTL.manipulation.mat', 'soundTTL')

cell_metrics = ProcessCellMetrics('session', session);



%MD236-237

basenames = {'temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\PFC_BAopto_SomaticStim\Arhgef6_cre_ChrimsonR\MD236_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_BAopto_SomaticStim\Arhgef6_cre_ChrimsonR\MD237_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);

% MD324
basenames = {'temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_shock_run_rest\MD233_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_shock_run_rest\MD234_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_shock_run_rest\MD235_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);




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

%% BA run-rest oep 6.3
basenames = {'temp_wh','temp_wh','temp_wh', 'temp_wh'};
basepaths = {...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\AMY_shok_rest_run\MD171_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\AMY_shok_rest_run\MD192_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\AMY_shok_rest_run\MD193_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\AMY_shok_rest_run\MD194_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);

%%PAG shock run rest OEPv063
basenames = {'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\PAG_shock_rest_run\MD195_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PAG_shock_rest_run\MD196_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PAG_shock_rest_run\MD197_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PAG_shock_rest_run\MD198_kilosort\kilosort25preprocess'};
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
basenames = {'temp_wh','temp_wh','temp_wh','temp_wh','temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'Z:\HajosLab\Dani\M2_shock_response\MD127_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\M2_shock_response\MD128_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\M2_shock_response\MD129_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\M2_shock_response\MD181_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\M2_shock_response\MD182_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\M2_shock_response\MD183_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\M2_shock_response\MD187_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\M2_shock_response\MD188_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);

basenames = {'temp_wh','temp_wh','temp_wh','temp_wh','temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'Z:\HajosLab\Dani\Magyar_Daniel\recordings\MD144_20221115_001\MD144_kilosort_waveform', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\recordings\MD145_20221116_001\MD145_kilosort_waveform', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\recordings\MD146_20221117_001\MD146_kilosort_waveform', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\recordings\MD147_20221118_001\MD147_kilosort_waveform', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\recordings\MD149_20221123_001\MD149_kilosort_waveform', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\recordings\MD150_20221202_001\MD150_kilosort_waveform', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\recordings\MD151_20221205_001\MD151_kilosort_waveform', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\recordings\MD152_20221207_001\MD152_kilosort_waveform'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);

% BA-fear cond
basenames = {'temp_wh','temp_wh','temp_wh', 'temp_wh', 'temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD243_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD250_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD252_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD253_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD254_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

basenames = {'temp_wh','temp_wh','temp_wh', 'temp_wh', 'temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD266_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD267_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD268_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD269_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

basenames = {'temp_wh','temp_wh','temp_wh', 'temp_wh', 'temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD275_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD276_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD277_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD278_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);




cell_metrics.general.TTL_shocks = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks;
    cell_metrics.general.TTL_shocks(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics.general.TTL_tone_habit_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_habit_first;
    cell_metrics.general.TTL_tone_habit_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_habit_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_habit_first;
    cell_metrics.general.TTL_noise_habit_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_tone_cond_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_cond_first;
    cell_metrics.general.TTL_tone_cond_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_cond_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_cond_first;
    cell_metrics.general.TTL_noise_cond_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_tone_recall_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_recall_first;
    cell_metrics.general.TTL_tone_recall_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_recall_first = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_recall_first;
    cell_metrics.general.TTL_noise_recall_first(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_tone_habit_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_habit_all;
    cell_metrics.general.TTL_tone_habit_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_habit_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_habit_all;
    cell_metrics.general.TTL_noise_habit_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_tone_cond_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_cond_all;
    cell_metrics.general.TTL_tone_cond_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_cond_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_cond_all;
    cell_metrics.general.TTL_noise_cond_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_tone_recall_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.tone_recall_all;
    cell_metrics.general.TTL_tone_recall_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_noise_recall_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.noise_recall_all;
    cell_metrics.general.TTL_noise_recall_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics = CellExplorer('metrics',cell_metrics);



% M2 Arch
basenames = {'temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\M2_Arch\MD255_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_Arch\MD256_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_Arch\MD257_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

cell_metrics.general.TTL_shocks_only = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks_only;
    cell_metrics.general.TTL_shocks_only(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_shocks_light = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks_light;
    cell_metrics.general.TTL_shocks_light(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL


cell_metrics = CellExplorer('metrics',cell_metrics);


% M2 SST Jaws
basenames = {'temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\M2_SST_Jaws\MD270_001_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_SST_Jaws\MD270_002_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_SST_Jaws\MD271_002_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\M2_SST_Jaws\MD271_003_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

cell_metrics.general.TTL_shocks_only = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks_only;
    cell_metrics.general.TTL_shocks_only(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_shocks_light = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks_light;
    cell_metrics.general.TTL_shocks_light(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL

cell_metrics = CellExplorer('metrics',cell_metrics);

%M2 optotag Charmine
session = sessionTemplate(cd);

load('TTLsKS.mat', 'shocks')
shocksTTL.timestamps = shocks;
shocksTTL.timestamps(:,2) = shocks + 0.5;
save('temp_wh.shocksTTL.manipulation.mat', 'shocksTTL')
load('TTLsKS.mat', 'pulses')
pulsesTTL.timestamps = pulses;
pulsesTTL.timestamps(:,2) = pulses + 0.5;
save('temp_wh.pulsesTTL.manipulation.mat', 'pulsesTTL')
cell_metrics = ProcessCellMetrics('session', session);



basenames = {'temp_wh','temp_wh','temp_wh','temp_wh','temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_SST_Chrmine\MD258_002_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_SST_Chrmine\MD258_003_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_SST_Chrmine\MD259_001_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_SST_Chrmine\MD259_004_kilosort\kilosort25preprocess',...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_SST_Chrmine\MD260_001_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_SST_Chrmine\MD260_002_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_SST_Chrmine\MD261_001_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_SST_Chrmine\MD261_003_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

cell_metrics.general.TTL_shocks = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks;
    cell_metrics.general.TTL_shocks(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_pulses = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.pulses;
    cell_metrics.general.TTL_pulses(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL



cell_metrics = CellExplorer('metrics',cell_metrics);


% MD248, 262-265
basenames = {'temp_wh', 'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\PFC_BAopto_SomaticStim\Chrmine_stGtACR2\MD248_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_BAopto_SomaticStim\Chrmine_stGtACR2\MD262_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_BAopto_SomaticStim\Chrmine_stGtACR2\MD263_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_BAopto_SomaticStim\Chrmine_stGtACR2\MD264_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\PFC_BAopto_SomaticStim\Chrmine_stGtACR2\MD265_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

cell_metrics.general.TTL_BA_trains_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.BA_trains_all;
    cell_metrics.general.TTL_BA_trains_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_TO_trains_all = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.TO_trains_all;
    cell_metrics.general.TTL_TO_trains_all(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics = CellExplorer('metrics',cell_metrics);



% MD281,MD282, MD2823
basenames = {'temp_wh', 'temp_wh','temp_wh','temp_wh','temp_wh','temp_wh', 'temp_wh','temp_wh','temp_wh','temp_wh', 'temp_wh'};
basepaths = {...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_NDNF_VIP2R_Arch\MD281_001_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_NDNF_VIP2R_Arch\MD281_002_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_NDNF_VIP2R_Arch\MD281_003_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_NDNF_VIP2R_Arch\MD281_004_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_NDNF_VIP2R_Arch\MD281_005_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_NDNF_VIP2R_Arch\MD282_001_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_NDNF_VIP2R_Arch\MD282_002_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_NDNF_VIP2R_Arch\MD283_002_kilosort\kilosort25preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\M2_NDNF_VIP2R_Arch\MD283_003_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);



cell_metrics.general.TTL_shocks_only = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks_only;
    cell_metrics.general.TTL_shocks_only(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics.general.TTL_shocks_light = {};
for ii = 1:max(cell_metrics.batchIDs)
    allTTL = load([cell_metrics.general.basepaths{ii} '\TTLsKS.mat']);
    TTL = allTTL.shocks_light;
    cell_metrics.general.TTL_shocks_light(end+1:end+sum(cell_metrics.batchIDs == ii)) = {TTL};
end
clear TTL AP ii allTTL
cell_metrics = CellExplorer('metrics',cell_metrics);





%MD288, 289, 290, 291


session = sessionTemplate(cd);

load('TTLsKS.mat', 'shocks')
shocksTTL.timestamps = shocks;
shocksTTL.timestamps(:,2) = shocks + 0.5;
save('temp_wh.shocksTTL.manipulation.mat', 'shocksTTL')

load('TTLsKS.mat', 'tone_habit_all') 
TTL_tone_habit_all.timestamps = tone_habit_all;
TTL_tone_habit_all.timestamps(:,2) = tone_habit_all + 0.5;
save('temp_wh.TTL_tone_habit_all.manipulation.mat', 'TTL_tone_habit_all')

load('TTLsKS.mat', 'tone_cond_all') 
TTL_tone_cond_all.timestamps = tone_cond_all;
TTL_tone_cond_all.timestamps(:,2) = tone_cond_all + 0.5;
save('temp_wh.TTL_tone_cond_all.manipulation.mat', 'TTL_tone_cond_all')

load('TTLsKS.mat', 'tone_recall_all') 
TTL_tone_recall_all.timestamps = tone_recall_all;
TTL_tone_recall_all.timestamps(:,2) = tone_recall_all + 0.5;
save('temp_wh.TTL_tone_recall_all.manipulation.mat', 'TTL_tone_recall_all')

load('TTLsKS.mat', 'tone_habit_first') 
TTL_tone_habit_first.timestamps = tone_habit_first;
TTL_tone_habit_first.timestamps(:,2) = tone_habit_first + 0.5;
save('temp_wh.TTL_tone_habit_first.manipulation.mat', 'TTL_tone_habit_first')

load('TTLsKS.mat', 'tone_cond_first') 
TTL_tone_cond_first.timestamps = tone_cond_first;
TTL_tone_cond_first.timestamps(:,2) = tone_cond_first + 0.5;
save('temp_wh.TTL_tone_cond_first.manipulation.mat', 'TTL_tone_cond_first')

load('TTLsKS.mat', 'tone_recall_first') 
TTL_tone_recall_first.timestamps = tone_recall_first;
TTL_tone_recall_first.timestamps(:,2) = tone_recall_first + 0.5;
save('temp_wh.TTL_tone_recall_first.manipulation.mat', 'TTL_tone_recall_first')

cell_metrics = ProcessCellMetrics('session', session);


basenames = {'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD288_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD289_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD290_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD291_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

cell_metrics = CellExplorer('metrics',cell_metrics);





% Tibi

session = sessionTemplate(cd);

load('TTLsKS.mat', 'shocks')
shocksTTL.timestamps = shocks;
shocksTTL.timestamps(:,2) = shocks + 0.5;
save('temp_wh.shocksTTL.manipulation.mat', 'shocksTTL')

load('TTLsKS.mat', 'pulses')
pulsesTTL.timestamps = pulses;
pulsesTTL.timestamps(:,2) = pulses + 0.5;
save('temp_wh.pulsesTTL.manipulation.mat', 'pulsesTTL')

load('TTLsKS.mat', 'firstlickUniq')
firstlickUniqTTL.timestamps = firstlickUniq;
firstlickUniqTTL.timestamps(:,2) = firstlickUniq + 0.5;
save('temp_wh.firstlickUniqTTL.manipulation.mat', 'firstlickUniqTTL')

load('TTLsKS.mat', 'nonRewardlick')
nonRewardlickTTL.timestamps = nonRewardlick;
nonRewardlickTTL.timestamps(:,2) = nonRewardlick + 0.5;
save('temp_wh.nonRewardlickTTL.manipulation.mat', 'nonRewardlickTTL')

cell_metrics = ProcessCellMetrics('session', session);






%MD294


session = sessionTemplate(cd);

load('TTLsKS.mat', 'shocks')
shocksTTL.timestamps = shocks;
shocksTTL.timestamps(:,2) = shocks + 0.5;
save('temp_wh.shocksTTL.manipulation.mat', 'shocksTTL')

load('TTLsKS.mat', 'shocks_nonpredicted')
shocks_nonpredictedTTL.timestamps = shocks_nonpredicted;
shocks_nonpredictedTTL.timestamps(:,2) = shocks_nonpredicted + 0.5;
save('temp_wh.shocks_nonpredictedTTL.manipulation.mat', 'shocks_nonpredictedTTL')

load('TTLsKS.mat', 'shocks_predicted')
shocks_predictedTTL.timestamps = shocks_predicted;
shocks_predictedTTL.timestamps(:,2) = shocks_predicted + 0.5;
save('temp_wh.shocks_predictedTTL.manipulation.mat', 'shocks_predictedTTL')

load('TTLsKS.mat', 'tone_habit_all') 
TTL_tone_habit_all.timestamps = tone_habit_all;
TTL_tone_habit_all.timestamps(:,2) = tone_habit_all + 0.5;
save('temp_wh.TTL_tone_habit_all.manipulation.mat', 'TTL_tone_habit_all')

load('TTLsKS.mat', 'tone_cond_all') 
TTL_tone_cond_all.timestamps = tone_cond_all;
TTL_tone_cond_all.timestamps(:,2) = tone_cond_all + 0.5;
save('temp_wh.TTL_tone_cond_all.manipulation.mat', 'TTL_tone_cond_all')

load('TTLsKS.mat', 'tone_recall_all') 
TTL_tone_recall_all.timestamps = tone_recall_all;
TTL_tone_recall_all.timestamps(:,2) = tone_recall_all + 0.5;
save('temp_wh.TTL_tone_recall_all.manipulation.mat', 'TTL_tone_recall_all')

load('TTLsKS.mat', 'tone_habit_first') 
TTL_tone_habit_first.timestamps = tone_habit_first;
TTL_tone_habit_first.timestamps(:,2) = tone_habit_first + 0.5;
save('temp_wh.TTL_tone_habit_first.manipulation.mat', 'TTL_tone_habit_first')

load('TTLsKS.mat', 'tone_cond_first') 
TTL_tone_cond_first.timestamps = tone_cond_first;
TTL_tone_cond_first.timestamps(:,2) = tone_cond_first + 0.5;
save('temp_wh.TTL_tone_cond_first.manipulation.mat', 'TTL_tone_cond_first')

load('TTLsKS.mat', 'tone_recall_first') 
TTL_tone_recall_first.timestamps = tone_recall_first;
TTL_tone_recall_first.timestamps(:,2) = tone_recall_first + 0.5;
save('temp_wh.TTL_tone_recall_first.manipulation.mat', 'TTL_tone_recall_first')

load('TTLsKS.mat', 'triptest_shocks_only')
triptest_shocks_onlyTTL.timestamps = triptest_shocks_only;
triptest_shocks_onlyTTL.timestamps(:,2) = triptest_shocks_only + 0.5;
save('temp_wh.triptest_shocks_onlyTTL.manipulation.mat', 'triptest_shocks_onlyTTL')

load('TTLsKS.mat', 'triptest_sound_only')
triptest_sound_onlyTTL.timestamps = triptest_sound_only;
triptest_sound_onlyTTL.timestamps(:,2) = triptest_sound_only + 0.5;
save('temp_wh.triptest_sound_onlyTTL.manipulation.mat', 'triptest_sound_onlyTTL')

load('TTLsKS.mat', 'triptest_both')
triptest_bothTTL.timestamps = triptest_both;
triptest_bothTTL.timestamps(:,2) = triptest_both + 0.5;
save('temp_wh.triptest_bothTTL.manipulation.mat', 'triptest_bothTTL')

cell_metrics = ProcessCellMetrics('session', session);




basenames = {'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh', 'temp_wh'};
basepaths = {...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD293_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD294_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD295_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD296_kilosort\kilosort25preprocess', ...
    'C:\Users\dmagyar\Desktop\BA_fear_cond\MD297_kilosort\kilosort25preprocess'};
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

cell_metrics = CellExplorer('metrics',cell_metrics);

