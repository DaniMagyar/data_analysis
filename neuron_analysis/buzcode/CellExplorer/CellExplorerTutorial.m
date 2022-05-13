session = sessionTemplate(cd);



load('TTLsKS.mat', 'BA_25_5Hz')
laserBA_25_5Hz_short.timestamps = BA_25_5Hz;
laserBA_25_5Hz_short.timestamps(:,2) = BA_25_5Hz+0.1;
save('temp_wh.laserBA_25_5Hz_short.manipulation.mat', 'laserBA_25_5Hz_short')
cell_metrics = ProcessCellMetrics('session', session);

load('TTLsKS.mat', 'TO_25_5Hz')
laserTO_25_5Hz_short.timestamps = TO_25_5Hz;
laserTO_25_5Hz_short.timestamps(:,2) = TO_25_5Hz+0.1;
save('temp_wh.laserTO_25_5Hz_short.manipulation.mat', 'laserTO_25_5Hz_short')
cell_metrics = ProcessCellMetrics('session', session);

load('TTLsKS.mat', 'BA_250_5Hz')
laserBA_250_5Hz.timestamps = BA_250_5Hz;
save('temp_wh.laserBA_250_5Hz.manipulation.mat', 'laserBA_250_5Hz')
cell_metrics = ProcessCellMetrics('session', session);











basenames = {'temp_wh','temp_wh','temp_wh','temp_wh'};
basepaths = {...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december\Kilosort_v2\MD098_kilosort\kilosort3preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december\Kilosort_v2\MD099_kilosort\kilosort3preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december\Kilosort_v2\MD101_kilosort\kilosort3preprocess', ...
    'Z:\HajosLab\Dani\Magyar_Daniel\experiments\PFC_layers\Chrimson_stGtACR\2021_december\Kilosort_v2\MD108_kilosort\kilosort3preprocess'};
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
