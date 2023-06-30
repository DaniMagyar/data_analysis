function [preferences] = preferences_PsthExplorer
% Default params ----------------------------------------------------------
preferences.fs       = 30000; % Firing rate
preferences.int      = [-5 5]; % psth interval, even numbers are better for plot
preferences.norm     = 1; % Normalise data
preferences.psth_bin = 6000; % 600 = 20ms
preferences.average = 1; % plot averages, 1 if yes
preferences.mainFolder = 'C:\Users\dmagyar\Desktop\2022_august\kilosort';
disp(preferences.mainFolder)
preferences.plottingMethod = 'PCA'; % 'PCA' or 'Zscore'

% Zscore analysis settings
preferences.Wcx_win = [-5 2]; % Wilcoxon window in seconds
preferences.testBins = 10; % bins included for Z-score ordering 

% PCA analysis setting
preferences.PCA_num  = 2; % Number of principal components to use
preferences.clustnum = 5; % Number of clusters
preferences.LKmethod = 'complete';
preferences.LKmetric = 'mahalanobis';

% Time axis setup
preferences.time = (preferences.int(1)*preferences.fs:preferences.psth_bin:preferences.int(2)*preferences.fs)';
preferences.time = preferences.time(1:end-1);

% Colormap blue:white:red -------------------------------------------------
c1 = 1/255*[0,115,185];
c2 = 1/255*[239,239,239];
c3 = 1/255*[217,84,26];
preferences.mycolormap(1:20,1)  = linspace(c1(1),c2(1),20);
preferences.mycolormap(21:64,1) = linspace(c2(1),c3(1),44);
preferences.mycolormap(1:20,2)  = linspace(c1(2),c2(2),20);
preferences.mycolormap(21:64,2) = linspace(c2(2),c3(2),44);
preferences.mycolormap(1:20,3)  = linspace(c1(3),c2(3),20);
preferences.mycolormap(21:64,3) = linspace(c2(3),c3(3),44);