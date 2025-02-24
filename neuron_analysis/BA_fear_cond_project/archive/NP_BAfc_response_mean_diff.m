function [respmean_diff] = NP_BAfc_response_mean_diff


int = [-10 15];
Wcx_win = [-5 5];
psth_bin = 6000;
fs = 30000;

[PSTHall_hab] = Zscore_analysis_CellExp_v3('mainFolder','C:\Users\dmagyar\Desktop\BA_fear_cond', ...
    'Stim','tone_habit_first', 'Recordings', {'MD288', 'MD289', 'MD290', 'MD291'},...
    'psth_bin',psth_bin,'int', int, 'Wcx_win', Wcx_win, 'plotType', 'none', 'smooth', 1, 'smoothvalue', 5, 'Sorted', 0);


[PSTHall_rec] = Zscore_analysis_CellExp_v3('mainFolder','C:\Users\dmagyar\Desktop\BA_fear_cond', ...
    'Stim','tone_recall_first', 'Recordings', {'MD288', 'MD289', 'MD290', 'MD291'},...
    'psth_bin',psth_bin,'int', int, 'Wcx_win', Wcx_win, 'plotType', 'none', 'smooth', 1, 'smoothvalue', 5, 'Sorted', 0);


respmean_habit = mean(PSTHall_hab(:,abs(int(1))*fs/psth_bin+1:(abs(int(1))+abs(Wcx_win(2)))*fs/psth_bin),2);

respmean_recall = mean(PSTHall_rec(:,abs(int(1))*fs/psth_bin+1:(abs(int(1))+abs(Wcx_win(2)))*fs/psth_bin),2);

respmean_diff = respmean_recall - respmean_habit;