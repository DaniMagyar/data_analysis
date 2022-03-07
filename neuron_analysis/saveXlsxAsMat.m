clear all;
name = 'PFC_Chrimson_stGtACR_2021_dec';
filename = [name '.xlsx'];
opts = detectImportOptions(name);
opts = setvartype(opts, 'Group', 'char');
PFC_Chrimson_stGtACR_2021_dec = readtable(name, opts);
save ([name '.mat'], 'PFC_Chrimson_stGtACR_2021_dec')
clear all;